#!/usr/bin/env python3
"""Resolve BAM inputs for the unified low-pass per-sample variant workflow."""
import argparse
import csv
import os
import re
import sys
from pathlib import Path


def sanitize_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value).replace(" ", "_"))


def sample_from_bam(path: str) -> str:
    b = os.path.basename(path)
    if b.endswith(".bam"):
        b = b[:-4]
    for suffix in [".sorted", ".markdup", ".dedup", ".bqsr", ".preprocessed"]:
        if b.endswith(suffix):
            b = b[: -len(suffix)]
    return sanitize_id(b)


def norm_path(path: str) -> str:
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))


def parse_manifest(path: str):
    out = []
    with open(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if "\t" in line:
                fields = line.split("\t")
            elif "," in line:
                fields = line.split(",")
            else:
                fields = [line]
            for field in fields:
                field = field.strip().strip('"').strip("'")
                if field.endswith(".bam"):
                    out.append(norm_path(field))
                    break
    return out


def resolve_bams(input_value: str, recursive: bool):
    input_value = input_value.strip()
    if os.path.isdir(input_value):
        root = Path(input_value)
        pattern = "**/*.bam" if recursive else "*.bam"
        paths = [norm_path(str(p)) for p in sorted(root.glob(pattern)) if p.is_file() or p.is_symlink()]
    elif os.path.isfile(input_value) and not input_value.endswith(".bam"):
        paths = parse_manifest(input_value)
    elif "," in input_value:
        paths = [norm_path(x.strip()) for x in input_value.split(",") if x.strip()]
    else:
        paths = [norm_path(input_value)]

    unique, seen = [], set()
    for p in paths:
        if p not in seen:
            unique.append(p)
            seen.add(p)

    missing = [p for p in unique if not (os.path.exists(p) and os.path.getsize(p) > 0)]
    if missing:
        raise SystemExit("ERROR: missing or empty BAM(s):\n" + "\n".join(missing))
    if not unique:
        raise SystemExit(f"ERROR: no BAMs resolved from input: {input_value}")
    return unique


def selected_from_metadata(metadata_tsv: str, cohort: str):
    if cohort == "all":
        return None
    if not metadata_tsv:
        raise SystemExit("ERROR: --database-tsv is required when --cohort is lymphoma or brain")
    if not os.path.exists(metadata_tsv):
        raise SystemExit(f"ERROR: metadata TSV not found: {metadata_tsv}")

    lymphoma_re = re.compile(
        r"lymphoma|hodgkin|dlbcl|diffuse large|follicular|mantle|small lymphocytic|"
        r"chronic lymphocytic|\bcll\b|high-grade b-cell|high grade b-cell|"
        r"b-cell non-hodgkin|t-cell lymphoma|sézary|sezary", re.I)
    brain_re = re.compile(
        r"brain|enceph|central nervous|\bcns\b|glioma|glioblastoma|meningioma|"
        r"arachnoid|cerebral|calvarium|dura|dural|metastatic.*brain|brain metast", re.I)
    regex = lymphoma_re if cohort == "lymphoma" else brain_re
    selected = set()
    with open(metadata_tsv, newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sid = (row.get("illumina_sample_id") or "").strip()
            if not sid:
                continue
            seq = (row.get("illumina_sequenced") or "").strip().lower()
            if seq and seq not in {"yes", "true", "1", "y", "si", "sí"}:
                continue
            text = " | ".join(str(row.get(k, "")) for k in [
                "specimen_organ", "anatomical_site_1", "anatomical_site_2",
                "clinical_diagnosis_1", "clinical_diagnosis_2",
                "diagnosis_category_1", "diagnosis_category_2", "final_diagnosis",
                "microscopic_summary_en"])
            if regex.search(text):
                selected.add(sid)
    return selected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--samples", default="")
    ap.add_argument("--recursive", action="store_true")
    ap.add_argument("--database-tsv", default="")
    ap.add_argument("--cohort", default="all", choices=["all", "lymphoma", "brain"])
    args = ap.parse_args()

    bams = resolve_bams(args.input, args.recursive)
    allowed = selected_from_metadata(args.database_tsv, args.cohort)
    if args.samples.strip():
        sample_set = {x.strip() for x in args.samples.replace(";", ",").split(",") if x.strip()}
        allowed = sample_set if allowed is None else (allowed & sample_set)

    rows = []
    for bam in bams:
        sample = sample_from_bam(bam)
        if allowed is None or sample in allowed:
            rows.append((sample, bam))
    if not rows:
        raise SystemExit("ERROR: no BAMs remain after --samples/--cohort filtering")

    by_sample = {}
    for sample, bam in rows:
        by_sample.setdefault(sample, []).append(bam)
    duplicates = {sample: paths for sample, paths in by_sample.items() if len(paths) > 1}
    if duplicates:
        details = []
        for sample in sorted(duplicates):
            details.append(f"  {sample}: " + ", ".join(duplicates[sample]))
        raise SystemExit(
            "ERROR: multiple BAMs resolve to the same logical sample ID; "
            "rename inputs or select one BAM per sample:\n" + "\n".join(details)
        )

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "bam"])
        writer.writerows(rows)
    sys.stderr.write(f"Resolved {len(rows)} selected BAM(s).\n")


if __name__ == "__main__":
    main()
