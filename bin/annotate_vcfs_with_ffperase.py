#!/usr/bin/env python3
"""
Annotate final per-sample VCFs with FFPErase classification tables.

Version: 2026-06-01-v2-safe-samefile-copy

Creates a clean folder, by default:
  <root>/ffperase_classification/final_vcf/*.final_variants.ffperase_annotated.vcf.gz
  <root>/ffperase_classification/reports/all_samples.ffperase_annotation_summary.tsv

This script does not re-run FFPErase. It only joins the existing FFPErase TSV
classification output back into your pipeline's final_vcf files by
CHROM/POS/REF/ALT.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

TAB = "\t"

INFO_HEADERS = [
    '##INFO=<ID=FFPERASE_EVALUATED,Number=1,Type=Integer,Description="1 if this allele was found in the FFPErase classification table, 0 otherwise">',
    '##INFO=<ID=FFPERASE_DECISION,Number=1,Type=String,Description="FFPErase simplified decision: REAL, ARTIFACT, or NOT_EVALUATED. REAL means not classified as FFPE artifact; it does not prove somatic origin">',
    '##INFO=<ID=FFPERASE_MUTATION_TYPE,Number=1,Type=String,Description="FFPErase mutation type used for classification: snvs or indels">',
    '##INFO=<ID=FFPERASE_ARTIFACT,Number=1,Type=Integer,Description="1 if FFPErase predicts likely FFPE artifact, 0 if not">',
    '##INFO=<ID=FFPERASE_PROB,Number=1,Type=Float,Description="FFPErase raw predicted probability/score for artifact class">',
    '##INFO=<ID=FFPERASE_VAF,Number=1,Type=Float,Description="Variant allele fraction reported by FFPErase pileup feature extraction">',
    '##INFO=<ID=FFPERASE_DEPTH,Number=1,Type=Integer,Description="Depth reported by FFPErase pileup feature extraction">',
    '##INFO=<ID=FFPERASE_VARIANT_READS,Number=1,Type=Integer,Description="Variant-supporting reads reported by FFPErase">',
    '##INFO=<ID=FFPERASE_STRAND_BIAS,Number=1,Type=Float,Description="Strand-bias feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_AVG_BQ,Number=1,Type=Float,Description="Average base quality feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_AVG_ALT_BQ,Number=1,Type=Float,Description="Average alternate-read base quality feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_AVG_MQ,Number=1,Type=Float,Description="Average mapping quality feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_AVG_ALT_MQ,Number=1,Type=Float,Description="Average alternate-read mapping quality feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_LOG_IS_RATIO,Number=1,Type=Float,Description="Log insert-size ratio feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_LOG_ALT_IS_RATIO,Number=1,Type=Float,Description="Log alternate-read insert-size ratio feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_AVG_EDIT_DIST,Number=1,Type=Float,Description="Average edit-distance feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_AVG_READ_BAL,Number=1,Type=Float,Description="Average read-balance feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_LOG_DEPTH_RATIO,Number=1,Type=Float,Description="Log depth-ratio feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_PA_BASE_CHANGE_ERROR,Number=1,Type=Float,Description="Picard pre-adapter base-change error feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_PA_TRINUCLEO_ERROR,Number=1,Type=Float,Description="Picard pre-adapter trinucleotide error feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_BB_BASE_CHANGE_ERROR,Number=1,Type=Float,Description="Picard bait-bias base-change error feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_BB_TRINUCLEO_ERROR,Number=1,Type=Float,Description="Picard bait-bias trinucleotide error feature reported by FFPErase">',
    '##INFO=<ID=FFPERASE_INDEL_TYPE,Number=1,Type=String,Description="FFPErase indel type: I insertion or D deletion">',
    '##INFO=<ID=FFPERASE_INDEL_LENGTH,Number=1,Type=Integer,Description="Indel length reported by FFPErase">',
    '##INFO=<ID=FFPERASE_MHCOUNT,Number=1,Type=Integer,Description="Microhomology count reported by FFPErase for indels">',
    '##INFO=<ID=FFPERASE_REPCOUNT,Number=1,Type=Integer,Description="Repeat count reported by FFPErase for indels">',
    '##INFO=<ID=FFPERASE_INDEL_CLASSIFICATION,Number=1,Type=String,Description="Indel sequence-context classification reported by FFPErase, e.g. Repeat-mediated or None">',
    '##INFO=<ID=FFPERASE_INDEL_COUNT,Number=1,Type=Float,Description="Indel count feature reported by FFPErase">',
]

FILTER_HEADER = '##FILTER=<ID=FFPERASE_ARTIFACT,Description="Variant classified as likely FFPE artifact by FFPErase">'

NUMERIC_FIELDS = {
    "FFPERASE_PROB": "float",
    "FFPERASE_VAF": "float",
    "FFPERASE_DEPTH": "int",
    "FFPERASE_VARIANT_READS": "int",
    "FFPERASE_STRAND_BIAS": "float",
    "FFPERASE_AVG_BQ": "float",
    "FFPERASE_AVG_ALT_BQ": "float",
    "FFPERASE_AVG_MQ": "float",
    "FFPERASE_AVG_ALT_MQ": "float",
    "FFPERASE_LOG_IS_RATIO": "float",
    "FFPERASE_LOG_ALT_IS_RATIO": "float",
    "FFPERASE_AVG_EDIT_DIST": "float",
    "FFPERASE_AVG_READ_BAL": "float",
    "FFPERASE_LOG_DEPTH_RATIO": "float",
    "FFPERASE_PA_BASE_CHANGE_ERROR": "float",
    "FFPERASE_PA_TRINUCLEO_ERROR": "float",
    "FFPERASE_BB_BASE_CHANGE_ERROR": "float",
    "FFPERASE_BB_TRINUCLEO_ERROR": "float",
    "FFPERASE_INDEL_LENGTH": "int",
    "FFPERASE_MHCOUNT": "int",
    "FFPERASE_REPCOUNT": "int",
    "FFPERASE_INDEL_COUNT": "float",
}

FIELD_MAP = {
    "VAF": "FFPERASE_VAF",
    "STRAND_BIAS": "FFPERASE_STRAND_BIAS",
    "AVG_BQ": "FFPERASE_AVG_BQ",
    "AVG_ALT_BQ": "FFPERASE_AVG_ALT_BQ",
    "AVG_MQ": "FFPERASE_AVG_MQ",
    "AVG_ALT_MQ": "FFPERASE_AVG_ALT_MQ",
    "LOG_IS_RATIO": "FFPERASE_LOG_IS_RATIO",
    "LOG_ALT_IS_RATIO": "FFPERASE_LOG_ALT_IS_RATIO",
    "AVG_EDIT_DIST": "FFPERASE_AVG_EDIT_DIST",
    "AVG_READ_BAL": "FFPERASE_AVG_READ_BAL",
    "DEPTH": "FFPERASE_DEPTH",
    "LOG_DEPTH_RATIO": "FFPERASE_LOG_DEPTH_RATIO",
    "VARIANT_READS": "FFPERASE_VARIANT_READS",
    "PA_BASE_CHANGE_ERROR": "FFPERASE_PA_BASE_CHANGE_ERROR",
    "PA_TRINUCLEO_ERROR": "FFPERASE_PA_TRINUCLEO_ERROR",
    "BB_BASE_CHANGE_ERROR": "FFPERASE_BB_BASE_CHANGE_ERROR",
    "BB_TRINUCLEO_ERROR": "FFPERASE_BB_TRINUCLEO_ERROR",
    "INDEL_LENGTH": "FFPERASE_INDEL_LENGTH",
    "INDEL_TYPE": "FFPERASE_INDEL_TYPE",
    "MHCOUNT": "FFPERASE_MHCOUNT",
    "REPCOUNT": "FFPERASE_REPCOUNT",
    "CLASSIFICATION": "FFPERASE_INDEL_CLASSIFICATION",
    "INDEL_COUNT": "FFPERASE_INDEL_COUNT",
}

RAW_PREDICT_COLS = ["ARTIFACT_raw_predicts", "FFPErase_raw_predicts", "FFPErase_snvs_raw_predicts", "FFPErase_indels_raw_predicts"]
PREDICT_COLS = ["ARTIFACT_predicts", "FFPErase_predicts", "FFPErase_snvs_predicts", "FFPErase_indels_predicts"]


@dataclass
class VcfJob:
    sample: str
    vcf: Path


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def norm_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    c = c.upper()
    if c == "M":
        c = "MT"
    return c


def bool_to_int(value: object) -> Optional[int]:
    s = str(value).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return 1
    if s in {"false", "f", "0", "no", "n"}:
        return 0
    return None


def is_missing(value: object) -> bool:
    if value is None:
        return True
    s = str(value).strip()
    return s == "" or s == "." or s.lower() in {"nan", "none", "null"}


def clean_numeric(value: object, kind: str) -> Optional[str]:
    if is_missing(value):
        return None
    s = str(value).strip()
    try:
        if kind == "int":
            return str(int(float(s)))
        return f"{float(s):.6g}"
    except Exception:
        return None


def clean_string(value: object) -> Optional[str]:
    if is_missing(value):
        return None
    s = str(value).strip()
    # Keep VCF INFO values safe without over-encoding everything.
    s = s.replace("%", "%25")
    s = s.replace(";", "%3B").replace("=", "%3D").replace(",", "%2C")
    s = re.sub(r"\s+", "_", s)
    return s or None


def split_header_line(line: str) -> List[str]:
    line = line.rstrip("\n\r")
    line = line.replace("BB_TRINUCLEO_ERRORARTIFACT_raw_predicts", "BB_TRINUCLEO_ERROR\tARTIFACT_raw_predicts")
    line = line.replace("BB_TRINUCLEO_ERRORFFPErase", "BB_TRINUCLEO_ERROR\tFFPErase")
    if "\t" in line:
        return line.split("\t")
    return re.split(r"\s+", line.strip())


def split_data_line(line: str, header_len: int) -> List[str]:
    raw = line.rstrip("\n\r")
    if not raw:
        return []
    if "\t" in raw:
        f = raw.split("\t")
    else:
        f = re.split(r"\s+", raw.strip())
    # Repair a common display/copy artifact where the final probability and boolean
    # are pasted as 0.1769False with no delimiter. This should rarely be needed
    # for real TSV files, but it makes the parser tolerant of pasted outputs.
    if len(f) == header_len - 1 and f:
        m = re.match(r"^(.+?)(True|False|true|false)$", f[-1])
        if m:
            f = f[:-1] + [m.group(1), m.group(2)]
    return f


def read_tsv_rows(path: Path, sample_hint: Optional[str], mutation_type_hint: Optional[str]) -> List[dict]:
    rows: List[dict] = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with open_text(path) as fh:
        header: Optional[List[str]] = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            header = split_header_line(line)
            break
        if not header:
            return rows
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = split_data_line(line, len(header))
            if not f:
                continue
            if len(f) < len(header):
                f += [""] * (len(header) - len(f))
            elif len(f) > len(header):
                # If a pasted table has too many fields, keep the leading columns and
                # merge extras into the final column rather than crashing.
                f = f[: len(header) - 1] + ["_".join(f[len(header) - 1 :])]
            d = dict(zip(header, f))
            if sample_hint and not d.get("SAMPLE"):
                d["SAMPLE"] = sample_hint
            if mutation_type_hint and not d.get("MUTATION_TYPE"):
                d["MUTATION_TYPE"] = mutation_type_hint
            rows.append(d)
    return rows


def infer_sample_from_classified_path(path: Path) -> Optional[str]:
    name = path.name
    m = re.match(r"(.+?)\.classified_df_(snvs|indels)\.tsv(?:\.gz)?$", name)
    if m:
        return m.group(1)
    # nested: <outroot>/<sample>/<snvs|indels>/classify/classified_df_snvs.tsv
    parts = path.parts
    for i, p in enumerate(parts):
        if p in {"snvs", "indels"} and i > 0:
            return parts[i - 1]
    return None


def infer_type_from_classified_path(path: Path) -> Optional[str]:
    s = path.name.lower()
    if "snv" in s:
        return "snvs"
    if "indel" in s:
        return "indels"
    for p in reversed(path.parts):
        if p in {"snvs", "indels"}:
            return p
    return None


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    for p in paths:
        if p and p.exists() and p.stat().st_size > 0:
            return p
    return None


def default_ffperase_dir(root: Path) -> Optional[Path]:
    candidates = [
        root / "ffperase_classification",
        root / "ffperase_post_v5_dedup_skip",
        root / "ffperase_post_v5_dedup",
        root / "ffperase_post_v4_plotfix",
        root / "ffperase_post",
    ]
    for c in candidates:
        if c.exists():
            if (c / "classified").exists() or list(c.glob("**/classified_df_snvs.tsv")) or list(c.glob("**/classified_df_indels.tsv")):
                return c
    return None


def collect_classified_files(base: Path, mutation_type: str, explicit: Optional[Path] = None) -> List[Path]:
    if explicit:
        return [explicit]
    filename = f"all_samples.ffperase_{mutation_type}.tsv"
    all_samples_candidates = [
        base / "classified" / filename,
        base / filename,
    ]
    found = first_existing(all_samples_candidates)
    if found:
        return [found]

    # Fall back to per-sample files.
    globs = [
        f"classified/*.classified_df_{mutation_type}.tsv",
        f"*.classified_df_{mutation_type}.tsv",
        f"*/{mutation_type}/classify/classified_df_{mutation_type}.tsv",
        f"*/{mutation_type}/classified_df_{mutation_type}.tsv",
        f"**/*.classified_df_{mutation_type}.tsv",
        f"**/classified_df_{mutation_type}.tsv",
    ]
    files: List[Path] = []
    seen = set()
    for pat in globs:
        for p in base.glob(pat):
            if p.exists() and p.stat().st_size > 0 and p.resolve() not in seen:
                seen.add(p.resolve())
                files.append(p)
    return sorted(files)


def choose_pred_col(row: dict, candidates: Sequence[str]) -> Optional[str]:
    for c in candidates:
        if c in row and not is_missing(row.get(c)):
            return c
    # Also match model-specific columns such as FFPErase_snvs_C10174_raw_predicts.
    for k in row:
        if k.endswith("_raw_predicts") and candidates is RAW_PREDICT_COLS and not is_missing(row.get(k)):
            return k
        if k.endswith("_predicts") and not k.endswith("_raw_predicts") and candidates is PREDICT_COLS and not is_missing(row.get(k)):
            return k
    return None


def row_probability(row: dict) -> Optional[float]:
    col = choose_pred_col(row, RAW_PREDICT_COLS)
    if not col:
        return None
    try:
        return float(str(row[col]).strip())
    except Exception:
        return None


def row_artifact(row: dict) -> Optional[int]:
    col = choose_pred_col(row, PREDICT_COLS)
    if not col:
        return None
    return bool_to_int(row[col])


def key_from_row(row: dict) -> Optional[Tuple[str, str, int, str, str]]:
    sample = row.get("SAMPLE") or row.get("sample")
    chrom = row.get("CHR") or row.get("CHROM") or row.get("chrom")
    pos = row.get("START") or row.get("POS") or row.get("pos")
    ref = row.get("REF") or row.get("ref")
    alt = row.get("ALT") or row.get("alt")
    if is_missing(sample) or is_missing(chrom) or is_missing(pos) or is_missing(ref) or is_missing(alt):
        return None
    try:
        pos_i = int(float(str(pos).strip()))
    except Exception:
        return None
    return (str(sample).strip(), norm_chrom(str(chrom)), pos_i, str(ref).strip().upper(), str(alt).strip().upper())


def build_lookup(paths: List[Path], mutation_type: str) -> Tuple[Dict[Tuple[str, str, int, str, str], dict], int, List[Path]]:
    lookup: Dict[Tuple[str, str, int, str, str], dict] = {}
    total = 0
    used_files: List[Path] = []
    for p in paths:
        sample_hint = infer_sample_from_classified_path(p)
        type_hint = infer_type_from_classified_path(p) or mutation_type
        rows = read_tsv_rows(p, sample_hint, type_hint)
        if rows:
            used_files.append(p)
        for r in rows:
            if not r.get("MUTATION_TYPE"):
                r["MUTATION_TYPE"] = mutation_type
            k = key_from_row(r)
            if not k:
                continue
            total += 1
            # Duplicate rows are resolved conservatively: keep artifact=True over
            # artifact=False; otherwise keep the higher raw probability row.
            if k in lookup:
                old = lookup[k]
                old_art = row_artifact(old)
                new_art = row_artifact(r)
                old_prob = row_probability(old)
                new_prob = row_probability(r)
                keep_new = False
                if old_art != 1 and new_art == 1:
                    keep_new = True
                elif old_art == new_art and new_prob is not None and (old_prob is None or new_prob > old_prob):
                    keep_new = True
                if keep_new:
                    lookup[k] = r
            else:
                lookup[k] = r
    return lookup, total, used_files


def parse_samples(samples_csv: str) -> Optional[set]:
    if not samples_csv:
        return None
    vals = [s.strip() for s in re.split(r"[,\s]+", samples_csv) if s.strip()]
    return set(vals)


def sample_from_vcf(path: Path) -> Optional[str]:
    n = path.name
    suffixes = [
        ".final_variants.vcf.gz",
        ".final_variants.vcf",
        ".vcf.gz",
        ".vcf",
    ]
    for suf in suffixes:
        if n.endswith(suf):
            return n[: -len(suf)]
    return None


def collect_vcfs(final_vcf_dir: Path, samples_csv: str) -> List[VcfJob]:
    selected = parse_samples(samples_csv)
    jobs_by_sample: Dict[str, Path] = {}

    # Prefer bgzipped VCFs over duplicate uncompressed VCFs.
    for p in sorted(final_vcf_dir.glob("*.final_variants.vcf.gz")):
        sample = sample_from_vcf(p)
        if sample and (selected is None or sample in selected):
            jobs_by_sample[sample] = p
    for p in sorted(final_vcf_dir.glob("*.final_variants.vcf")):
        sample = sample_from_vcf(p)
        if sample and (selected is None or sample in selected) and sample not in jobs_by_sample:
            jobs_by_sample[sample] = p

    return [VcfJob(s, jobs_by_sample[s]) for s in sorted(jobs_by_sample)]


def existing_info_ids(header_lines: List[str]) -> set:
    ids = set()
    for line in header_lines:
        m = re.match(r"##INFO=<ID=([^,>]+)", line)
        if m:
            ids.add(m.group(1))
    return ids


def existing_filter_ids(header_lines: List[str]) -> set:
    ids = set()
    for line in header_lines:
        m = re.match(r"##FILTER=<ID=([^,>]+)", line)
        if m:
            ids.add(m.group(1))
    return ids


def parse_info(info: str) -> List[str]:
    if info == "." or not info:
        return []
    return [x for x in info.split(";") if x]


def row_to_info(row: dict) -> List[str]:
    out: List[str] = []
    art = row_artifact(row)
    prob = row_probability(row)
    out.append("FFPERASE_EVALUATED=1")
    if art == 1:
        out.append("FFPERASE_DECISION=ARTIFACT")
        out.append("FFPERASE_ARTIFACT=1")
    elif art == 0:
        out.append("FFPERASE_DECISION=REAL")
        out.append("FFPERASE_ARTIFACT=0")
    if row.get("MUTATION_TYPE"):
        out.append(f"FFPERASE_MUTATION_TYPE={clean_string(row.get('MUTATION_TYPE'))}")
    if prob is not None:
        out.append(f"FFPERASE_PROB={prob:.6g}")

    for src, dest in FIELD_MAP.items():
        if src not in row:
            continue
        if dest in NUMERIC_FIELDS:
            val = clean_numeric(row.get(src), NUMERIC_FIELDS[dest])
        else:
            val = clean_string(row.get(src))
        if val is not None:
            out.append(f"{dest}={val}")
    return out


def append_filter(old_filter: str, new_filter: str) -> str:
    if not old_filter or old_filter == "." or old_filter == "PASS":
        return new_filter
    parts = old_filter.split(";")
    if new_filter not in parts:
        parts.append(new_filter)
    return ";".join(parts)


def annotate_one_vcf(
    sample: str,
    in_vcf: Path,
    out_vcf_plain: Path,
    lookup: Dict[Tuple[str, str, int, str, str], dict],
    annotate_missing: bool,
    filter_artifacts: bool,
) -> dict:
    stats = defaultdict(int)
    out_vcf_plain.parent.mkdir(parents=True, exist_ok=True)

    header_lines: List[str] = []
    chrom_line: Optional[str] = None
    # Read header first.
    with open_text(in_vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line.rstrip("\n\r"))
            elif line.startswith("#CHROM"):
                chrom_line = line.rstrip("\n\r")
                break
    if chrom_line is None:
        raise RuntimeError(f"No #CHROM header found in {in_vcf}")

    info_ids = existing_info_ids(header_lines)
    filter_ids = existing_filter_ids(header_lines)
    added_info = []
    for h in INFO_HEADERS:
        m = re.match(r"##INFO=<ID=([^,>]+)", h)
        if m and m.group(1) not in info_ids:
            added_info.append(h)
    added_filter = []
    if filter_artifacts and "FFPERASE_ARTIFACT" not in filter_ids:
        added_filter.append(FILTER_HEADER)

    with open_text(in_vcf) as fh, open(out_vcf_plain, "wt", encoding="utf-8", newline="") as out:
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                for h in added_info:
                    out.write(h + "\n")
                for h in added_filter:
                    out.write(h + "\n")
                out.write(line)
                continue
            if not line.strip():
                continue
            stats["vcf_records"] += 1
            f = line.rstrip("\n\r").split("\t")
            if len(f) < 8:
                out.write(line)
                continue
            chrom, pos_s, vid, ref, alts, qual, filt, info = f[:8]
            try:
                pos_i = int(pos_s)
            except Exception:
                out.write(line)
                continue
            alt_list = alts.split(",")
            # Final VCFs should already be split by bcftools norm -m -any. If not,
            # exact multi-allelic matching is attempted before falling back to first ALT.
            key_exact = (sample, norm_chrom(chrom), pos_i, ref.upper(), alts.upper())
            row = lookup.get(key_exact)
            if row is None and len(alt_list) == 1:
                key = (sample, norm_chrom(chrom), pos_i, ref.upper(), alt_list[0].upper())
                row = lookup.get(key)
            elif row is None and len(alt_list) > 1:
                key = (sample, norm_chrom(chrom), pos_i, ref.upper(), alt_list[0].upper())
                row = lookup.get(key)
                stats["multiallelic_records"] += 1

            extra_info = parse_info(info)
            if row is not None:
                stats["ffperase_evaluated"] += 1
                art = row_artifact(row)
                if art == 1:
                    stats["ffperase_artifact"] += 1
                    if filter_artifacts:
                        f[6] = append_filter(f[6], "FFPERASE_ARTIFACT")
                elif art == 0:
                    stats["ffperase_real"] += 1
                else:
                    stats["ffperase_unknown_prediction"] += 1
                extra_info.extend(row_to_info(row))
            else:
                stats["ffperase_missing"] += 1
                if annotate_missing:
                    extra_info.extend(["FFPERASE_EVALUATED=0", "FFPERASE_DECISION=NOT_EVALUATED"])
            f[7] = ";".join(extra_info) if extra_info else "."
            out.write("\t".join(f) + "\n")

    return dict(stats)


def bgzip_and_index(plain_vcf: Path, gz_vcf: Path, threads: int, keep_plain: bool) -> Optional[Path]:
    bgzip = shutil.which("bgzip")
    tabix = shutil.which("tabix")
    if not bgzip:
        eprint(f"[WARN] bgzip not found. Leaving uncompressed VCF: {plain_vcf}")
        return None
    gz_vcf.parent.mkdir(parents=True, exist_ok=True)
    tmp_gz = gz_vcf.with_suffix(gz_vcf.suffix + ".tmp")
    with open(tmp_gz, "wb") as out:
        cmd = [bgzip, "-@", str(max(1, threads)), "-c", str(plain_vcf)]
        subprocess.run(cmd, stdout=out, check=True)
    tmp_gz.replace(gz_vcf)
    if tabix:
        subprocess.run([tabix, "-f", "-p", "vcf", str(gz_vcf)], check=True)
    else:
        eprint(f"[WARN] tabix not found. No index produced for {gz_vcf}")
    if not keep_plain:
        try:
            plain_vcf.unlink()
        except FileNotFoundError:
            pass
    return gz_vcf


def copy_classified_files(ffperase_dir: Path, outdir: Path) -> None:
    dest = outdir / "classified"
    dest.mkdir(parents=True, exist_ok=True)
    candidates = []
    for mt in ["snvs", "indels"]:
        p = first_existing([
            ffperase_dir / "classified" / f"all_samples.ffperase_{mt}.tsv",
            ffperase_dir / f"all_samples.ffperase_{mt}.tsv",
        ])
        if p:
            candidates.append(p)
    for p in candidates:
        target = dest / p.name
        try:
            if target.exists() and p.resolve() == target.resolve():
                continue
        except OSError:
            pass
        shutil.copy2(p, target)


def write_summary(path: Path, rows: List[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample", "input_vcf", "annotated_vcf", "annotated_vcf_index",
        "vcf_records", "ffperase_evaluated", "ffperase_real", "ffperase_artifact",
        "ffperase_missing", "ffperase_unknown_prediction", "multiallelic_records",
    ]
    count_fields = set(fields[4:])
    with open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            fixed = {}
            for k in fields:
                if k in row:
                    fixed[k] = row[k]
                elif k in count_fields:
                    fixed[k] = 0
                else:
                    fixed[k] = ""
            w.writerow(fixed)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Annotate final_vcf/*.final_variants.vcf(.gz) with FFPErase TSV classification output.",
        epilog="""
Example
-------
python3 annotate_vcfs_with_ffperase.py \
  --root /path/to/pipeline_results \
  --ffperase-dir /path/to/pipeline_results/ffperase_post \
  --outdir /path/to/pipeline_results/ffperase_classification \
  --threads 4

Conservative artifact-filter tagging, while preserving all records:
python3 annotate_vcfs_with_ffperase.py --root ROOT --filter-artifacts
""",
    )
    p.add_argument("--root", required=True, help="Pipeline output root containing final_vcf/.")
    p.add_argument("--final-vcf-dir", default="", help="Final VCF directory. Default: <root>/final_vcf")
    p.add_argument("--ffperase-dir", default="", help="Existing FFPErase result dir. Default: auto-detect under <root>.")
    p.add_argument("--outdir", default="", help="Output dir. Default: <root>/ffperase_classification")
    p.add_argument("--samples", default="", help="Comma/space-separated sample IDs to annotate. Default: all final_vcf samples.")
    p.add_argument("--snvs-tsv", default="", help="Explicit all_samples or per-sample SNV classification TSV.")
    p.add_argument("--indels-tsv", default="", help="Explicit all_samples or per-sample indel classification TSV.")
    p.add_argument("--threads", type=int, default=2, help="Threads for bgzip. Default: 2")
    p.add_argument("--keep-uncompressed", action="store_true", help="Keep uncompressed annotated VCFs in addition to .vcf.gz.")
    p.add_argument("--filter-artifacts", action="store_true", help="Append FILTER=FFPERASE_ARTIFACT for variants predicted artifact. Default: preserve original FILTER.")
    p.add_argument("--no-annotate-missing", action="store_true", help="Do not add FFPERASE_EVALUATED=0 to variants missing from FFPErase TSVs.")
    p.add_argument("--no-copy-classified", action="store_true", help="Do not copy all_samples FFPErase TSVs into <outdir>/classified/.")
    p.add_argument("--force", action="store_true", help="Overwrite existing annotated VCFs.")
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    root = Path(args.root).resolve()
    final_vcf_dir = Path(args.final_vcf_dir).resolve() if args.final_vcf_dir else root / "final_vcf"
    outdir = Path(args.outdir).resolve() if args.outdir else root / "ffperase_classification"

    if not final_vcf_dir.exists():
        raise SystemExit(f"ERROR: final VCF directory not found: {final_vcf_dir}")

    if args.ffperase_dir:
        ffperase_dir = Path(args.ffperase_dir).resolve()
    else:
        detected = default_ffperase_dir(root)
        if not detected:
            raise SystemExit(
                "ERROR: could not auto-detect FFPErase result directory. "
                "Pass --ffperase-dir /path/to/ffperase_post_v5_dedup_skip"
            )
        ffperase_dir = detected.resolve()

    out_vcf_dir = outdir / "final_vcf"
    report_dir = outdir / "reports"
    out_vcf_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    snv_paths = collect_classified_files(ffperase_dir, "snvs", Path(args.snvs_tsv).resolve() if args.snvs_tsv else None)
    indel_paths = collect_classified_files(ffperase_dir, "indels", Path(args.indels_tsv).resolve() if args.indels_tsv else None)

    snv_lookup, snv_total, snv_used = build_lookup(snv_paths, "snvs")
    indel_lookup, indel_total, indel_used = build_lookup(indel_paths, "indels")
    lookup = dict(snv_lookup)
    for k, v in indel_lookup.items():
        lookup[k] = v

    jobs = collect_vcfs(final_vcf_dir, args.samples)
    if not jobs:
        raise SystemExit(f"ERROR: no final VCFs found in {final_vcf_dir}")

    eprint(f"Root             : {root}")
    eprint(f"Final VCF dir    : {final_vcf_dir}")
    eprint(f"FFPErase dir     : {ffperase_dir}")
    eprint(f"Output dir       : {outdir}")
    eprint(f"Samples          : {len(jobs)}")
    eprint(f"SNV TSV files    : {len(snv_used)}; rows indexed={len(snv_lookup)} / parsed={snv_total}")
    for p in snv_used:
        eprint(f"  - {p}")
    eprint(f"Indel TSV files  : {len(indel_used)}; rows indexed={len(indel_lookup)} / parsed={indel_total}")
    for p in indel_used:
        eprint(f"  - {p}")

    if not args.no_copy_classified:
        copy_classified_files(ffperase_dir, outdir)

    summary_rows: List[dict] = []
    for job in jobs:
        out_plain = out_vcf_dir / f"{job.sample}.final_variants.ffperase_annotated.vcf"
        out_gz = out_vcf_dir / f"{job.sample}.final_variants.ffperase_annotated.vcf.gz"
        if out_gz.exists() and not args.force:
            eprint(f"SKIP existing: {out_gz}")
            summary_rows.append({
                "sample": job.sample,
                "input_vcf": str(job.vcf),
                "annotated_vcf": str(out_gz),
                "annotated_vcf_index": str(out_gz) + ".tbi" if Path(str(out_gz) + ".tbi").exists() else "",
            })
            continue
        eprint(f"Annotating {job.sample}: {job.vcf}")
        stats = annotate_one_vcf(
            sample=job.sample,
            in_vcf=job.vcf,
            out_vcf_plain=out_plain,
            lookup=lookup,
            annotate_missing=not args.no_annotate_missing,
            filter_artifacts=args.filter_artifacts,
        )
        gz = bgzip_and_index(out_plain, out_gz, args.threads, args.keep_uncompressed)
        annotated_path = gz if gz else out_plain
        idx_path = Path(str(annotated_path) + ".tbi")
        row = {
            "sample": job.sample,
            "input_vcf": str(job.vcf),
            "annotated_vcf": str(annotated_path),
            "annotated_vcf_index": str(idx_path) if idx_path.exists() else "",
        }
        row.update(stats)
        summary_rows.append(row)
        eprint(
            f"  records={stats.get('vcf_records',0)} evaluated={stats.get('ffperase_evaluated',0)} "
            f"artifact={stats.get('ffperase_artifact',0)} real={stats.get('ffperase_real',0)} missing={stats.get('ffperase_missing',0)}"
        )

    summary_path = report_dir / "all_samples.ffperase_annotation_summary.tsv"
    write_summary(summary_path, summary_rows)
    eprint(f"Summary: {summary_path}")
    eprint(f"Done. Annotated VCFs: {out_vcf_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
