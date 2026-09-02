#!/usr/bin/env python3
import argparse
import csv
import math
import os
import sys
from typing import Optional, Tuple

TAB = "\t"


def truthy(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    return str(v).strip().lower() in {"true", "1", "yes", "y", "si", "sí"}


def safe_float(x) -> Optional[float]:
    try:
        if x is None or str(x).strip() in {"", "."}:
            return None
        return float(x)
    except Exception:
        return None


def weighted_mean_from_samtools_coverage(path: str) -> Tuple[Optional[float], str]:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, "missing_genome_coverage"
    total_bases = 0.0
    weighted_depth = 0.0
    used = 0
    with open(path, "rt") as fh:
        header = None
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                # samtools coverage header usually starts with #rname
                header = [x.lstrip("#") for x in line.split()]
                continue
            fields = line.split()
            if len(fields) < 7:
                continue
            # samtools coverage columns: rname startpos endpos numreads covbases coverage meandepth ...
            chrom = fields[0]
            if chrom == "*":
                continue
            try:
                start = int(float(fields[1])); end = int(float(fields[2])); meandepth = float(fields[6])
            except Exception:
                continue
            length = max(0, end - start + 1)
            if length <= 0:
                continue
            total_bases += length
            weighted_depth += meandepth * length
            used += 1
    if total_bases <= 0:
        return None, "no_usable_genome_coverage_rows"
    return weighted_depth / total_bases, f"samtools_coverage_weighted_contigs={used}"


def weighted_mean_from_bedtools_mean(path: str) -> Tuple[Optional[float], str]:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, "missing_target_coverage"
    total_bases = 0.0
    weighted_depth = 0.0
    used = 0
    with open(path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                continue
            try:
                start = int(float(fields[1])); end = int(float(fields[2])); mean_depth = float(fields[-1])
            except Exception:
                continue
            length = max(0, end - start)
            if length <= 0:
                continue
            total_bases += length
            weighted_depth += mean_depth * length
            used += 1
    if total_bases <= 0:
        return None, "no_usable_target_coverage_rows"
    return weighted_depth / total_bases, f"bedtools_coverage_mean_weighted_intervals={used}"


def choose_thresholds(
    mean_depth: Optional[float],
    manual_min_dp: int,
    manual_alt: int,
    manual_af: float,
    manual_vote: int,
    manual_ffpe_vaf: Optional[float],
    manual_ffpe_alt: Optional[int],
    auto: bool,
    sample_mode: str,
) -> dict:
    ffpe_mode = sample_mode == "ffpe"
    if not auto or mean_depth is None or mean_depth <= 0:
        thresholds = {
            "threshold_mode": "manual" if not auto else "manual_fallback_no_coverage",
            "effective_min_dp": int(manual_min_dp),
            "effective_min_alt_reads": int(manual_alt),
            "effective_min_af": float(manual_af),
            "effective_vote_threshold": int(manual_vote),
            "heuristic_note": "using_user_parameters",
        }
        if ffpe_mode:
            thresholds.update(
                effective_ffpe_ct_vaf_keep=float(manual_ffpe_vaf),
                effective_ffpe_ct_min_alt=int(manual_ffpe_alt),
            )
        return thresholds

    d = float(mean_depth)
    # Common low-pass thresholds adapt to coverage in either sample mode. FFPE
    # guard thresholds are populated only for FFPE samples.
    if d < 2.0:
        min_dp, min_alt, min_af, note = 3, 2, 0.35, "ultra_low_coverage_high_AF_guard"
        if ffpe_mode:
            ffpe_vaf, ffpe_alt = 0.55, 3
    elif d < 5.0:
        min_dp, min_alt, min_af, note = 4, 2, 0.30, "very_low_coverage_high_AF_guard"
        if ffpe_mode:
            ffpe_vaf, ffpe_alt = 0.50, 3
    elif d < 10.0:
        min_dp, min_alt, min_af, note = 6, 3, 0.25, "low_coverage_moderate_AF_guard"
        if ffpe_mode:
            ffpe_vaf, ffpe_alt = 0.45, 5
    elif d < 20.0:
        min_dp, min_alt, min_af, note = 8, 4, 0.20, "moderate_lowpass_default_like"
        if ffpe_mode:
            ffpe_vaf, ffpe_alt = 0.35, 8
    elif d < 40.0:
        min_dp, min_alt, min_af, note = 10, 5, 0.15, "higher_depth_more_sensitive_AF"
        if ffpe_mode:
            ffpe_vaf, ffpe_alt = 0.30, 10
    else:
        min_dp, min_alt, min_af, note = 12, 6, 0.10, "high_depth_more_sensitive_AF"
        if ffpe_mode:
            ffpe_vaf, ffpe_alt = 0.25, 12

    thresholds = {
        "threshold_mode": "auto_coverage_heuristic",
        "effective_min_dp": min_dp,
        "effective_min_alt_reads": min_alt,
        "effective_min_af": min_af,
        "effective_vote_threshold": int(manual_vote),
        "heuristic_note": note,
    }
    if ffpe_mode:
        thresholds.update(
            effective_ffpe_ct_vaf_keep=ffpe_vaf,
            effective_ffpe_ct_min_alt=ffpe_alt,
        )
    return thresholds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sample-mode", choices=("ffpe", "fresh"), required=True)
    ap.add_argument("--target-mode", required=True)
    ap.add_argument("--genome-coverage", default="")
    ap.add_argument("--target-coverage", default="")
    ap.add_argument("--auto-thresholds", default="true")
    ap.add_argument("--manual-min-dp", type=int, required=True)
    ap.add_argument("--manual-min-alt-reads", type=int, required=True)
    ap.add_argument("--manual-min-af", type=float, required=True)
    ap.add_argument("--manual-vote-threshold", type=int, required=True)
    ap.add_argument("--manual-ffpe-ct-vaf-keep", type=float)
    ap.add_argument("--manual-ffpe-ct-min-alt", type=int)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if args.sample_mode == "ffpe":
        if args.manual_ffpe_ct_vaf_keep is None or args.manual_ffpe_ct_min_alt is None:
            ap.error("FFPE mode requires both --manual-ffpe-ct-vaf-keep and --manual-ffpe-ct-min-alt")

    if args.target_mode != "genome_wide":
        mean_depth, coverage_source = weighted_mean_from_bedtools_mean(args.target_coverage)
        if mean_depth is None:
            mean_depth, coverage_source2 = weighted_mean_from_samtools_coverage(args.genome_coverage)
            coverage_source = coverage_source + ";fallback_" + coverage_source2
    else:
        mean_depth, coverage_source = weighted_mean_from_samtools_coverage(args.genome_coverage)

    thresholds = choose_thresholds(
        mean_depth,
        args.manual_min_dp,
        args.manual_min_alt_reads,
        args.manual_min_af,
        args.manual_vote_threshold,
        args.manual_ffpe_ct_vaf_keep,
        args.manual_ffpe_ct_min_alt,
        truthy(args.auto_thresholds),
        args.sample_mode,
    )
    mean_depth_str = "." if mean_depth is None else f"{mean_depth:.6f}"
    fields = [
        "sample", "target_mode", "coverage_source", "mean_depth", "threshold_mode",
        "effective_min_dp", "effective_min_alt_reads", "effective_min_af", "effective_vote_threshold",
        "heuristic_note",
    ]
    if args.sample_mode == "ffpe":
        fields.extend(["effective_ffpe_ct_vaf_keep", "effective_ffpe_ct_min_alt"])
    row = {
        "sample": args.sample,
        "target_mode": args.target_mode,
        "coverage_source": coverage_source,
        "mean_depth": mean_depth_str,
        **thresholds,
    }
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, delimiter=TAB, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerow(row)
    message = (
        f"[CALIBRATE_THRESHOLDS] sample={args.sample} target_mode={args.target_mode} "
        f"sample_mode={args.sample_mode} "
        f"mean_depth={mean_depth_str} mode={row['threshold_mode']} "
        f"min_dp={row['effective_min_dp']} min_alt_reads={row['effective_min_alt_reads']} "
        f"min_af={row['effective_min_af']} vote_threshold={row['effective_vote_threshold']} "
        f"note={row['heuristic_note']}"
    )
    if args.sample_mode == "ffpe":
        message += (
            f" ffpe_ct_vaf_keep={row['effective_ffpe_ct_vaf_keep']}"
            f" ffpe_ct_min_alt={row['effective_ffpe_ct_min_alt']}"
        )
    print(message)


if __name__ == "__main__":
    main()
