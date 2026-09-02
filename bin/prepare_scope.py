#!/usr/bin/env python3
"""Prepare genome-wide/custom BED/gene target files for the low-pass variant workflow."""
import argparse
import gzip
import os
import re
import subprocess
from collections import defaultdict
from pathlib import Path


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def parse_attrs(attr_string):
    out = {}
    for part in attr_string.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r'(\S+)\s+"([^"]+)"', part)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def chrom_sort_key(chrom):
    c = str(chrom)
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    order = {str(i): i for i in range(1, 23)}
    order.update({"X": 23, "Y": 24, "MT": 25})
    return (order.get(c, 999), c)


def merge_intervals(rows3):
    rows3 = sorted(rows3, key=lambda x: (chrom_sort_key(x[0]), x[0], int(x[1]), int(x[2])))
    merged = []
    for chrom, start, end in rows3:
        start, end = int(start), int(end)
        if end <= start:
            continue
        if not merged or chrom != merged[-1][0] or start > merged[-1][2]:
            merged.append([chrom, start, end])
        else:
            merged[-1][2] = max(merged[-1][2], end)
    return merged


def write_rows(path, rows):
    with open(path, "w") as out:
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")


def read_custom_bed(path):
    rows = []
    with open(path, "rt", errors="replace") as handle:
        n = 0
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except Exception:
                continue
            if end <= start:
                continue
            chrom = fields[0]
            region_id = fields[3] if len(fields) > 3 and fields[3] else f"custom_region_{n+1}"
            region_type = fields[4] if len(fields) > 4 and fields[4] else "custom_bed"
            gene = fields[5] if len(fields) > 5 and fields[5] else region_id
            rows.append((chrom, start, end, region_id, region_type, gene))
            n += 1
    if not rows:
        raise SystemExit(f"ERROR: no valid intervals found in BED: {path}")
    return sorted(rows, key=lambda x: (chrom_sort_key(x[0]), x[0], int(x[1]), int(x[2]), x[3]))


def download_gencode(gtf_path: Path, version: str):
    url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{version}/gencode.v{version}.annotation.gtf.gz"
    gtf_path.parent.mkdir(parents=True, exist_ok=True)
    if gtf_path.exists() and gtf_path.stat().st_size > 0:
        return str(gtf_path)
    subprocess.run(["wget", "-O", str(gtf_path), url], check=True)
    return str(gtf_path)


def build_gene_targets(gtf, genes_csv, upstream, downstream):
    genes = [x.strip() for x in genes_csv.replace(";", ",").split(",") if x.strip()]
    if not genes:
        raise SystemExit("ERROR: --genes was provided but no valid gene symbols were found")
    want = set(genes)
    gene_info = {}
    exons = defaultdict(list)
    with open_maybe_gzip(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attrs = f
            attr = parse_attrs(attrs)
            gene = attr.get("gene_name")
            if not gene or gene not in want:
                continue
            start, end = int(start), int(end)
            if feature == "gene":
                gene_info[gene] = (chrom, start, end, strand)
            elif feature == "exon":
                exons[gene].append((chrom, start, end, strand))
    missing = sorted(want - set(gene_info))
    if missing:
        print("WARNING: genes missing in GTF: " + ",".join(missing))
    rows = []
    for gene in genes:
        if gene not in gene_info:
            continue
        chrom, start, end, strand = gene_info[gene]
        if strand == "-":
            tss = end
            p_start = max(0, tss - int(downstream) - 1)
            p_end = tss + int(upstream)
        else:
            tss = start
            p_start = max(0, tss - int(upstream) - 1)
            p_end = tss + int(downstream)
        rows.append((chrom, p_start, p_end, f"{gene}|promoter", "promoter", gene))
        exon_3 = [(c, max(0, s - 1), e) for c, s, e, _st in exons.get(gene, [])]
        for idx, (c, s, e) in enumerate(merge_intervals(exon_3), start=1):
            rows.append((c, s, e, f"{gene}|exon_{idx}", "exon", gene))
    if not rows:
        raise SystemExit("ERROR: no target rows were created from --genes")
    return sorted(rows, key=lambda x: (chrom_sort_key(x[0]), x[0], int(x[1]), int(x[2]), x[3]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="")
    ap.add_argument("--bed", default="")
    ap.add_argument("--gtf", default="")
    ap.add_argument("--gencode-version", default="45")
    ap.add_argument("--promoter-upstream", type=int, default=2000)
    ap.add_argument("--promoter-downstream", type=int, default=500)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    annot = outdir / "targets.annotated.bed"
    merged = outdir / "targets.merged.bed"
    scope = outdir / "scope.env"
    if args.genes and args.bed:
        raise SystemExit("ERROR: use only one of --genes or --bed")
    if args.bed:
        mode = "custom_bed"
        rows = read_custom_bed(args.bed)
    elif args.genes:
        mode = "genes"
        if args.gtf:
            if not os.path.exists(args.gtf):
                raise SystemExit(f"ERROR: GTF not found: {args.gtf}")
            gtf = args.gtf
        else:
            gtf = download_gencode(outdir / f"gencode.v{args.gencode_version}.annotation.gtf.gz", args.gencode_version)
        rows = build_gene_targets(gtf, args.genes, args.promoter_upstream, args.promoter_downstream)
    else:
        mode = "genome_wide"
        rows = []
    write_rows(annot, rows)
    merged_rows = merge_intervals([(r[0], r[1], r[2]) for r in rows]) if rows else []
    write_rows(merged, merged_rows)
    with open(scope, "w") as out:
        out.write(f"TARGET_MODE={mode}\n")
        out.write(f"TARGET_ANNOT_BED={annot.name}\n")
        out.write(f"TARGET_MERGED_BED={merged.name}\n")
    print(f"Prepared target scope: {mode}; intervals={len(rows)}; merged_intervals={len(merged_rows)}")


if __name__ == "__main__":
    main()
