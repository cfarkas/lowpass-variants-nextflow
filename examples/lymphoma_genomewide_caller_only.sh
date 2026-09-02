#!/usr/bin/env bash
set -Eeuo pipefail

# Portable fresh-material template. Export the required paths before running:
#   INPUT_DIR=/data/bams OUTDIR=/data/results REF=/refs/genome.fa \
#   KNOWN_SITES=/refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz \
#   bash examples/lymphoma_genomewide_caller_only.sh

EXAMPLE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${EXAMPLE_DIR}/.." && pwd -P)"

: "${INPUT_DIR:?Set INPUT_DIR to the BAM input path}"
: "${OUTDIR:?Set OUTDIR to the output directory}"
: "${REF:?Set REF to the reference FASTA}"
: "${KNOWN_SITES:?Set KNOWN_SITES to comma-separated indexed VCFs}"

WORK_DIR="${WORK_DIR:-${OUTDIR}.work}"
SAMPLES="${SAMPLES:-SAMPLE_01,SAMPLE_02,SAMPLE_03}"

"${PACKAGE_DIR}/bin/run_pipeline.sh" \
  -profile conda \
  -resume \
  -work-dir "${WORK_DIR}" \
  --fresh \
  --input "${INPUT_DIR}" \
  --outdir "${OUTDIR}" \
  --samples "${SAMPLES}" \
  --ref "${REF}" \
  --known_sites "${KNOWN_SITES}"
