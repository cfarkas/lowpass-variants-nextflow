#!/usr/bin/env bash
set -Eeuo pipefail

# Portable FFPE template. Export the required paths before running:
#   INPUT_DIR=/data/bams OUTDIR=/data/results REF=/refs/genome.fa \
#   KNOWN_SITES=/refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz \
#   bash examples/lymphoma_genomewide_simplified.sh

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
  --ffpe \
  --input "${INPUT_DIR}" \
  --outdir "${OUTDIR}" \
  --samples "${SAMPLES}" \
  --ref "${REF}" \
  --known_sites "${KNOWN_SITES}" \
  --ffperase_threads "${FFPERASE_THREADS:-16}" \
  --ffperase_classify_threads "${FFPERASE_CLASSIFY_THREADS:-2}" \
  --ffperase_sample_jobs "${FFPERASE_SAMPLE_JOBS:-1}" \
  --ffperase_picard_jobs "${FFPERASE_PICARD_JOBS:-1}" \
  --ffperase_split_pileup "${FFPERASE_SPLIT_PILEUP:-5000}" \
  --ffperase_split_reads "${FFPERASE_SPLIT_READS:-7500000}" \
  --ffperase_engine "${FFPERASE_ENGINE:-auto}" \
  --ffperase_filter_artifacts "${FFPERASE_FILTER_ARTIFACTS:-true}" \
  --ffperase_fail_on_error "${FFPERASE_FAIL_ON_ERROR:-true}"
