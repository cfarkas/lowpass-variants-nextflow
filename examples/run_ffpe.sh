#!/usr/bin/env bash
set -Eeuo pipefail

# Complete FFPE run with BQSR and FFPErase.
#
# Required environment values:
#   INPUT_PATH=/data/bams
#   OUTDIR=/data/results/ffpe
#   REF=/refs/genome.fa
#   KNOWN_SITES=/refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz
#
# Optional values:
#   SAMPLES=SAMPLE_01,SAMPLE_02
#   WORK_DIR=/scratch/project/nextflow-work
#
# -work-dir and -resume are optional Nextflow controls. This production
# template intentionally includes them so an interrupted run can reuse completed
# tasks. Remove those two command lines to use Nextflow's default work directory
# and start without resume.

EXAMPLE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${EXAMPLE_DIR}/.." && pwd -P)"

: "${INPUT_PATH:?Set INPUT_PATH to a BAM directory, BAM, manifest, or BAM list}"
: "${OUTDIR:?Set OUTDIR to the output directory}"
: "${REF:?Set REF to the reference FASTA}"
: "${KNOWN_SITES:?Set KNOWN_SITES to comma-separated indexed VCFs}"

WORK_DIR="${WORK_DIR:-${OUTDIR}.work}"
SAMPLE_ARGS=()
if [[ -n "${SAMPLES:-}" ]]; then
  SAMPLE_ARGS=(--samples "${SAMPLES}")
fi

"${PACKAGE_DIR}/bin/run_pipeline.sh" \
  -profile conda \
  -work-dir "${WORK_DIR}" \
  -resume \
  --ffpe \
  --input "${INPUT_PATH}" \
  --outdir "${OUTDIR}" \
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
  --ffperase_fail_on_error "${FFPERASE_FAIL_ON_ERROR:-true}" \
  "${SAMPLE_ARGS[@]}"
