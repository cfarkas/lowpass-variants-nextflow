#!/usr/bin/env bash
set -Eeuo pipefail

# Complete fresh-material run with BQSR.
#
# Required environment values:
#   INPUT_PATH=/data/bams
#   OUTDIR=/data/results/fresh
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
  --fresh \
  --input "${INPUT_PATH}" \
  --outdir "${OUTDIR}" \
  --ref "${REF}" \
  --known_sites "${KNOWN_SITES}" \
  "${SAMPLE_ARGS[@]}"
