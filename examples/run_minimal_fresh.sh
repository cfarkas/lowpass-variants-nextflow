#!/usr/bin/env bash
set -Eeuo pipefail

EXAMPLE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${EXAMPLE_DIR}/.." && pwd -P)"
DATA_DIR="${EXAMPLE_DIR}/tiny"
OUTDIR="${OUTDIR:-${PACKAGE_DIR}/example-results/minimal-fresh}"
NEXTFLOW_PROFILE="${NEXTFLOW_PROFILE-conda}"

PROFILE_ARGS=()
if [[ -n "$NEXTFLOW_PROFILE" ]]; then
  PROFILE_ARGS=(-profile "$NEXTFLOW_PROFILE")
fi

printf 'Running bundled one-variant fresh quickstart\n'
printf '  input:  %s\n' "$DATA_DIR/SAMPLE.bam"
printf '  output: %s\n' "$OUTDIR"

"${PACKAGE_DIR}/bin/run_pipeline.sh" \
  "${PROFILE_ARGS[@]}" \
  --fresh \
  --input "${DATA_DIR}/SAMPLE.bam" \
  --outdir "$OUTDIR" \
  --ref "${DATA_DIR}/tiny.fa" \
  --known_sites "${DATA_DIR}/known-sites.vcf.gz" \
  --overwrite true \
  --skip_mutect2 true \
  --auto_thresholds false \
  --min_dp 10 \
  --min_alt_reads 4 \
  --min_af 0.20 \
  --vote_threshold 2 \
  --memory '4 GB' \
  --mutect2_memory '4 GB' \
  --freebayes_memory '2 GB' \
  --bcftools_memory '2 GB' \
  --normalize_memory '2 GB' \
  --finalize_memory '2 GB' \
  --threads 1 \
  --prepare_cpus 1 \
  --mutect2_cpus 1 \
  --freebayes_cpus 1 \
  --bcftools_cpus 1 \
  --normalize_cpus 1 \
  --finalize_cpus 1 \
  --max_samples_parallel 1 \
  --max_prepare_parallel 1 \
  --max_mutect2_parallel 1 \
  --max_freebayes_parallel 1 \
  --max_bcftools_parallel 1 \
  --max_normalize_parallel 1 \
  --max_finalize_parallel 1

printf 'Fresh quickstart result (one SNV at chr1:25): %s\n' \
  "$OUTDIR/final_vcf/SAMPLE.final_variants.vcf.gz"
