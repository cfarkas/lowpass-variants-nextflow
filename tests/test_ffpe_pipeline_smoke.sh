#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${TEST_DIR}/.." && pwd -P)"
EXAMPLE="${PACKAGE_DIR}/examples/run_minimal_ffpe.sh"

fail() {
  printf 'not ok - %s\n' "$*" >&2
  exit 1
}

skip() {
  printf 'ok - SKIP FFPE empty-variant pipeline smoke: %s\n' "$*"
  exit 0
}

if [[ -n "${LOWPASS_TEST_TOOL_ENV:-}" ]]; then
  [[ -d "${LOWPASS_TEST_TOOL_ENV}/bin" ]] ||
    fail "LOWPASS_TEST_TOOL_ENV has no bin directory: ${LOWPASS_TEST_TOOL_ENV}"
  export PATH="${LOWPASS_TEST_TOOL_ENV}/bin:${PATH}"
fi

for executable in gatk samtools bcftools bgzip tabix python nextflow timeout; do
  command -v "$executable" >/dev/null 2>&1 ||
    skip "required executable not found: ${executable}"
done

engine_found=false
for engine in apptainer singularity docker; do
  if command -v "$engine" >/dev/null 2>&1; then
    engine_found=true
    break
  fi
done
[[ "$engine_found" == true ]] || skip "no Apptainer, Singularity, or Docker executable found"
[[ -x "$EXAMPLE" ]] || fail "FFPE example is missing or not executable: ${EXAMPLE}"

TEST_TMP="$(mktemp -d)"
trap 'rm -rf -- "$TEST_TMP"' EXIT
OUTDIR="${TEST_TMP}/results"
RUN_LOG="${TEST_TMP}/ffpe-example.log"

status=0
set +e
(
  cd "$TEST_TMP"
  timeout "${FFPE_SMOKE_TIMEOUT:-300}" \
    env \
      NEXTFLOW_PROFILE= \
      OUTDIR="$OUTDIR" \
      NXF_OFFLINE=true \
      NXF_DISABLE_CHECK_LATEST=true \
      bash "$EXAMPLE"
) > "$RUN_LOG" 2>&1
status=$?
set -e

if [[ "$status" -ne 0 ]]; then
  sed -n '1,320p' "$RUN_LOG" >&2
  fail "FFPE empty-variant example failed with status ${status}"
fi

DONE_FILE="${OUTDIR}/logs/ffperase.done.txt"
ANNOTATED_VCF="${OUTDIR}/ffperase_classification/final_vcf/SAMPLE.final_variants.ffperase_annotated.vcf.gz"
STATUS_DIR="${OUTDIR}/ffperase_status"

[[ -s "$DONE_FILE" ]] || fail "missing FFPE completion summary"
grep -Fqx $'status\tcompleted' "$DONE_FILE" || fail "FFPE completion status is not completed"
grep -Fqx $'failed_sample_types\t0' "$DONE_FILE" || fail "FFPE smoke reported failures"
[[ -s "$ANNOTATED_VCF" && -s "${ANNOTATED_VCF}.tbi" ]] ||
  fail "missing indexed annotated FFPE VCF"
bcftools index --stats "$ANNOTATED_VCF" >/dev/null 2>&1 ||
  fail "annotated FFPE VCF index is invalid"

mapfile -t status_files < <(find "$STATUS_DIR" -maxdepth 1 -type f -name 'SAMPLE.*.status.tsv' | sort)
[[ "${#status_files[@]}" -eq 2 ]] || fail "expected two FFPE sample/type status files"
for status_file in "${status_files[@]}"; do
  awk -F '\t' 'NR == 2 && $3 == "SKIP_NO_VARIANTS" {ok=1} END {exit(ok ? 0 : 1)}' "$status_file" ||
    fail "unexpected FFPE status in ${status_file}"
done

printf 'ok - bundled FFPE example completes BQSR, Picard metrics, empty-variant statuses, and annotation\n'
