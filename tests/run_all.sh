#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
export PYTHONDONTWRITEBYTECODE=1

overall_status=0

run_test() {
  local label="$1"
  shift
  printf '\n== %s ==\n' "$label"
  if "$@"; then
    printf 'PASS: %s\n' "$label"
  else
    printf 'FAIL: %s\n' "$label" >&2
    overall_status=1
  fi
}

run_test 'parser launcher unit tests' "${TEST_DIR}/test_run_pipeline.sh"
run_test 'package static rules' python3 "${TEST_DIR}/test_package_rules.py"
run_test 'BAM/sample resolution' python3 "${TEST_DIR}/test_resolve_bams.py"
run_test 'bundled synthetic examples' "${TEST_DIR}/test_example_data.sh"
run_test 'fresh versus FFPE unit behavior' python3 "${TEST_DIR}/test_sample_mode.py"
run_test 'VCF left-alignment integration' "${TEST_DIR}/test_left_alignment.sh"

if [[ "${SKIP_REAL_PARSER_SMOKE:-false}" == "true" ]]; then
  printf '\nSKIP: real Nextflow parser smoke (SKIP_REAL_PARSER_SMOKE=true)\n'
else
  run_test 'real Nextflow v1/v2 parser smoke' "${TEST_DIR}/test_parser_smoke.sh"
fi

if [[ "${RUN_FRESH_PIPELINE_SMOKE:-false}" == "true" ]]; then
  run_test 'optional real fresh/BQSR pipeline smoke' "${TEST_DIR}/test_fresh_pipeline_smoke.sh"
else
  printf '\nSKIP: optional fresh/BQSR workflow smoke (set RUN_FRESH_PIPELINE_SMOKE=true)\n'
fi

if [[ "${RUN_FFPE_PIPELINE_SMOKE:-false}" == "true" ]]; then
  run_test 'optional FFPE empty-variant pipeline smoke' "${TEST_DIR}/test_ffpe_pipeline_smoke.sh"
else
  printf '\nSKIP: optional FFPE empty-variant smoke (set RUN_FFPE_PIPELINE_SMOKE=true)\n'
fi

exit "$overall_status"
