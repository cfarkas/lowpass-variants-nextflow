#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${TEST_DIR}/.." && pwd -P)"
LAUNCHER="${PACKAGE_DIR}/bin/run_pipeline.sh"
FAKE_BIN_DIR="${TEST_DIR}/fixtures"
PIPELINE="${PACKAGE_DIR}/main.nf"
TEST_TMP="$(mktemp -d)"
trap 'rm -rf -- "$TEST_TMP"' EXIT

fail() {
  printf 'not ok - %s\n' "$*" >&2
  exit 1
}

assert_contains() {
  local text="$1"
  local expected="$2"
  [[ "$text" == *"$expected"* ]] || fail "expected output to contain: ${expected}"
}

assert_line_count() {
  local file="$1"
  local expected_line="$2"
  local expected_count="$3"
  local actual_count
  actual_count="$(grep -Fxc -- "$expected_line" "$file" || true)"
  [[ "$actual_count" == "$expected_count" ]] || fail \
    "expected ${expected_count} occurrence(s) of '${expected_line}', found ${actual_count}"
}

run_count() {
  awk -F'|' '$2 == "run" { count++ } END { print count + 0 }' "$1"
}

test_autodetects_v2_and_forwards_arguments() {
  local log="${TEST_TMP}/autodetect-v2.log"
  local output
  output="$(
    env -u NXF_SYNTAX_PARSER \
      NEXTFLOW_BIN="${FAKE_BIN_DIR}/nextflow" \
      PATH="${FAKE_BIN_DIR}:${PATH}" \
      FAKE_NF_LOG="$log" \
      FAKE_NF_MODE=v2_only \
      "$LAUNCHER" --sample 'sample one' 'literal*glob' 2>&1
  )" || fail "v2 autodetection launch failed"

  assert_line_count "$log" "v2|run|${PIPELINE}|--help" 1
  assert_line_count "$log" "v2|run|${PIPELINE}|--sample|sample\\ one|literal\\*glob" 1
  [[ "$(run_count "$log")" == 2 ]] || fail "v2 detection should perform one probe and one launch"
  assert_contains "$output" 'nextflow version 99.1-test build 1'
  assert_contains "$output" 'Syntax parser: v2 (autodetected)'
  printf 'ok - autodetects v2 and forwards arguments without word splitting\n'
}

test_falls_back_to_v1() {
  local log="${TEST_TMP}/fallback-v1.log"
  local output
  output="$(
    env -u NXF_SYNTAX_PARSER \
      NEXTFLOW_BIN="${FAKE_BIN_DIR}/nextflow" \
      PATH="${FAKE_BIN_DIR}:${PATH}" \
      FAKE_NF_LOG="$log" \
      FAKE_NF_MODE=v1_only \
      "$LAUNCHER" --fresh 2>&1
  )" || fail "v1 fallback launch failed"

  assert_line_count "$log" "v2|run|${PIPELINE}|--help" 1
  assert_line_count "$log" "v1|run|${PIPELINE}|--help" 1
  assert_line_count "$log" "v1|run|${PIPELINE}|--fresh" 1
  [[ "$(run_count "$log")" == 3 ]] || fail "v1 fallback should perform two probes and one launch"
  assert_contains "$output" 'Syntax parser: v1 (autodetected)'
  printf 'ok - falls back from v2 to v1\n'
}

test_respects_explicit_valid_parser() {
  local log="${TEST_TMP}/explicit-v1.log"
  local output
  output="$(
    env \
      NEXTFLOW_BIN="${FAKE_BIN_DIR}/nextflow" \
      PATH="${FAKE_BIN_DIR}:${PATH}" \
      FAKE_NF_LOG="$log" \
      FAKE_NF_MODE=both \
      NXF_SYNTAX_PARSER=v1 \
      "$LAUNCHER" --ffpe 2>&1
  )" || fail "explicit v1 launch failed"

  assert_line_count "$log" "v1|run|${PIPELINE}|--ffpe" 1
  [[ "$(run_count "$log")" == 1 ]] || fail "an explicit valid parser must bypass probing"
  assert_contains "$output" 'Syntax parser: v1 (environment)'
  printf 'ok - respects an explicit valid parser\n'
}

test_invalid_parser_is_replaced_by_autodetection() {
  local log="${TEST_TMP}/invalid-parser.log"
  local output
  output="$(
    env \
      NEXTFLOW_BIN="${FAKE_BIN_DIR}/nextflow" \
      PATH="${FAKE_BIN_DIR}:${PATH}" \
      FAKE_NF_LOG="$log" \
      FAKE_NF_MODE=v2_only \
      NXF_SYNTAX_PARSER=legacy \
      "$LAUNCHER" --fresh 2>&1
  )" || fail "autodetection after invalid parser failed"

  assert_contains "$output" 'ignoring invalid NXF_SYNTAX_PARSER=legacy'
  assert_contains "$output" 'Syntax parser: v2 (autodetected)'
  assert_line_count "$log" "v2|run|${PIPELINE}|--fresh" 1
  printf 'ok - replaces an invalid parser value by autodetection\n'
}

test_reports_when_neither_parser_works() {
  local log="${TEST_TMP}/neither-parser.log"
  local output
  local status=0
  output="$(
    env -u NXF_SYNTAX_PARSER \
      NEXTFLOW_BIN="${FAKE_BIN_DIR}/nextflow" \
      PATH="${FAKE_BIN_DIR}:${PATH}" \
      FAKE_NF_LOG="$log" \
      FAKE_NF_MODE=neither \
      "$LAUNCHER" --fresh 2>&1
  )" || status=$?

  [[ "$status" -ne 0 ]] || fail "launcher unexpectedly succeeded when neither parser worked"
  [[ "$(run_count "$log")" == 2 ]] || fail "failed detection should stop after the two probes"
  assert_contains "$output" 'did not pass --help with either syntax parser'
  printf 'ok - stops cleanly when neither parser can compile the package help\n'
}

test_autodetects_v2_and_forwards_arguments
test_falls_back_to_v1
test_respects_explicit_valid_parser
test_invalid_parser_is_replaced_by_autodetection
test_reports_when_neither_parser_works
