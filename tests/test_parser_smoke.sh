#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${TEST_DIR}/.." && pwd -P)"
PIPELINE="${PACKAGE_DIR}/main.nf"

fail() {
  printf 'not ok - %s\n' "$*" >&2
  exit 1
}

resolve_nextflow() {
  if [[ -n "${NEXTFLOW_BIN:-}" ]]; then
    [[ -x "$NEXTFLOW_BIN" ]] || fail "NEXTFLOW_BIN is not executable: ${NEXTFLOW_BIN}"
    printf '%s\n' "$NEXTFLOW_BIN"
    return
  fi
  if [[ -x "${PACKAGE_DIR}/nextflow" ]]; then
    printf '%s\n' "${PACKAGE_DIR}/nextflow"
    return
  fi
  command -v nextflow 2>/dev/null || fail "Nextflow was not found"
}

NEXTFLOW="$(resolve_nextflow)"
command -v timeout >/dev/null 2>&1 || fail "timeout was not found"
[[ -f "$PIPELINE" ]] || fail "pipeline is missing: ${PIPELINE}"

VERSION_OUTPUT="$("$NEXTFLOW" -version 2>&1)" || fail "could not query Nextflow version"
NEXTFLOW_VERSION="$(
  printf '%s\n' "$VERSION_OUTPUT" |
    sed -nE 's/.*[Vv]ersion[[:space:]]+([^[:space:]]+).*/\1/p' |
    head -n 1
)"

ORIGINAL_NXF_HOME="${NXF_HOME:-${HOME}/.nextflow}"
CACHED_JAR=""
if [[ -n "$NEXTFLOW_VERSION" ]]; then
  jar_candidate="${ORIGINAL_NXF_HOME}/framework/${NEXTFLOW_VERSION}/nextflow-${NEXTFLOW_VERSION}-one.jar"
  [[ -s "$jar_candidate" ]] && CACHED_JAR="$jar_candidate"
fi

TEST_TMP="$(mktemp -d)"
trap 'rm -rf -- "$TEST_TMP"' EXIT

PLUGIN_DIR="${NXF_PLUGINS_DIR:-/opt/nextflow/plugins}"
if [[ ! -d "$PLUGIN_DIR" ]]; then
  PLUGIN_DIR="${TEST_TMP}/plugins"
  mkdir -p "$PLUGIN_DIR"
fi

run_parser() {
  local parser="$1"
  local parser_home="${TEST_TMP}/nxf-${parser}"
  local run_dir="${TEST_TMP}/run-${parser}"
  local work_dir="${TEST_TMP}/work-${parser}"
  local output="${TEST_TMP}/${parser}.log"
  local status=0
  local -a runtime_env=(
    "NXF_HOME=${parser_home}"
    "NXF_SYNTAX_PARSER=${parser}"
    "NXF_PLUGINS_DIR=${PLUGIN_DIR}"
    "NXF_OFFLINE=true"
    "NXF_DISABLE_CHECK_LATEST=true"
  )

  [[ -z "$NEXTFLOW_VERSION" ]] || runtime_env+=("NXF_VER=${NEXTFLOW_VERSION}")
  [[ -z "$CACHED_JAR" ]] || runtime_env+=("NXF_BIN=${CACHED_JAR}")
  mkdir -p "$parser_home" "$run_dir" "$work_dir"

  set +e
  (
    cd "$run_dir"
    timeout "${NFX_SMOKE_TIMEOUT:-90}" \
      env "${runtime_env[@]}" \
      "$NEXTFLOW" run "$PIPELINE" \
        -ansi-log false \
        -work-dir "$work_dir" \
        --help
  ) > "$output" 2>&1
  status=$?
  set -e

  if [[ "$status" -ne 0 ]]; then
    printf '%s\n' "--- Nextflow ${parser} smoke output ---" >&2
    sed -n '1,240p' "$output" >&2
    fail "Nextflow ${parser} --help smoke failed with status ${status}"
  fi
  [[ -s "$output" ]] || fail "Nextflow ${parser} --help produced no output"
  printf 'ok - real Nextflow parser %s compiles and renders --help\n' "$parser"
}

run_parser v1
run_parser v2
