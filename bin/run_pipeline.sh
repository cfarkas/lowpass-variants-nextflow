#!/usr/bin/env bash
set -Eeuo pipefail

# Parser-compatible entry point for the packaged Nextflow workflow.
# NEXTFLOW_BIN may name an executable on PATH or an explicit executable path.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd -P)"
PIPELINE="${PACKAGE_DIR}/main.nf"

die() {
  printf '[lowpass-variants] ERROR: %s\n' "$*" >&2
  exit 1
}

resolve_command() {
  local requested="$1"
  local resolved=""

  if [[ "$requested" == */* ]]; then
    [[ -f "$requested" && -x "$requested" ]] || return 1
    printf '%s\n' "$requested"
    return 0
  fi

  resolved="$(command -v -- "$requested" 2>/dev/null || true)"
  [[ -n "$resolved" && -x "$resolved" ]] || return 1
  printf '%s\n' "$resolved"
}

find_nextflow() {
  local resolved=""

  if [[ -n "${NEXTFLOW_BIN:-}" ]]; then
    resolved="$(resolve_command "$NEXTFLOW_BIN" || true)"
    [[ -n "$resolved" ]] || die "NEXTFLOW_BIN is not executable or was not found: ${NEXTFLOW_BIN}"
    printf '%s\n' "$resolved"
    return 0
  fi

  # Prefer a package-local executable when present, then the user's PATH.
  if [[ -x "${PACKAGE_DIR}/nextflow" ]]; then
    printf '%s\n' "${PACKAGE_DIR}/nextflow"
    return 0
  fi

  resolved="$(resolve_command nextflow || true)"
  [[ -n "$resolved" ]] || die "Nextflow was not found; install it or set NEXTFLOW_BIN"
  printf '%s\n' "$resolved"
}

probe_parser() {
  local parser="$1"
  NXF_SYNTAX_PARSER="$parser" "$NEXTFLOW" run "$PIPELINE" --help \
    >/dev/null 2>&1
}

[[ -f "$PIPELINE" ]] || die "pipeline entry point is missing: ${PIPELINE}"
NEXTFLOW="$(find_nextflow)"

case "${NXF_SYNTAX_PARSER:-}" in
  v1|v2)
    SELECTED_PARSER="$NXF_SYNTAX_PARSER"
    PARSER_SOURCE="environment"
    ;;
  *)
    if [[ -n "${NXF_SYNTAX_PARSER:-}" ]]; then
      printf '[lowpass-variants] WARNING: ignoring invalid NXF_SYNTAX_PARSER=%q; expected v1 or v2\n' \
        "$NXF_SYNTAX_PARSER" >&2
    fi
    unset NXF_SYNTAX_PARSER
    SELECTED_PARSER=""
    PARSER_SOURCE="autodetected"
    for candidate in v2 v1; do
      printf '[lowpass-variants] Probing Nextflow syntax parser %s\n' "$candidate" >&2
      if probe_parser "$candidate"; then
        SELECTED_PARSER="$candidate"
        break
      fi
    done
    [[ -n "$SELECTED_PARSER" ]] || die \
      "the packaged pipeline did not pass --help with either syntax parser (tried v2, then v1)"
    ;;
esac

VERSION_OUTPUT="$("$NEXTFLOW" -version 2>&1 || true)"
VERSION_SUMMARY="$(
  printf '%s\n' "$VERSION_OUTPUT" |
    awk 'NF { gsub(/^[[:space:]]+|[[:space:]]+$/, ""); printf "%s%s", separator, $0; separator=" " } END { print "" }'
)"
[[ -n "$VERSION_SUMMARY" ]] || VERSION_SUMMARY="version unavailable"

printf '[lowpass-variants] Nextflow: %s (%s)\n' "$NEXTFLOW" "$VERSION_SUMMARY" >&2
printf '[lowpass-variants] Syntax parser: %s (%s)\n' "$SELECTED_PARSER" "$PARSER_SOURCE" >&2

export NXF_SYNTAX_PARSER="$SELECTED_PARSER"
exec "$NEXTFLOW" run "$PIPELINE" "$@"
