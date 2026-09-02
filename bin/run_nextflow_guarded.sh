#!/usr/bin/env bash
set -Eeuo pipefail

# Guarded runner for this package. It starts the parent Nextflow command in its
# own process group. If Ctrl+C/TERM/HUP is received, the full process group,
# including any nested FFPErase workflow in --ffpe mode, is terminated.

OUTDIR=""
PROJECT_DIR=""
GRACE_SECONDS=25

usage() {
  cat <<'EOF'
run_nextflow_guarded.sh

Usage:
  bash bin/run_nextflow_guarded.sh --outdir OUTDIR --project-dir PACKAGE_DIR -- PACKAGE_DIR/bin/run_pipeline.sh [args...]

Options:
  --outdir DIR           Parent pipeline output directory.
  --project-dir DIR      Package directory. Default: inferred from this script.
  --grace-seconds N      TERM-to-KILL grace on abort. Default: 25.
  --                     End wrapper options; remaining arguments are the Nextflow command.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --project-dir) PROJECT_DIR="$2"; shift 2 ;;
    --grace-seconds|--grace) GRACE_SECONDS="$2"; shift 2 ;;
    --) shift; break ;;
    -h|--help) usage; exit 0 ;;
    *) break ;;
  esac
done

if (( $# == 0 )); then
  echo "ERROR: missing Nextflow command after --" >&2
  usage >&2
  exit 1
fi

[[ "$GRACE_SECONDS" =~ ^[0-9]+$ ]] || { echo "ERROR: --grace-seconds must be integer" >&2; exit 1; }

if [[ -z "$PROJECT_DIR" ]]; then
  PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
PROJECT_DIR="$(readlink -m "$PROJECT_DIR")"
if [[ -n "$OUTDIR" ]]; then OUTDIR="$(readlink -m "$OUTDIR")"; fi

PARENT_PID=""
PARENT_PGID=""
SELF_PGID="$(ps -o pgid= -p "$$" 2>/dev/null | tr -d '[:space:]' || true)"

kill_parent_group() {
  local sig="${1:-TERM}"
  if [[ -n "$PARENT_PID" && -z "$PARENT_PGID" ]]; then
    PARENT_PGID="$(ps -o pgid= -p "$PARENT_PID" 2>/dev/null | tr -d '[:space:]' || true)"
  fi
  if [[ -n "$PARENT_PGID" && "$PARENT_PGID" != "0" && "$PARENT_PGID" != "$SELF_PGID" ]]; then
    kill -"$sig" "-${PARENT_PGID}" 2>/dev/null || true
  elif [[ -n "$PARENT_PID" ]]; then
    pkill -"$sig" -P "$PARENT_PID" 2>/dev/null || true
    kill -"$sig" "$PARENT_PID" 2>/dev/null || true
  fi
}

on_signal() {
  local sig="$1"
  echo "[guarded-nextflow] received ${sig}; killing parent Nextflow process group" >&2
  kill_parent_group TERM
  sleep "$GRACE_SECONDS" || true
  kill_parent_group KILL
  exit 143
}

trap 'on_signal INT' INT
trap 'on_signal TERM' TERM
trap 'on_signal HUP' HUP

if command -v setsid >/dev/null 2>&1; then
  setsid "$@" &
else
  "$@" &
fi
PARENT_PID=$!
sleep 0.2 || true
PARENT_PGID="$(ps -o pgid= -p "$PARENT_PID" 2>/dev/null | tr -d '[:space:]' || true)"
echo "[guarded-nextflow] parent_pid=${PARENT_PID} parent_pgid=${PARENT_PGID:-unknown}"

set +e
wait "$PARENT_PID"
STATUS=$?
set -e

exit "$STATUS"
