#!/usr/bin/env bash
set -Eeuo pipefail

PACKAGE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

case "${1:-}" in
  fresh)
    exec "${PACKAGE_DIR}/examples/run_minimal_fresh.sh"
    ;;
  ffpe)
    exec "${PACKAGE_DIR}/examples/run_minimal_ffpe.sh"
    ;;
  ""|-h|--help)
    printf 'Usage: %s fresh|ffpe\n' "$0"
    printf 'Runs the bundled synthetic smoke example for the selected mode.\n'
    ;;
  *)
    printf 'Unknown example mode: %s\n' "$1" >&2
    printf 'Usage: %s fresh|ffpe\n' "$0" >&2
    exit 2
    ;;
esac
