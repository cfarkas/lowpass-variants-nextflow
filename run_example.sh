#!/usr/bin/env bash
set -Eeuo pipefail
PACKAGE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
exec "${PACKAGE_DIR}/bin/run_pipeline.sh" --help
