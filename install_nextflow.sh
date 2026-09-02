#!/usr/bin/env bash
set -Eeuo pipefail

if command -v nextflow >/dev/null 2>&1; then
  echo "Nextflow already available: $(command -v nextflow)"
  nextflow -version || true
  exit 0
fi

if [[ ! -d "$HOME/bin" ]]; then
  mkdir -p "$HOME/bin"
fi

curl -s https://get.nextflow.io | bash
mv -f nextflow "$HOME/bin/nextflow"
chmod +x "$HOME/bin/nextflow"

if ! grep -q 'export PATH="$HOME/bin:$PATH"' "$HOME/.bashrc" 2>/dev/null; then
  echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
fi

export PATH="$HOME/bin:$PATH"
nextflow -version
