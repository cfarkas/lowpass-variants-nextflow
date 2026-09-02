SHELL := /bin/bash
PYTHON ?= python3
NEXTFLOW ?= nextflow

.PHONY: help lint test install-nextflow

help:
	@printf '%s\n' \
	  'Available targets:' \
	  '  help              Show this target list.' \
	  '  lint              Check Bash, Python, and Nextflow syntax.' \
	  '  test              Run the normal aggregate test suite.' \
	  '  install-nextflow  Run the included Nextflow installer.'

lint:
	@find bin examples tests -type f -name '*.sh' -exec bash -n {} +
	@bash -n install_nextflow.sh run_example.sh tests/fixtures/nextflow
	@$(PYTHON) -c 'from pathlib import Path; paths = sorted(Path("bin").glob("*.py")); assert paths, "no Python files found"; [compile(path.read_bytes(), str(path), "exec") for path in paths]'
	@command -v "$(NEXTFLOW)" >/dev/null 2>&1 || { printf 'Nextflow not found; run make install-nextflow\n' >&2; exit 127; }
	@NXF_SYNTAX_PARSER=v2 "$(NEXTFLOW)" lint .

test:
	@./tests/run_all.sh

install-nextflow:
	@./install_nextflow.sh
