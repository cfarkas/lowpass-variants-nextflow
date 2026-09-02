# Release notes

## 1.0.0 (unreleased)

- Unified the fresh-material and FFPE low-pass workflows behind
  `./bin/run_pipeline.sh` with exactly one material mode required per run.
- Enabled BQSR by default with comma-separated indexed `--known_sites` inputs
  and an explicit `--skip_bqsr` escape hatch.
- Added strict `bcftools norm` splitting and left alignment for every caller VCF
  and each final per-sample VCF.
- Isolated FFPE-only behavior: fresh outputs contain no FFPE assessment fields,
  while FFPE runs add the C>T/G>A review heuristic and external FFPErase
  classification.
- Added Nextflow syntax parser v2/v1 selection, portable examples, full help,
  third-party notices, and guarded publication of existing outputs.
- Added CI plus synthetic tests for launcher behavior, sample resolution,
  material-mode separation, package rules, VCF normalization, and both Nextflow
  parsers. An optional smoke test exercises a real fresh run with BQSR.
- Added a committed artificial BAM/reference/known-sites fixture and one-command
  minimal fresh and empty-variant FFPE smoke examples.

The FFPErase integration remains subject to its upstream terms. See
[THIRD_PARTY.md](THIRD_PARTY.md).
