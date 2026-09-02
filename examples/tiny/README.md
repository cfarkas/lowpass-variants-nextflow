# Tiny synthetic example data

This directory contains a 200-base artificial reference, twelve artificial
reads supporting one SNV, and one artificial known site. It contains no human,
patient, or observed biological sequence data.

The committed BAM and indexes make the two smoke examples runnable immediately
after the software requirements are installed:

```bash
bash examples/run_minimal_fresh.sh
bash examples/run_minimal_ffpe.sh
```

Equivalent entry points are `make example-fresh`, `make example-ffpe`,
`./run_example.sh fresh`, and `./run_example.sh ffpe`. Results default to
`example-results/minimal-fresh/` or `example-results/minimal-ffpe/`. Override
that location with `OUTDIR=/another/path`. The scripts use the `conda` profile
by default; set `NEXTFLOW_PROFILE=` only when every required tool is already on
`PATH`.

These are installation and workflow-wiring checks, not variant-analysis
benchmarks. The fresh quickstart runs real BQSR, skips Mutect2, and calls with
FreeBayes and BCFtools. Its final VCF contains exactly `chr1:25 A>C`, with both
callers contributing (`VOTE_COUNT=2`). The FFPE example runs BQSR but skips all
callers, so its final VCF is header-only. It exercises the FFPE branch, Picard
artifact metrics, status handling, and VCF annotation; its SNV and indel
classifications are recorded as `SKIP_NO_VARIANTS`. It does not validate
non-empty FFPErase model inference.

Files:

- `tiny.fa`, `tiny.fa.fai`, and `tiny.dict`: artificial reference and sidecars.
- `SAMPLE.sam`: readable source for the artificial alignment.
- `SAMPLE.bam` and `SAMPLE.bam.bai`: sorted/indexed runnable alignment.
- `known-sites.vcf`, `known-sites.vcf.gz`, and
  `known-sites.vcf.gz.tbi`: artificial BQSR resource and indexed runnable copy.

To reproduce the generated files in a separate directory, use:

```bash
bash examples/tiny/rebuild.sh /tmp/lowpass-tiny-data
```

The rebuild command requires `samtools`, `bcftools`, `bgzip`, and `tabix` on
`PATH`.
