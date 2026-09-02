# lowpass-variants-nextflow

A unified Nextflow workflow for low-pass per-sample SNV/indel calling from either FFPE or fresh-material BAMs. The workflow runs Mutect2, FreeBayes, and BCFtools; recalibrates BAM base qualities by default; left-aligns caller and final VCFs; and applies the external FFPErase workflow only in --ffpe mode.

## Entry point

Run the package through its launcher from the repository root:

~~~bash
./bin/run_pipeline.sh --help
~~~

The launcher probes the installed Nextflow with syntax parser v2 and then v1, selects the first parser that can compile the packaged workflow, and forwards all arguments unchanged. An explicitly valid NXF_SYNTAX_PARSER=v1 or NXF_SYNTAX_PARSER=v2 is respected.

Do not use nextflow run main.nf in normal usage because that bypasses parser auto-detection.

## Required sample mode

Every real run must select exactly one mode:

- --fresh: run the common low-pass calling and normalization workflow. Fresh mode performs no FFPE assessment: it does not calculate FFPE-specific thresholds or C>T/G>A risk, does not run FFPErase, and emits no FFPE-specific TSV columns or VCF INFO fields.
- --ffpe: run the same common workflow, apply the FFPE C>T/G>A heuristic during finalization, then run external FFPErase and annotate the final VCFs.

Supplying both modes, or neither mode, is an error.

## BQSR and known sites

Base Quality Score Recalibration is enabled by default in both modes. Supply one or more known-sites VCFs as a comma-separated value:

~~~text
--known_sites /refs/dbsnp.vcf.gz,/refs/Mills_and_1000G_gold_standard.indels.vcf.gz
~~~

Each known-sites VCF must match --ref and have an adjacent index. A compressed VCF may use .tbi, .csi, or .idx; an uncompressed VCF requires .idx. If BQSR is intentionally inappropriate for a dataset, use the explicit --skip_bqsr flag. A run without --known_sites is rejected unless --skip_bqsr is set.

## Input staging and safety

The workflow stages BAMs, the reference, known-sites VCFs, and optional Mutect2 germline-resource/panel-of-normals VCFs into task work directories. It copies --ref and creates its .fai and .dict there. If an input BAM needs an index, the index is also created from a private task copy; files beside the user-supplied BAM and reference are not modified.

An optional Mutect2 resource VCF must have an adjacent .tbi, .csi, or .idx. Resources staged into the same step must have distinct basenames; colliding VCF or index basenames are rejected.

## VCF normalization

Both modes use bcftools norm against --ref to split multiallelic records and left-align:

- Mutect2 VCFs
- FreeBayes VCFs
- BCFtools VCFs
- the merged final per-sample VCF

Normalization or indexing failure is treated as a workflow failure; an unnormalized VCF is not substituted as a successful result. Every normalized compressed VCF must support a .tbi index. References or coordinates outside TBI limits are currently unsupported.

## Fresh example

~~~bash
./bin/run_pipeline.sh \
  -profile conda \
  -resume \
  -work-dir /path/to/nextflow-work \
  --fresh \
  --input /path/to/bams \
  --outdir /path/to/fresh-results \
  --samples SAMPLE_01,SAMPLE_02 \
  --ref /path/to/reference.fa \
  --known_sites /path/to/dbsnp.vcf.gz,/path/to/known_indels.vcf.gz
~~~

## FFPE example

~~~bash
./bin/run_pipeline.sh \
  -profile conda \
  -resume \
  -work-dir /path/to/nextflow-work \
  --ffpe \
  --input /path/to/bams \
  --outdir /path/to/ffpe-results \
  --samples SAMPLE_01,SAMPLE_02 \
  --ref /path/to/reference.fa \
  --known_sites /path/to/dbsnp.vcf.gz,/path/to/known_indels.vcf.gz \
  --ffperase_threads 16 \
  --ffperase_classify_threads 2 \
  --ffperase_sample_jobs 1 \
  --ffperase_picard_jobs 1
~~~

For a deliberate no-BQSR run, replace --known_sites ... with --skip_bqsr.

## Existing outputs

Published pipeline products and logs are not replaced by default. Use --overwrite true to replace them deliberately. Nextflow execution metadata under <outdir>/nextflow_info is refreshed independently.

## FFPErase controls

These parameters are used only with --ffpe:

| Parameter | Purpose |
|---|---|
| --ffperase_threads | Parallelism for FFPErase preprocessing/pileup tasks. |
| --ffperase_classify_threads | Parallelism for FFPErase classification. |
| --ffperase_sample_jobs | Concurrent FFPErase sample/type jobs. |
| --ffperase_picard_jobs | Concurrent Picard artifact-metrics jobs. |
| --ffperase_split_pileup | Variants per FFPErase pileup split. |
| --ffperase_split_reads | Reads per FFPErase split. |
| --ffperase_container | FFPErase container path or URI. |
| --ffperase_engine | Container engine: auto, apptainer, singularity, or docker. |
| --ffperase_repository | Upstream nf-ffperase repository. |
| --ffperase_revision | Pinned upstream nf-ffperase revision. |
| --ffperase_filter_artifacts | Add FFPErase artifact FILTER tags to annotated VCF records. |
| --ffperase_fail_on_error | Fail the workflow when an FFPErase sample/type fails. |

FFPErase may need container and model downloads unless the requested container, revision, and models are already cached. For reproducibility, pin --ffperase_revision and use an immutable container reference where possible.

## Main outputs

Both modes emit independent, left-aligned per-sample VCFs:

~~~text
<outdir>/final_vcf/<sample>.final_variants.vcf
<outdir>/final_vcf/<sample>.final_variants.vcf.gz
<outdir>/final_vcf/<sample>.final_variants.vcf.gz.tbi
<outdir>/reports/per_sample/<sample>.final_variants.retained.tsv
~~~

FFPE mode additionally emits curated FFPErase classifications, status, Picard metric summaries, and annotated VCFs:

~~~text
<outdir>/ffperase_classification/classified/
<outdir>/ffperase_classification/final_vcf/<sample>.final_variants.ffperase_annotated.vcf.gz
<outdir>/ffperase_classification/final_vcf/<sample>.final_variants.ffperase_annotated.vcf.gz.tbi
<outdir>/ffperase_classification/reports/all_samples.ffperase_annotation_summary.tsv
<outdir>/ffperase_status/
<outdir>/ffperase_picard_metrics_summary/
<outdir>/logs/ffperase.done.txt
~~~

Transient FFPErase work trees, feature caches, models, staged BAMs, and reference copies remain in the Nextflow work directory for -resume and are not published as results.

## Portable configuration

Paths in commands and parameter files are examples. Replace them with paths available on the machine running Nextflow. The package does not assume a particular storage mount, home directory, Nextflow cache, container cache, or reference location.

The examples/ directory contains fresh and FFPE shell templates plus a parameter-file template.

## Validation

Run the static, unit, left-alignment, launcher, and real parser checks with:

~~~bash
bash tests/run_all.sh
~~~

An opt-in synthetic run also exercises real GATK BQSR and the complete fresh
workflow with callers skipped. If the tools are in a separate environment,
set `LOWPASS_TEST_TOOL_ENV` to that environment prefix:

~~~bash
RUN_FRESH_PIPELINE_SMOKE=true \
LOWPASS_TEST_TOOL_ENV=/path/to/tool-environment \
bash tests/run_all.sh
~~~

Before a production run, test one representative sample and review the Nextflow report, trace, normalization logs, caller-status files, and—when using --ffpe—FFPErase status tables and annotation summary.

## Third-party software and terms

This package invokes third-party tools, workflows, containers, models, and reference datasets. Their upstream licenses and terms apply independently. See [THIRD_PARTY.md](THIRD_PARTY.md) before redistribution or deployment. This repository does not add a project-wide license.
