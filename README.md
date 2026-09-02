# lowpass-variants-nextflow

[![CI](https://github.com/cfarkas/lowpass-variants-nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/cfarkas/lowpass-variants-nextflow/actions/workflows/ci.yml)

A portable workflow for reproducible low-pass SNV/indel calling from fresh or FFPE BAMs.

## Getting Started

Clone the repository and run the bundled synthetic fresh example:

~~~bash
git clone https://github.com/cfarkas/lowpass-variants-nextflow.git
cd lowpass-variants-nextflow

bash examples/run_minimal_fresh.sh
~~~

The repository includes the artificial BAM, BAM index, 200-base reference,
reference sidecars, known-sites VCF, and VCF index under `examples/tiny/`.
The first run may take longer while Nextflow creates the supplied Conda
environment. This quickstart performs BQSR, skips Mutect2, and calls with
FreeBayes and BCFtools. Its main result is:

~~~text
example-results/minimal-fresh/final_vcf/SAMPLE.final_variants.vcf.gz
~~~

That VCF contains exactly one synthetic call, `chr1:25 A>C`, supported by both
enabled callers (`VOTE_COUNT=2`).

Run the bundled FFPE branch check with:

~~~bash
bash examples/run_minimal_ffpe.sh
~~~

The FFPE check also runs real BQSR, but deliberately skips all callers and
therefore produces a header-only final VCF. It reaches FFPE orchestration,
Picard metrics, empty-variant status handling, and annotation without performing
non-empty model inference. A host Apptainer, Singularity, or Docker executable
is still required. These examples verify installation and workflow wiring, not
biological calling accuracy. See [examples/tiny/README.md](examples/tiny/README.md).

Neither bundled command needs `-work-dir` or `-resume`; Nextflow uses `./work` in
the launch directory. Use the real-data commands below for an analysis.

## Table of Contents

- [Getting Started](#getting-started)
- [Users' Guide](#users-guide)
  - [Why this repository exists](#why-this-repository-exists)
  - [Installation](#installation)
  - [General usage](#general-usage)
  - [Mandatory parameters](#mandatory-parameters)
  - [Material modes](#material-modes)
  - [Minimal real-data runs](#minimal-real-data-runs)
  - [Complete run examples](#complete-run-examples)
  - [Main outputs](#main-outputs)
  - [Getting help](#getting-help)
- [Repository Guide](#repository-guide)
- [Tests](#tests)
- [Limitations](#limitations)
- [Third-party software and terms](#third-party-software-and-terms)

## Users' Guide

### Why this repository exists

Low-pass variant analysis is sensitive to inconsistent BAM preparation, caller-specific representations, and sample-material artifacts. This repository provides one auditable per-sample workflow that performs BQSR by default, combines evidence from Mutect2, FreeBayes, and BCFtools, and strictly normalizes and left-aligns caller and final VCFs.

The material mode is deliberately explicit. Fresh samples receive only the common calling workflow and no FFPE assessment. FFPE samples follow the same common path and then add the C>T/G>A review heuristic and external FFPErase classification. Requiring exactly one mode prevents FFPE assumptions from being silently applied to fresh material.

The package also stages inputs without modifying source BAM/reference locations and supplies a launcher that selects a compatible Nextflow syntax parser (v2 first, then v1).

### Installation

Requirements:

- Linux with Java 17 or newer.
- Nextflow 24.10 or newer.
- Conda or Mamba when using the supplied `conda` profile.
- Docker with access to a local daemon when using the supplied `docker` profile.
- Python 3 on the host for input discovery when using the `docker` profile.
- For FFPErase: Apptainer, Singularity, or Docker, plus access to the upstream workflow, container, and models.

Install from the private GitHub repository:

~~~bash
git clone https://github.com/cfarkas/lowpass-variants-nextflow.git
cd lowpass-variants-nextflow
~~~

If Nextflow is not already installed, the included helper installs it under `$HOME/bin`:

~~~bash
bash install_nextflow.sh
~~~

Verify the package and parser-selection launcher:

~~~bash
./bin/run_pipeline.sh --help
./bin/run_pipeline.sh --help_full
~~~

The recommended execution profile is `-profile conda`, which lets Nextflow create and reuse the declared environments in `envs/`. If all required tools are supplied by another environment or site profile, the Conda profile may be omitted.

The alternative `-profile docker` runs pipeline tasks in the versioned image
`ghcr.io/cfarkas/lowpass-variants-nextflow:1.0.0`. Pull it with:

~~~bash
docker pull ghcr.io/cfarkas/lowpass-variants-nextflow:1.0.0
~~~

Before the registry tag is available, build the same tag from this checkout:

~~~bash
docker build -t ghcr.io/cfarkas/lowpass-variants-nextflow:1.0.0 .
~~~

The launcher and Nextflow still run on the host; the image supplies the tools
used by pipeline tasks. Docker will pull the configured image automatically if
it is not already local.

Run both bundled checks through Docker with:

~~~bash
NEXTFLOW_PROFILE=docker bash examples/run_minimal_fresh.sh
NEXTFLOW_PROFILE=docker bash examples/run_minimal_ffpe.sh
~~~

### General usage

Run the package through its launcher from the repository root:

~~~bash
./bin/run_pipeline.sh --help
~~~

The launcher probes the installed Nextflow with syntax parser v2 and then v1, selects the first parser that can compile the packaged workflow, and forwards all arguments unchanged. An explicitly valid NXF_SYNTAX_PARSER=v1 or NXF_SYNTAX_PARSER=v2 is respected.

Do not use nextflow run main.nf in normal usage because that bypasses parser auto-detection.

### Mandatory parameters

Every analysis run must provide:

| Parameter | Requirement |
|---|---|
| `--fresh` or `--ffpe` | Choose exactly one. Both or neither is an error. |
| `--input PATH` | BAM directory, one BAM, a manifest, or a comma-separated BAM list. |
| `--outdir PATH` | Destination for published results and execution reports. |
| `--ref FASTA` | Reference used for BAM processing, calling, and normalization. |
| `--known_sites VCF1,VCF2` or `--skip_bqsr` | Supply indexed known-sites VCFs for default BQSR, or deliberately disable BQSR. |

`--samples` is optional. Without it, every BAM resolved from `--input` is selected.

`--dry_run true` is an input-resolution check: it still requires the mode,
input, output, and reference parameters above, but it does not start the
analysis and therefore does not require known sites or `--skip_bqsr`.

### Optional Nextflow run controls

These are Nextflow options, not mandatory pipeline parameters:

- `-work-dir DIR` places task state in a selected directory. Without it, Nextflow uses `./work` in the launch directory.
- `-resume` reuses completed tasks from the same work directory after an interruption or rerun.
- `-profile conda` selects the supplied dependency environments. It is recommended unless tools are provided another way.
- `-profile docker` runs tasks in `ghcr.io/cfarkas/lowpass-variants-nextflow:1.0.0` and requires a local Docker daemon.

For large BAM workflows, use a persistent work directory with adequate free space and keep it unchanged when resuming.

### Material modes

Every real run must select exactly one mode:

- --fresh: run the common low-pass calling and normalization workflow. Fresh mode performs no FFPE assessment: it does not calculate FFPE-specific thresholds or C>T/G>A risk, does not run FFPErase, and emits no FFPE-specific TSV columns or VCF INFO fields.
- --ffpe: run the same common workflow, apply the FFPE C>T/G>A heuristic during finalization, then run external FFPErase and annotate the final VCFs.

Supplying both modes, or neither mode, is an error.

### BQSR and known sites

Base Quality Score Recalibration is enabled by default in both modes. Supply one or more known-sites VCFs as a comma-separated value:

~~~text
--known_sites /refs/dbsnp.vcf.gz,/refs/Mills_and_1000G_gold_standard.indels.vcf.gz
~~~

Each known-sites VCF must match --ref and have an adjacent index. A compressed VCF may use .tbi, .csi, or .idx; an uncompressed VCF requires .idx. If BQSR is intentionally inappropriate for a dataset, use the explicit --skip_bqsr flag. A run without --known_sites is rejected unless --skip_bqsr is set.

### Input staging and safety

The workflow stages BAMs, the reference, known-sites VCFs, and optional Mutect2 germline-resource/panel-of-normals VCFs into task work directories. It copies --ref and creates its .fai and .dict there. If an input BAM needs an index, the index is also created from a private task copy; files beside the user-supplied BAM and reference are not modified.

An optional Mutect2 resource VCF must have an adjacent .tbi, .csi, or .idx. Resources staged into the same step must have distinct basenames; colliding VCF or index basenames are rejected.

### VCF normalization

Both modes use bcftools norm against --ref to split multiallelic records and left-align:

- Mutect2 VCFs
- FreeBayes VCFs
- BCFtools VCFs
- the merged final per-sample VCF

Normalization or indexing failure is treated as a workflow failure; an unnormalized VCF is not substituted as a successful result. Every normalized compressed VCF must support a .tbi index. References or coordinates outside TBI limits are currently unsupported.

### Minimal real-data runs

These commands contain only the required pipeline inputs plus the recommended
Conda profile. `--samples`, `-work-dir`, `-resume`, Mutect2 resources, and
resource-tuning flags are optional.

#### Minimal fresh analysis

~~~bash
./bin/run_pipeline.sh \
  -profile conda \
  --fresh \
  --input /data/SAMPLE_01.bam \
  --outdir /data/results/fresh \
  --ref /refs/reference.fa \
  --known_sites /refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz
~~~

#### Minimal FFPE analysis

~~~bash
./bin/run_pipeline.sh \
  -profile conda \
  --ffpe \
  --input /data/SAMPLE_01.bam \
  --outdir /data/results/ffpe \
  --ref /refs/reference.fa \
  --known_sites /refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz
~~~

For either mode, replace `--known_sites ...` with `--skip_bqsr` only when BQSR
is deliberately inappropriate. A non-empty FFPE analysis may download the
pinned upstream workflow, container image, and models unless the same task is
resumed from populated caches.

#### Minimal fresh analysis with Docker

~~~bash
./bin/run_pipeline.sh \
  -profile docker \
  --fresh \
  --input /data/SAMPLE_01.bam \
  --outdir /data/results/fresh \
  --ref /refs/reference.fa \
  --known_sites /refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz
~~~

#### Minimal FFPE analysis with Docker

~~~bash
./bin/run_pipeline.sh \
  -profile docker \
  --ffpe \
  --input /data/SAMPLE_01.bam \
  --outdir /data/results/ffpe \
  --ref /refs/reference.fa \
  --known_sites /refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz
~~~

Use absolute host paths for Docker runs. In FFPE mode, the profile mounts the
local Docker daemon socket into the FFPErase orchestration task so it can launch
the separate upstream runtime image. Socket access is privileged-equivalent
host access; use this profile only with trusted workflow code and images. The
main pipeline image does not bundle the upstream FFPErase workflow, its runtime
image, or its models. A non-empty FFPE run fetches those assets separately
unless they are already cached.

### Complete run examples

The first two examples show resumable production commands. `-work-dir` and
`-resume` are optional Nextflow controls, not mandatory pipeline flags. Remove
both lines to use Nextflow's default `./work` directory and start without
resume. `-profile conda` may also be omitted when another execution environment
provides all required tools.

#### Fresh run with BQSR and all three callers

~~~bash
./bin/run_pipeline.sh \
  -profile conda \
  -work-dir /scratch/project/fresh-nextflow-work \
  -resume \
  --fresh \
  --input /data/bams \
  --outdir /data/results/fresh \
  --samples SAMPLE_01,SAMPLE_02 \
  --ref /refs/reference.fa \
  --known_sites /refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz \
  --duplicate_mode mark \
  --germline_resource /refs/af-only-gnomad.vcf.gz \
  --panel_of_normals /refs/panel-of-normals.vcf.gz \
  --prepare_cpus 8 \
  --mutect2_cpus 8 \
  --freebayes_cpus 2 \
  --bcftools_cpus 8 \
  --normalize_cpus 4 \
  --finalize_cpus 4 \
  --max_samples_parallel 4
~~~

`--samples`, the two Mutect2 resource VCFs, and the resource-tuning flags are
optional. Each supplied resource VCF needs an adjacent index. Fresh mode stops
after the common calling and normalization path; none of the FFPE thresholds,
heuristics, jobs, or output fields are used.

#### FFPE run with BQSR and FFPErase

~~~bash
./bin/run_pipeline.sh \
  -profile conda \
  -work-dir /scratch/project/ffpe-nextflow-work \
  -resume \
  --ffpe \
  --input /data/bams \
  --outdir /data/results/ffpe \
  --samples SAMPLE_01,SAMPLE_02 \
  --ref /refs/reference.fa \
  --known_sites /refs/dbsnp.vcf.gz,/refs/known_indels.vcf.gz \
  --duplicate_mode mark \
  --germline_resource /refs/af-only-gnomad.vcf.gz \
  --panel_of_normals /refs/panel-of-normals.vcf.gz \
  --prepare_cpus 8 \
  --mutect2_cpus 8 \
  --freebayes_cpus 2 \
  --bcftools_cpus 8 \
  --normalize_cpus 4 \
  --finalize_cpus 4 \
  --max_samples_parallel 4 \
  --ffperase_threads 16 \
  --ffperase_classify_threads 2 \
  --ffperase_sample_jobs 1 \
  --ffperase_picard_jobs 1 \
  --ffperase_split_pileup 5000 \
  --ffperase_split_reads 7500000 \
  --ffperase_container /containers/ffperase.sif \
  --ffperase_engine apptainer \
  --ffperase_repository papaemmelab/nf-ffperase \
  --ffperase_revision b0dd56cbd0a939896a966b9ce30c4d719b158170 \
  --ffperase_filter_artifacts true \
  --ffperase_fail_on_error true
~~~

Replace the example container with a local file or immutable URI available on
the execution host. The upstream revision shown is the package default.

#### Fresh run that deliberately skips BQSR

~~~bash
./bin/run_pipeline.sh \
  -profile conda \
  --fresh \
  --input /data/example/SAMPLE_01.bam \
  --outdir ./no-bqsr-results \
  --ref /refs/reference.fa \
  --skip_bqsr
~~~

Use the no-BQSR form only when recalibration is inappropriate for the data.

### Existing outputs

Published pipeline products and logs are not replaced by default. Use --overwrite true to replace them deliberately. Nextflow execution metadata under <outdir>/nextflow_info is refreshed independently.

### FFPErase controls

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

### Main outputs

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

### Portable configuration

Paths in commands and parameter files are examples. Replace them with paths available on the machine running Nextflow. The package does not assume a particular storage mount, home directory, Nextflow cache, container cache, or reference location.

The examples/ directory contains fresh and FFPE shell templates plus a parameter-file template.

### Getting help

Use the short help for mandatory inputs and the quick example, or full help for
all documented parameters, FFPErase controls, and complete commands:

~~~bash
./bin/run_pipeline.sh --help
./bin/run_pipeline.sh --help_full
make help
~~~

Operational answers, including recovery after an interrupted download, are in
[FAQ.md](FAQ.md).

## Repository Guide

| Path | Purpose |
|---|---|
| `README.md` | Installation, quickstart, complete usage, outputs, and limits. |
| `FAQ.md` | Operational questions and recovery guidance. |
| `NEWS.md` | User-visible release notes. |
| `Makefile` | Stable `help`, `lint`, `test`, and `install-nextflow` commands. |
| `main.nf`, `nextflow.config` | Workflow graph, defaults, profiles, and reports. |
| `Dockerfile`, `docker/` | Versioned pipeline-task image definition and container runtime dependencies. |
| `bin/` | Launcher, input resolver, finalization, and FFPErase helpers. |
| `envs/`, `environment.yml` | Reproducible Conda environments. |
| `examples/` | Ready-to-run synthetic smoke data plus production templates. |
| `assets/` | Short/full CLI help and packaged static inputs. |
| `tests/` | Launcher, unit, parser, normalization, and workflow smoke tests. |
| `FILES.txt` | Auditable package manifest. |
| `THIRD_PARTY.md` | Upstream software and redistribution terms. |

## Tests

Run the static, unit, left-alignment, launcher, and real parser checks with:

~~~bash
bash tests/run_all.sh
~~~

The normal suite validates the committed synthetic BAM/reference/resources and
checks that both minimal launchers forward the correct mode and inputs.

An opt-in synthetic run also exercises real GATK BQSR and the fresh workflow
with FreeBayes and BCFtools enabled. It asserts the single `chr1:25 A>C` call
and two-caller support. If the tools are in a separate environment, set
`LOWPASS_TEST_TOOL_ENV` to that environment prefix:

~~~bash
RUN_FRESH_PIPELINE_SMOKE=true \
LOWPASS_TEST_TOOL_ENV=/path/to/tool-environment \
bash tests/run_all.sh
~~~

The bundled FFPE empty-variant path has a separate opt-in end-to-end check:

~~~bash
RUN_FFPE_PIPELINE_SMOKE=true \
LOWPASS_TEST_TOOL_ENV=/path/to/tool-environment \
bash tests/run_all.sh
~~~

Before a production run, test one representative sample and review the Nextflow report, trace, normalization logs, caller-status files, and—when using --ffpe—FFPErase status tables and annotation summary.

The same normal suite runs in GitHub Actions on every push and pull request.

## Limitations

- The supplied profiles target Linux with local or Slurm execution. Other
  schedulers require a site-specific Nextflow profile.
- Final compressed VCFs require TBI indexes; references or coordinates beyond
  TBI limits are not supported.
- FFPErase is an external workflow with its own container, models, network/cache
  requirements, and upstream terms. Its availability is independent of the
  common fresh-material path.
- The bundled FFPE example has no variants and therefore does not test model
  classification. Validate at least one representative non-empty FFPE sample
  before production use.
- Low-pass data can have limited sensitivity and unstable allele fractions.
  Validate thresholds and calls for the intended assay; this workflow is not a
  substitute for clinical validation.

## Third-party software and terms

This package invokes third-party tools, workflows, containers, models, and reference datasets. Their upstream licenses and terms apply independently. See [THIRD_PARTY.md](THIRD_PARTY.md) before redistribution or deployment. This repository does not add a project-wide license.
