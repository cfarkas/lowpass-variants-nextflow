# Frequently asked questions

## Which command should I use to run the workflow?

Run `./bin/run_pipeline.sh` from the repository root. The launcher probes the
installed Nextflow with syntax parser v2 and then v1, selects a working parser,
and forwards the original arguments. Calling `nextflow run main.nf` directly
bypasses that selection.

## Are runnable example data included?

Yes. `examples/tiny/` contains a synthetic BAM and index, reference and
sidecars, and an indexed known-sites VCF. After installing Nextflow and Conda,
run:

```bash
bash examples/run_minimal_fresh.sh
bash examples/run_minimal_ffpe.sh
```

The fresh quickstart runs BQSR, skips Mutect2, and calls with FreeBayes and
BCFtools. Its final VCF contains exactly `chr1:25 A>C` with two supporting
callers. The FFPE check runs BQSR but skips all callers, producing a header-only
final VCF and verifying empty-variant FFPE handling. These are installation
checks, not biological benchmarks; meaningful FFPErase classification still
requires real variants and the external workflow, image, and models.

## Why does the workflow ask for `--fresh` or `--ffpe`?

Every non-help run must select exactly one material mode. Supplying both flags,
or neither flag, is an error.

- `--fresh` runs the common preparation, calling, merging, and normalization
  path. It performs no FFPE assessment and emits no FFPE-specific fields.
- `--ffpe` adds the C>T/G>A review heuristic, runs the external FFPErase
  workflow, and annotates the final VCFs with its classifications.

## Why are known-sites VCFs required?

BQSR is enabled by default. Pass one or more assembly-compatible, indexed VCFs
as a comma-separated `--known_sites` value. Use `--skip_bqsr` only when you
deliberately want to bypass recalibration. Compressed known-sites VCFs may use
an adjacent `.tbi`, `.csi`, or `.idx`; uncompressed VCFs require `.idx`.

## Is `-work-dir` mandatory?

No. Without it, Nextflow uses `./work` in the directory where the command is
launched. Set `-work-dir` when you want task files on a larger or persistent
filesystem. It is especially useful with `-resume`, but neither option is
required for a new run.

## How do I use the Docker profile?

The profile uses the public image
`ghcr.io/cfarkas/lowpass-variants-nextflow:1.0.0`. Pull it anonymously with:

```bash
docker pull ghcr.io/cfarkas/lowpass-variants-nextflow:1.0.0
```

Then use `-profile docker` in place of `-profile conda`. Docker also pulls the
configured image on demand. To build the same tag locally instead:

```bash
docker build -t ghcr.io/cfarkas/lowpass-variants-nextflow:1.0.0 .
```

The lightweight input-discovery step runs before BAMs are staged and requires
`python3` on the host. Run the bundled checks in Docker mode with:

```bash
NEXTFLOW_PROFILE=docker bash examples/run_minimal_fresh.sh
NEXTFLOW_PROFILE=docker bash examples/run_minimal_ffpe.sh
```

Fresh runs require normal access to the local Docker daemon. For FFPE runs, the
profile additionally mounts `/var/run/docker.sock` into the FFPErase task so
the external workflow can launch its own image. Access to that socket is
privileged-equivalent host access, so run only trusted workflow code and
images. The upstream FFPErase workflow, runtime image, and models are not in
the main image and are fetched separately when needed unless already cached.

## Will the workflow modify my input BAM or reference directory?

The workflow stages input data in Nextflow task directories. If a BAM or
reference index is needed, it is created from a private staged copy rather than
beside the user-supplied input.

## How do I resume after an interruption?

Repeat the same command and add `-resume`. If you set `-work-dir`, keep the same
path; if you used Nextflow's default, launch from the same directory so the same
`./work` is available. Do not remove that work directory before resuming.
Changes to inputs, parameters, or workflow code can cause affected tasks to run
again.

## What should I check after an FFPErase network failure?

FFPErase can fetch its pinned upstream workflow, a container, and model files on
first use. Confirm network access, free space, the selected container engine,
and any configured cache paths, then resume the original Nextflow command with
`-resume`. Pin `--ffperase_revision` and use an immutable container reference
when reproducible external artifacts are required.

## Where are the main results?

The common per-sample results are under `<outdir>/final_vcf/` and
`<outdir>/reports/per_sample/`. FFPE mode additionally publishes status under
`<outdir>/ffperase_status/`, Picard summaries under
`<outdir>/ffperase_picard_metrics_summary/`, and classifications plus annotated
VCFs below `<outdir>/ffperase_classification/`.

## Why was an existing output not replaced?

Published products are protected by default. Pass `--overwrite true` only when
you intend to replace files in the selected output directory.

## How can I check the package before a production run?

Run `make lint` for shell, Python, and Nextflow checks, followed by `make test`
for the normal test suite. The real fresh/BQSR workflow smoke test is optional:

```bash
RUN_FRESH_PIPELINE_SMOKE=true make test
```

That optional test skips itself when its local tool requirements are missing.
The normal suite does not launch FFPErase or fetch its container or models.

## What does `make install-nextflow` change?

It runs the included `install_nextflow.sh` helper. If Nextflow is absent, the
helper downloads it to `$HOME/bin` and may add that directory to `$HOME/.bashrc`.

## What terms apply to FFPErase and this repository?

Read [THIRD_PARTY.md](THIRD_PARTY.md) before using or redistributing the package.
FFPErase has restrictive upstream terms, and this repository does not add a
project-wide license.
