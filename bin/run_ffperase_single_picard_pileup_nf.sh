#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# run_ffperase_single_picard_pileup_nf.sh
# Version: 2026-06-01-v3-sanitize-ambiguous-snvs-for-ffperase
#
# Single-file FFPErase launcher for already generated final VCFs/BAMs.
#
# It does two logical phases:
#   1) PREPROCESS
#      1a. Precompute Picard/GATK CollectSequencingArtifactMetrics once per sample
#          outside nf-ffperase.
#      1b. Run nf-ffperase --step preprocess using --picardMetrics to compute
#          pileup/features.tsv, avoiding internal PICARD scatter.
#   2) CLASSIFY
#      Run nf-ffperase --step classify using the cached features.tsv produced in
#      phase 1.
#
# This file is self-contained only with respect to the local Picard precompute
# and orchestration helpers embedded below. It does not vendor or redistribute
# papaemmelab/nf-ffperase source. Nextflow retrieves and runs the pinned
# upstream revision at runtime. See THIRD_PARTY.md for its terms of use.
###############################################################################

SCRIPT_NAME="$(basename "$0")"
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"

ROOT=""
REF=""
SAMPLES=""
OUTROOT=""
PICARD_ROOT=""
FEATURES_ROOT=""
MODELS_DIR=""
TOOLS_DIR=""

THREADS=16
CLASSIFY_THREADS=2
SAMPLE_JOBS=1
PICARD_JOBS=1
SPLIT_PILEUP=5000
SPLIT_READS=7500000
JAVA_MEM="12g"
TIMEOUT_MIN=7200
PICARD_TIMEOUT_MIN=720
CLASSIFY_TIMEOUT_MIN=720
MIN_MAPQ=20
MIN_BASEQ=20
MIN_DEPTH=1

RUN_PICARD=true
RUN_PREPROCESS=true
RUN_CLASSIFY=true
FORCE_PICARD=false
FORCE_PREPROCESS=false
FORCE_CLASSIFY=false
FAIL_ON_ERROR=true
RECOVER_MISSING_FEATURES=true
UPPERCASE_REFERENCE=true
SANITIZE_SNVS=true
USE_RESUME=false
DOWNLOAD_MODELS=true
PREFER_LOCAL_MODEL=true

RUN_SNVS=true
RUN_INDELS=true

CONTAINER="auto"
ENGINE="auto"
REPOSITORY="papaemmelab/nf-ffperase"
REVISION="b0dd56cbd0a939896a966b9ce30c4d719b158170"

PICARD_BACKEND="gatk"
PICARD_CONDA_ENV="ffperase_picard_metrics"
PICARD_PROJECT_DIR=""
PICARD_GATK_BIN=""
PICARD_JAR=""
PICARD_NO_CONDA=false
PICARD_NO_CREATE_CONDA=false

LOG=""

usage() {
cat <<'EOF_HELP'
run_ffperase_single_picard_pileup_nf.sh

Purpose
-------
Single-file launcher to preprocess Picard + pileup/features, then run
nf-ffperase classification using the precomputed inputs.

Input root expected by default
------------------------------
  <root>/final_vcf/<sample>.final_variants.vcf.gz
  <root>/preprocessed_bam/<sample>.preprocessed.bam
  <root>/preprocessed_bam/<sample>.preprocessed.bam.bai

What it does
------------
  Phase 1 / preprocess:
    A. Picard/GATK CollectSequencingArtifactMetrics once per sample.
       Output:
         <picard-metrics-root>/<sample>/pre_adapter_metrics.tsv
         <picard-metrics-root>/<sample>/bait_bias_metrics.tsv

    B. nf-ffperase --step preprocess using --picardMetrics.
       This computes pileup + features.tsv but avoids nf-ffperase internal
       PICARD scatter.
       Output:
         <features-root>/<sample>/snvs/features.tsv
         <features-root>/<sample>/indels/features.tsv

  Phase 2 / classify:
    C. nf-ffperase --step classify using cached features.tsv from Phase 1.
       Output:
         <outroot>/<sample>/snvs/classify/classified_df_snvs.tsv
         <outroot>/<sample>/indels/classify/classified_df_indels.tsv
         <outroot>/classified/all_samples.ffperase_snvs.tsv
         <outroot>/classified/all_samples.ffperase_indels.tsv

Third-party boundary
--------------------
  This launcher contains local adapter code, not the upstream nf-ffperase
  workflow source. Nextflow obtains the pinned upstream revision at runtime.
  FFPErase and its outputs are subject to upstream terms; see THIRD_PARTY.md.

Required
--------
  --root PATH
      Variant-calling root with final_vcf/ and preprocessed_bam/.

  --ref FASTA
      Reference FASTA.

Recommended
-----------
  --samples CSV
      Comma-separated sample IDs. If omitted, samples are inferred from files.

Main output paths
-----------------
  --outroot PATH
      Default: <root>/ffperase_classification

  --picard-metrics-root PATH
      Default: <root>/ffperase_picard_metrics

  --features-root PATH
      Default: <outroot>/features_cache

Threading / parallelism
-----------------------
  --threads N
      Controls nested nf-ffperase maxForks during preprocess-only. This is the
      useful knob for PILEUP chunk parallelism. Default: 16.

  --classify-threads N
      Controls nested nf-ffperase maxForks during classify-only. Default: 2.

  --sample-jobs N, --jobs N
      Number of sample/type nf-ffperase jobs launched at the same time.
      Default: 1. Keep 1 for WGS unless disk I/O is very strong.

  --picard-jobs N
      Number of samples processed in parallel during Picard precompute.
      Default: 1. Picard itself is not effectively multithreaded per BAM.

  --split-pileup N
      Number of variants per nf-ffperase pileup split. Default: 5000.
      Try 10000 if task overhead dominates and memory/I/O behave well.

Stage control
-------------
  --skip-picard
      Do not precompute Picard; use existing --picard-metrics-root.

  --skip-preprocess
      Do not run nf-ffperase preprocess/pileup; use existing --features-root.

      In normal full mode, preprocess automatically skips sample/type pairs
      where <features-root>/<sample>/<type>/features.tsv already exists, unless
      --force-preprocess or --force is used.

  --skip-classify
      Do not run classification.

  --picard-only
      Run only Picard precompute.

  --preprocess-only
      Run Picard if enabled, then nf-ffperase preprocess/pileup/features only.

  --classify-only
      Run only classification from existing --features-root.

Force / resume
--------------
  --force
      Force all enabled stages.

  --force-picard
  --force-preprocess
  --force-classify

  --use-resume true|false
      Default: false for nested nf-ffperase to avoid stale locks.

  --fail-on-error true|false
      Default: true. This prevents silent partial runs such as indels-only
      classification when SNV features failed during preprocessing.

  --no-recover-missing-features
      Disable the automatic one-pass retry for missing features.tsv files
      after Phase 1B. Default: retry missing SNV/indel features once.

  --sanitize-snvs true|false
      Default: true. Before nf-ffperase sees an SNV VCF, remove records that
      cannot be indexed by FFPErase/Picard context tables: non-ACGT REF/ALT,
      multiallelic SNVs left after splitting, or SNVs whose +/-1 reference
      context contains N/non-ACGT. Skipped records are written to
      <outroot>/<sample>/input/<sample>.snvs.skipped_for_ffperase.tsv.

Container / nf-ffperase
-----------------------
  --container PATH_OR_URI
      Default: auto. Auto uses $APPTAINER_CACHEDIR/nf-ffperase_v1.0.0.sif
      if present, otherwise docker://papaemmelab/nf-ffperase:v1.0.0.

  --engine auto|apptainer|singularity|docker
      Default: auto.

  --repository REPO
      Default: papaemmelab/nf-ffperase

  --revision REVISION
      Default: b0dd56cbd0a939896a966b9ce30c4d719b158170

  --models-dir PATH
      Directory containing or receiving model.snvs.joblib/model.indels.joblib.
      Default: <outroot>/models

  --download-models true|false
      Default: true.

Variant type options
--------------------
  --snvs-only
  --indels-only
  --skip-snvs
  --skip-indels

Picard/GATK options
-------------------
  --java-mem MEM          Default: 12g
  --picard-timeout-min N  Default: 720
  --backend gatk|auto|jar|picard
  --conda-env NAME_OR_PATH
  --project-dir PATH      Optional project dir to help locate an existing gatk.
  --gatk-bin PATH
  --picard-jar PATH
  --no-conda
  --no-create-conda-env

Quality options
---------------
  --min-mapq N     Default: 20
  --min-baseq N    Default: 20
  --min-depth N    Default: 1
  --split-reads N  Default: 7500000

Examples
--------
  1. One-sample test:

    ROOT=/path/to/variant-calling-output
    REF=/path/to/reference/genome.fa

    bash bin/run_ffperase_single_picard_pileup_nf.sh \
      --root "$ROOT" \
      --ref "$REF" \
      --samples SAMPLE_01 \
      --threads 32 \
      --picard-jobs 1 \
      --sample-jobs 1 \
      --split-pileup 5000

  2. Full lymphoma WGS:

    LYMPHOMA_SAMPLES="SAMPLE_01,SAMPLE_02,SAMPLE_03"

    bash bin/run_ffperase_single_picard_pileup_nf.sh \
      --root "$ROOT" \
      --ref "$REF" \
      --samples "$LYMPHOMA_SAMPLES" \
      --threads 32 \
      --picard-jobs 1 \
      --sample-jobs 1 \
      --split-pileup 5000

  3. Classification only after features already exist:

    bash bin/run_ffperase_single_picard_pileup_nf.sh \
      --root "$ROOT" \
      --ref "$REF" \
      --samples SAMPLE_01 \
      --classify-only \
      --classify-threads 2
EOF_HELP
}

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*" >&2
  [[ -n "${LOG:-}" ]] && printf '[%(%F %T)T] %s\n' -1 "$*" >> "$LOG" || true
}

die() { log "ERROR: $*"; exit 1; }

as_bool() {
  local v="${1:-}"
  v="$(printf '%s' "$v" | tr '[:upper:]' '[:lower:]')"
  case "$v" in
    true|1|yes|y) echo true ;;
    false|0|no|n) echo false ;;
    *) die "Invalid boolean: $1" ;;
  esac
}

abs_path() {
  realpath -m "$1" 2>/dev/null || readlink -m "$1" 2>/dev/null || python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$1"
}

extract_embedded_tool() {
  local start_marker="$1"
  local end_marker="$2"
  local output_file="$3"
  awk -v start="$start_marker" -v end="$end_marker" '
    $0 == start {inside=1; next}
    $0 == end {inside=0; found=1; next}
    inside {print}
    END { if (!found) exit 42 }
  ' "$SCRIPT_PATH" > "$output_file" || die "Could not extract embedded tool to $output_file"
  chmod +x "$output_file"
}

extract_tools() {
  mkdir -p "$TOOLS_DIR"
  extract_embedded_tool "###__FFPERASE_EMBEDDED_PICARD_BEGIN__###" \
                        "###__FFPERASE_EMBEDDED_PICARD_END__###" \
                        "$TOOLS_DIR/precompute_ffperase_picard_metrics_by_sample_v3_for_ffperase_nf.sh"
  extract_embedded_tool "###__FFPERASE_EMBEDDED_NF_WRAPPER_BEGIN__###" \
                        "###__FFPERASE_EMBEDDED_NF_WRAPPER_END__###" \
                        "$TOOLS_DIR/run_ffperase_post_pipeline_v8_nf_threads.sh"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --outroot) OUTROOT="$2"; shift 2 ;;
    --picard-metrics-root|--picard-root) PICARD_ROOT="$2"; shift 2 ;;
    --features-root|--features-cache-root) FEATURES_ROOT="$2"; shift 2 ;;
    --models-dir) MODELS_DIR="$2"; shift 2 ;;
    --tools-dir) TOOLS_DIR="$2"; shift 2 ;;

    --threads) THREADS="$2"; shift 2 ;;
    --classify-threads) CLASSIFY_THREADS="$2"; shift 2 ;;
    --sample-jobs|--jobs) SAMPLE_JOBS="$2"; shift 2 ;;
    --picard-jobs) PICARD_JOBS="$2"; shift 2 ;;
    --split-pileup) SPLIT_PILEUP="$2"; shift 2 ;;
    --split-reads) SPLIT_READS="$2"; shift 2 ;;
    --java-mem) JAVA_MEM="$2"; shift 2 ;;
    --timeout-min) TIMEOUT_MIN="$2"; shift 2 ;;
    --picard-timeout-min) PICARD_TIMEOUT_MIN="$2"; shift 2 ;;
    --classify-timeout-min) CLASSIFY_TIMEOUT_MIN="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-baseq) MIN_BASEQ="$2"; shift 2 ;;
    --min-depth) MIN_DEPTH="$2"; shift 2 ;;

    --skip-picard) RUN_PICARD=false; shift ;;
    --skip-preprocess|--skip-pileup) RUN_PREPROCESS=false; shift ;;
    --skip-classify) RUN_CLASSIFY=false; shift ;;
    --picard-only) RUN_PICARD=true; RUN_PREPROCESS=false; RUN_CLASSIFY=false; shift ;;
    --preprocess-only) RUN_PICARD=true; RUN_PREPROCESS=true; RUN_CLASSIFY=false; shift ;;
    --classify-only) RUN_PICARD=false; RUN_PREPROCESS=false; RUN_CLASSIFY=true; shift ;;
    --force) FORCE_PICARD=true; FORCE_PREPROCESS=true; FORCE_CLASSIFY=true; shift ;;
    --force-picard) FORCE_PICARD=true; shift ;;
    --force-preprocess|--force-pileup) FORCE_PREPROCESS=true; shift ;;
    --force-classify) FORCE_CLASSIFY=true; shift ;;
    --use-resume) USE_RESUME="$(as_bool "$2")"; shift 2 ;;
    --fail-on-error) FAIL_ON_ERROR="$(as_bool "$2")"; shift 2 ;;
    --no-recover-missing-features) RECOVER_MISSING_FEATURES=false; shift ;;
    --uppercase-reference) UPPERCASE_REFERENCE="$(as_bool "$2")"; shift 2 ;;
    --sanitize-snvs) SANITIZE_SNVS="$(as_bool "$2")"; shift 2 ;;

    --container) CONTAINER="$2"; shift 2 ;;
    --engine) ENGINE="$2"; shift 2 ;;
    --repository|--repo) REPOSITORY="$2"; shift 2 ;;
    --revision|-r) REVISION="$2"; shift 2 ;;
    --download-models) DOWNLOAD_MODELS="$(as_bool "$2")"; shift 2 ;;
    --prefer-local-model) PREFER_LOCAL_MODEL="$(as_bool "$2")"; shift 2 ;;

    --snvs-only) RUN_SNVS=true; RUN_INDELS=false; shift ;;
    --indels-only) RUN_SNVS=false; RUN_INDELS=true; shift ;;
    --skip-snvs) RUN_SNVS=false; shift ;;
    --skip-indels) RUN_INDELS=false; shift ;;

    --backend) PICARD_BACKEND="$2"; shift 2 ;;
    --conda-env) PICARD_CONDA_ENV="$2"; shift 2 ;;
    --project-dir) PICARD_PROJECT_DIR="$2"; shift 2 ;;
    --gatk-bin) PICARD_GATK_BIN="$2"; shift 2 ;;
    --picard-jar) PICARD_JAR="$2"; shift 2 ;;
    --no-conda) PICARD_NO_CONDA=true; shift ;;
    --no-create-conda-env) PICARD_NO_CREATE_CONDA=true; shift ;;

    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$ROOT" ]] || { usage >&2; die "--root is required"; }
[[ -n "$REF" ]] || { usage >&2; die "--ref is required"; }
ROOT="$(abs_path "$ROOT")"
REF="$(abs_path "$REF")"
[[ -d "$ROOT" ]] || die "Root not found: $ROOT"
[[ -s "$REF" ]] || die "Reference not found/empty: $REF"

OUTROOT="${OUTROOT:-$ROOT/ffperase_classification}"
OUTROOT="$(abs_path "$OUTROOT")"
PICARD_ROOT="${PICARD_ROOT:-$ROOT/ffperase_picard_metrics}"
FEATURES_ROOT="${FEATURES_ROOT:-$OUTROOT/features_cache}"
MODELS_DIR="${MODELS_DIR:-$OUTROOT/models}"
TOOLS_DIR="${TOOLS_DIR:-$OUTROOT/.embedded_tools}"
PICARD_ROOT="$(abs_path "$PICARD_ROOT")"
FEATURES_ROOT="$(abs_path "$FEATURES_ROOT")"
MODELS_DIR="$(abs_path "$MODELS_DIR")"
TOOLS_DIR="$(abs_path "$TOOLS_DIR")"

mkdir -p "$OUTROOT" "$PICARD_ROOT" "$FEATURES_ROOT" "$MODELS_DIR" "$TOOLS_DIR"
LOG="$OUTROOT/ffperase_single_picard_pileup_nf.log"
: > "$LOG"

# Portable cache locations for the local helpers. NXF_VER is intentionally not
# set here: callers may select a Nextflow version, otherwise the installed
# executable chooses its own version. The legacy parser is forced only on the
# actual nested nf-ffperase commands below.
unset NXF_PLUGINS_DIR || true
unset NXF_WORK || true
export NXF_HOME="${NXF_HOME:-${HOME:?HOME must be set}/.nextflow}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${HOME}/.apptainer/cache}"
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-${HOME}/.singularity/cache}"
mkdir -p "$NXF_HOME/plugins" "$APPTAINER_CACHEDIR" "$SINGULARITY_CACHEDIR"

extract_tools
PICARD_SCRIPT="$TOOLS_DIR/precompute_ffperase_picard_metrics_by_sample_v3_for_ffperase_nf.sh"
FFPERASE_WRAPPER="$TOOLS_DIR/run_ffperase_post_pipeline_v8_nf_threads.sh"

log "Script             : $SCRIPT_NAME 2026-06-01-v3-sanitize-ambiguous-snvs-for-ffperase"
log "Root               : $ROOT"
log "Reference          : $REF"
log "Outroot            : $OUTROOT"
log "Picard metrics root: $PICARD_ROOT"
log "Features root      : $FEATURES_ROOT"
log "Models dir         : $MODELS_DIR"
log "Embedded tools dir : $TOOLS_DIR"
log "Samples            : ${SAMPLES:-auto}"
log "Threads            : $THREADS for nf-ffperase preprocess/PILEUP maxForks"
log "Classify threads   : $CLASSIFY_THREADS"
log "Picard jobs        : $PICARD_JOBS sample-level jobs"
log "Sample jobs        : $SAMPLE_JOBS sample/type jobs in wrapper"
log "Split pileup       : $SPLIT_PILEUP"
log "Stages             : picard=$RUN_PICARD preprocess/pileup=$RUN_PREPROCESS classify=$RUN_CLASSIFY"
log "Sanitize SNVs      : $SANITIZE_SNVS"
log "Fail on error      : $FAIL_ON_ERROR"
log "Recover missing    : $RECOVER_MISSING_FEATURES"
log "Nextflow           : $(nextflow -version 2>/dev/null | tr '\n' ' ' | sed 's/[[:space:]][[:space:]]*/ /g')"

sample_args=()
[[ -n "$SAMPLES" ]] && sample_args=(--samples "$SAMPLES")

type_args=()
if [[ "$RUN_SNVS" == "true" && "$RUN_INDELS" == "false" ]]; then type_args=(--snvs-only); fi
if [[ "$RUN_SNVS" == "false" && "$RUN_INDELS" == "true" ]]; then type_args=(--indels-only); fi
if [[ "$RUN_SNVS" == "false" && "$RUN_INDELS" == "false" ]]; then die "Both SNVs and indels are disabled"; fi

force_picard_args=()
[[ "$FORCE_PICARD" == "true" ]] && force_picard_args=(--force)
force_preprocess_args=()
[[ "$FORCE_PREPROCESS" == "true" ]] && force_preprocess_args=(--force)
force_classify_args=()
[[ "$FORCE_CLASSIFY" == "true" ]] && force_classify_args=(--force)

picard_tool_args=(--backend "$PICARD_BACKEND" --conda-env "$PICARD_CONDA_ENV")
[[ -n "$PICARD_PROJECT_DIR" ]] && picard_tool_args+=(--project-dir "$PICARD_PROJECT_DIR")
[[ -n "$PICARD_GATK_BIN" ]] && picard_tool_args+=(--gatk-bin "$PICARD_GATK_BIN")
[[ -n "$PICARD_JAR" ]] && picard_tool_args+=(--picard-jar "$PICARD_JAR")
[[ "$PICARD_NO_CONDA" == "true" ]] && picard_tool_args+=(--no-conda)
[[ "$PICARD_NO_CREATE_CONDA" == "true" ]] && picard_tool_args+=(--no-create-conda-env)

ffperase_common_args=(
  --container "$CONTAINER"
  --engine "$ENGINE"
  --repository "$REPOSITORY"
  --revision "$REVISION"
  --download-models "$DOWNLOAD_MODELS"
  --prefer-local-model "$PREFER_LOCAL_MODEL"
)

# Build a sample list for validation/recovery without relying on the nested wrapper.
write_validation_samples() {
  local out="$1"
  : > "$out"
  if [[ -n "$SAMPLES" ]]; then
    printf '%s' "$SAMPLES" | tr ',;' '\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | awk 'NF' | sort -u > "$out"
  else
    find "$ROOT/final_vcf" -maxdepth 1 -type f \( -name '*.final_variants.vcf.gz' -o -name '*.final_variants.vcf' \) \
      | sed 's#.*/##; s/\.final_variants\.vcf\.gz$//; s/\.final_variants\.vcf$//' \
      | sort -u > "$out"
  fi
  [[ -s "$out" ]] || die "Could not infer any samples for validation"
}

resolve_final_vcf_for_sample() {
  local sample="$1"
  local f=""
  for f in \
    "$ROOT/final_vcf/${sample}.final_variants.vcf.gz" \
    "$ROOT/final_vcf/${sample}.final_variants.vcf"; do
    if [[ -s "$f" ]]; then
      printf '%s\n' "$f"
      return 0
    fi
  done
  return 1
}

vcf_has_type_variants() {
  local vcf="$1"
  local mtype="$2"
  local kind="snps"
  [[ "$mtype" == "indels" ]] && kind="indels"
  local first=""
  set +e
  set +o pipefail
  if [[ "$mtype" == "snvs" && "$SANITIZE_SNVS" == "true" ]]; then
    # Only expect FFPErase SNV features if at least one biallelic A/C/G/T SNV exists.
    first="$(bcftools view -H -v snps "$vcf" 2>/dev/null | awk -F'\t' 'length($4)==1 && length($5)==1 && $4 ~ /^[ACGTacgt]$/ && $5 ~ /^[ACGTacgt]$/ {print; exit}' | head -n 1)"
  else
    first="$(bcftools view -H -v "$kind" "$vcf" 2>/dev/null | head -n 1)"
  fi
  local st=$?
  set -o pipefail
  set -e
  [[ "$st" -eq 0 && -n "$first" ]]
}
feature_present_or_cache() {
  local sample="$1"
  local mtype="$2"
  local cache="$FEATURES_ROOT/$sample/$mtype/features.tsv"
  local fallback="$OUTROOT/$sample/$mtype/preprocess/features.tsv"
  if [[ -s "$cache" ]]; then
    return 0
  fi
  if [[ -s "$fallback" ]]; then
    mkdir -p "$(dirname "$cache")"
    cp -f "$fallback" "$cache"
    log "Recovered features cache from nested preprocess output: sample=$sample type=$mtype cache=$cache"
    return 0
  fi
  return 1
}

write_missing_features_report() {
  local report="$1"
  local samples_file="$OUTROOT/.validation_samples.txt"
  write_validation_samples "$samples_file"
  printf 'sample\ttype\tfinal_vcf\texpected\tfeatures\n' > "$report"
  local sample vcf mtype feature_path
  while IFS= read -r sample; do
    [[ -n "$sample" ]] || continue
    if ! vcf="$(resolve_final_vcf_for_sample "$sample")"; then
      printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "." "MISSING_FINAL_VCF" "yes" "." >> "$report"
      continue
    fi
    for mtype in snvs indels; do
      if [[ "$mtype" == "snvs" && "$RUN_SNVS" != "true" ]]; then continue; fi
      if [[ "$mtype" == "indels" && "$RUN_INDELS" != "true" ]]; then continue; fi
      if vcf_has_type_variants "$vcf" "$mtype"; then
        feature_path="$FEATURES_ROOT/$sample/$mtype/features.tsv"
        if ! feature_present_or_cache "$sample" "$mtype"; then
          printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "$vcf" "yes" "$feature_path" >> "$report"
        fi
      fi
    done
  done < "$samples_file"
  # Return non-zero if there are rows after the header.
  [[ "$(awk 'NR>1{n++} END{print n+0}' "$report")" -eq 0 ]]
}

run_ffperase_preprocess_for_type() {
  local sample="$1"
  local mtype="$2"
  local local_type_args=()
  if [[ "$mtype" == "snvs" ]]; then
    local_type_args=(--snvs-only)
  elif [[ "$mtype" == "indels" ]]; then
    local_type_args=(--indels-only)
  else
    die "Internal error: invalid mtype for recovery: $mtype"
  fi

  log "RECOVERY: rerunning nf-ffperase preprocess-only for sample=$sample type=$mtype"
  bash "$FFPERASE_WRAPPER" \
    --root "$ROOT" \
    --outroot "$OUTROOT" \
    --ref "$REF" \
    --samples "$sample" \
    --picard-metrics-root "$PICARD_ROOT" \
    --features-cache-root "$FEATURES_ROOT" \
    --preprocess-only true \
    --jobs 1 \
    --threads "$THREADS" \
    --split-pileup "$SPLIT_PILEUP" \
    --split-reads "$SPLIT_READS" \
    --timeout-min "$TIMEOUT_MIN" \
    --min-mapq "$MIN_MAPQ" \
    --min-baseq "$MIN_BASEQ" \
    --min-depth "$MIN_DEPTH" \
    --uppercase-reference "$UPPERCASE_REFERENCE" \
    --sanitize-snvs "$SANITIZE_SNVS" \
    --ffperase-model-name ARTIFACT \
    --models-dir "$MODELS_DIR" \
    --use-resume false \
    --fail-on-error true \
    "${local_type_args[@]}" \
    "${ffperase_common_args[@]}" \
    --force
}

validate_and_recover_features() {
  local tag="$1"
  local report="$OUTROOT/missing_features_${tag}.tsv"
  if write_missing_features_report "$report"; then
    log "Feature validation OK after $tag: all expected SNV/indel features.tsv files exist"
    return 0
  fi

  log "Feature validation found missing features after $tag: $report"
  cat "$report" >&2 || true

  if [[ "$RECOVER_MISSING_FEATURES" == "true" && "$RUN_PREPROCESS" == "true" ]]; then
    log "Attempting one recovery pass for missing features.tsv files"
    local sample mtype vcf expected feature
    tail -n +2 "$report" | while IFS=$'\t' read -r sample mtype vcf expected feature; do
      [[ -n "$sample" && "$sample" != "." && "$mtype" != "." ]] || continue
      run_ffperase_preprocess_for_type "$sample" "$mtype"
    done

    local report2="$OUTROOT/missing_features_${tag}_after_recovery.tsv"
    if write_missing_features_report "$report2"; then
      log "Feature validation OK after recovery"
      return 0
    fi
    log "Missing features remain after recovery: $report2"
    cat "$report2" >&2 || true
    die "Cannot continue to classification because at least one expected features.tsv is missing. Inspect the corresponding console log under $OUTROOT/<sample>/<type>/"
  fi

  if [[ "$FAIL_ON_ERROR" == "true" ]]; then
    die "Cannot continue because expected FFPErase features are missing. Use --no-recover-missing-features only for diagnostics, or rerun the failed type with --force-preprocess --snvs-only/--indels-only."
  fi

  log "WARNING: continuing despite missing features because --fail-on-error false"
  return 0
}

classified_present() {
  local sample="$1"
  local mtype="$2"
  [[ -s "$OUTROOT/$sample/$mtype/classify/classified_df_${mtype}.tsv" || -s "$OUTROOT/classified/${sample}.classified_df_${mtype}.tsv" ]]
}

validate_classified_outputs() {
  local report="$OUTROOT/missing_classified_after_classify.tsv"
  local samples_file="$OUTROOT/.validation_samples.txt"
  write_validation_samples "$samples_file"
  printf 'sample\ttype\texpected_features\tclassified\n' > "$report"
  local sample mtype
  while IFS= read -r sample; do
    [[ -n "$sample" ]] || continue
    for mtype in snvs indels; do
      if [[ "$mtype" == "snvs" && "$RUN_SNVS" != "true" ]]; then continue; fi
      if [[ "$mtype" == "indels" && "$RUN_INDELS" != "true" ]]; then continue; fi
      if feature_present_or_cache "$sample" "$mtype"; then
        if ! classified_present "$sample" "$mtype"; then
          printf '%s\t%s\t%s\t%s\n' "$sample" "$mtype" "$FEATURES_ROOT/$sample/$mtype/features.tsv" "$OUTROOT/$sample/$mtype/classify/classified_df_${mtype}.tsv" >> "$report"
        fi
      fi
    done
  done < "$samples_file"
  if [[ "$(awk 'NR>1{n++} END{print n+0}' "$report")" -ne 0 ]]; then
    log "Classification validation found missing classified outputs: $report"
    cat "$report" >&2 || true
    if [[ "$FAIL_ON_ERROR" == "true" ]]; then
      die "Classification ended with missing outputs. Inspect console logs under $OUTROOT/<sample>/<type>/"
    fi
  else
    log "Classification validation OK: all cached feature sets have classified outputs"
  fi
}

if [[ "$RUN_PICARD" == "true" ]]; then
  log "PHASE 1A: Picard/GATK CollectSequencingArtifactMetrics precompute"
  bash "$PICARD_SCRIPT" \
    --root "$ROOT" \
    --ref "$REF" \
    --outroot "$PICARD_ROOT" \
    "${sample_args[@]}" \
    --jobs "$PICARD_JOBS" \
    --java-mem "$JAVA_MEM" \
    --min-mapq "$MIN_MAPQ" \
    --min-baseq "$MIN_BASEQ" \
    --timeout-min "$PICARD_TIMEOUT_MIN" \
    --uppercase-reference "$UPPERCASE_REFERENCE" \
    --keep-going \
    "${picard_tool_args[@]}" \
    "${force_picard_args[@]}"
fi

if [[ "$RUN_PREPROCESS" == "true" ]]; then
  log "PHASE 1B: nf-ffperase preprocess-only: PILEUP + ANNOTATE_VARIANTS + cached features.tsv"
  bash "$FFPERASE_WRAPPER" \
    --root "$ROOT" \
    --outroot "$OUTROOT" \
    --ref "$REF" \
    "${sample_args[@]}" \
    --picard-metrics-root "$PICARD_ROOT" \
    --features-cache-root "$FEATURES_ROOT" \
    --preprocess-only true \
    --jobs "$SAMPLE_JOBS" \
    --threads "$THREADS" \
    --split-pileup "$SPLIT_PILEUP" \
    --split-reads "$SPLIT_READS" \
    --timeout-min "$TIMEOUT_MIN" \
    --min-mapq "$MIN_MAPQ" \
    --min-baseq "$MIN_BASEQ" \
    --min-depth "$MIN_DEPTH" \
    --uppercase-reference "$UPPERCASE_REFERENCE" \
    --sanitize-snvs "$SANITIZE_SNVS" \
    --ffperase-model-name ARTIFACT \
    --models-dir "$MODELS_DIR" \
    --use-resume "$USE_RESUME" \
    --fail-on-error "$FAIL_ON_ERROR" \
    "${type_args[@]}" \
    "${ffperase_common_args[@]}" \
    "${force_preprocess_args[@]}"
fi

if [[ "$RUN_PREPROCESS" == "true" && "$RUN_CLASSIFY" == "true" ]]; then
  validate_and_recover_features "after_preprocess"
fi

if [[ "$RUN_CLASSIFY" == "true" ]]; then
  if [[ "$RUN_PREPROCESS" != "true" ]]; then
    validate_and_recover_features "before_classify"
  fi
  log "PHASE 2: nf-ffperase classify-only from cached features.tsv"
  bash "$FFPERASE_WRAPPER" \
    --root "$ROOT" \
    --outroot "$OUTROOT" \
    --ref "$REF" \
    "${sample_args[@]}" \
    --features-root "$FEATURES_ROOT" \
    --classify-only true \
    --jobs "$SAMPLE_JOBS" \
    --threads "$CLASSIFY_THREADS" \
    --timeout-min "$CLASSIFY_TIMEOUT_MIN" \
    --uppercase-reference "$UPPERCASE_REFERENCE" \
    --sanitize-snvs "$SANITIZE_SNVS" \
    --ffperase-model-name ARTIFACT \
    --models-dir "$MODELS_DIR" \
    --use-resume "$USE_RESUME" \
    --fail-on-error "$FAIL_ON_ERROR" \
    "${type_args[@]}" \
    "${ffperase_common_args[@]}" \
    "${force_classify_args[@]}"
fi

if [[ "$RUN_CLASSIFY" == "true" ]]; then
  validate_classified_outputs
fi

log "Done. Main classified summaries, if generated:"
find "$OUTROOT/classified" -maxdepth 1 -type f -name 'all_samples.ffperase_*.tsv' -print 2>/dev/null | sort | tee -a "$LOG" || true

exit 0

###__FFPERASE_EMBEDDED_PICARD_BEGIN__###
#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# precompute_ffperase_picard_metrics_by_sample_v3.sh
#
# Purpose:
#   Precompute Picard/GATK CollectSequencingArtifactMetrics once per BAM/sample
#   and write FFPErase-ready metric tables:
#
#     <outroot>/<sample>/pre_adapter_metrics.tsv
#     <outroot>/<sample>/bait_bias_metrics.tsv
#
# Why:
#   For genome-wide VCFs, nf-ffperase can explode into hundreds of thousands of
#   tiny Picard jobs. This script runs one Picard/GATK job per sample, controlled
#   by --jobs, and prepares metrics for later use with:
#
#     run_ffperase_post_pipeline.sh --picard-metrics-root <outroot>
#
# Key fixes in v3:
#   - Creates a dedicated conda env if it does not exist.
#   - Defaults to GATK CollectSequencingArtifactMetrics, not fragile `picard`.
#   - Builds/uses an uppercase reference by default to avoid FFPErase context
#     errors from lowercase reference sequence.
#   - Prints useful stderr tails when Picard/GATK aborts.
#   - Ctrl+C / TERM cleanup kills child process groups.
#
# Version: 2026-05-31-v3-conda-autocreate
###############################################################################

ROOT=""
BAM_DIR=""
BAM_SUFFIX=".preprocessed.bam"
OUTROOT=""
REF=""
DBSNP=""
SAMPLES_CSV=""
SAMPLES_FILE_IN=""
SAMPLES_FROM_DIR=""
JOBS=1
JAVA_MEM="12g"
MIN_MAPQ=20
MIN_BASEQ=20
SAMTOOLS_THREADS=4
PICARD_JAR=""
PICARD_BACKEND="gatk"
GATK_BIN=""
PICARD_BIN=""
PROJECT_DIR=""
CONDA_ENV="ffperase_picard_metrics"
USE_CONDA=true
CREATE_CONDA_ENV=true
CONDA_FRONTEND_FORCED=""
FORCE=false
SKIP_MISSING=false
KEEP_GOING=false
DRY_RUN=false
TIMEOUT_MIN=0
TMP_ROOT=""
UPPERCASE_REFERENCE=true

SCRIPT_NAME="$(basename "$0")"

usage() {
cat <<'EOF_HELP'
precompute_ffperase_picard_metrics_by_sample_v3.sh

Purpose
-------
Precompute Picard/GATK CollectSequencingArtifactMetrics once per sample and
format outputs for nf-ffperase / run_ffperase_post_pipeline.sh.

Expected default input layout
-----------------------------
  <root>/preprocessed_bam/<sample>.preprocessed.bam
  <root>/preprocessed_bam/<sample>.preprocessed.bam.bai

Main output per sample
----------------------
  <outroot>/<sample>/pre_adapter_metrics.tsv
  <outroot>/<sample>/bait_bias_metrics.tsv
  <outroot>/<sample>/logs/<sample>.picard.stderr.log
  <outroot>/picard_metrics_manifest.tsv

Required
--------
  --root PATH
      Variant-calling output directory. By default the script looks for BAMs in:
      <root>/preprocessed_bam

  --ref FASTA
      Reference FASTA used for alignment.

Recommended setup
-----------------
  Default behavior creates/uses this conda env if missing:

      ffperase_picard_metrics

  The env includes gatk4, samtools, htslib, python, and openjdk. Then the script
  runs GATK CollectSequencingArtifactMetrics.

Conda options
-------------
  --conda-env NAME_OR_PATH
      Conda environment name or full prefix path.
      Default: ffperase_picard_metrics

      Examples:
        --conda-env ffperase_picard_metrics
        --conda-env /path/to/conda/envs/ffperase_picard_metrics

  --create-conda-env true|false
      Create the conda env if it is absent. Default: true

  --no-create-conda-env
      Do not create the env if absent; fail if required tools are missing.

  --no-conda
      Do not activate or create a conda env. Use tools already in PATH or those
      supplied by --gatk-bin / --picard-jar.

  --conda-frontend conda|mamba
      Force conda or mamba for env creation. Default: auto, prefers mamba if
      available.

Sample selection options
------------------------
  --samples LIST
      Comma-separated or semicolon-separated sample IDs.
      Example: --samples SAMPLE_01,SAMPLE_02,SAMPLE_03

  --samples-file FILE
      Text file with one sample ID per line. Blank lines and # comments ignored.
      TSV/CSV is accepted; the first column is interpreted as sample ID or BAM.

  --samples-from-dir DIR
      Infer sample IDs from immediate children of a directory. Suffixes such as
      .final_variants.vcf.gz, .preprocessed.bam, .bam, .vcf.gz, .vcf are stripped.

  If none of --samples, --samples-file, or --samples-from-dir is provided,
  samples are inferred from <bam-dir>/*<bam-suffix>.

Input directory options
-----------------------
  --bam-dir PATH
      Directory containing BAMs. Default: <root>/preprocessed_bam

  --bam-suffix STR
      BAM suffix used to resolve sample IDs. Default: .preprocessed.bam

Backend/tool options
--------------------
  --backend gatk|auto|jar|picard
      gatk  : force GATK CollectSequencingArtifactMetrics. Default and safest.
      auto  : prefer GATK, then picard.jar, then picard wrapper.
      jar   : force java -jar picard.jar.
      picard: force picard wrapper.

  --gatk-bin FILE
      Explicit GATK executable path.

  --picard-jar FILE
      Explicit picard.jar path. Required if --backend jar and auto-detection fails.

  --project-dir PATH
      Optional Nextflow project directory. If supplied, the script also searches:
      <project-dir>/.conda/env-*/bin/gatk

Output/options
--------------
  --outroot PATH
      Output directory for FFPErase-ready Picard metrics.
      Default: <root>/ffperase_picard_metrics

  --uppercase-reference true|false
      Build and use an uppercase copy of the reference for Picard metrics.
      This avoids FFPErase trinucleotide context errors from lowercase genome
      sequence. Default: true

  --dbsnp VCF[.gz]
      Optional dbSNP VCF for Picard/GATK DB_SNP. If supplied, Picard/GATK ignores
      known polymorphic positions when estimating sequencing artifacts.

  --jobs N
      Number of samples processed in parallel. Recommended for WGS: 1-2.
      Default: 1

  --java-mem MEM
      Java memory per Picard/GATK job. Default: 12g

  --min-mapq N
      Picard MINIMUM_MAPPING_QUALITY. Default: 20

  --min-baseq N
      Picard MINIMUM_QUALITY_SCORE. Default: 20

  --samtools-threads N
      Threads used only for samtools indexing if BAM index is missing. Default: 4

  --timeout-min N
      Optional timeout per sample Picard/GATK job in minutes. 0 disables timeout.
      Default: 0

  --tmp-root PATH
      Optional temporary root. Default: <outroot>/<sample>/tmp

  --force
      Recompute metrics even if outputs already exist.

  --skip-missing
      Skip samples whose BAM is missing instead of failing.

  --keep-going
      Continue processing other samples if one sample fails.

  --dry-run
      Build and print the sample/BAM manifest only; do not run Picard/GATK.

  -h, --help
      Show this help.

Examples
--------
  1. Test one sample and create conda env automatically if absent:

    bash precompute_ffperase_picard_metrics_by_sample_v3.sh \
      --root /path/to/variant-calling-output \
      --ref /path/to/reference/genome.fa \
      --samples SAMPLE_01 \
      --jobs 1 \
      --java-mem 12g \
      --timeout-min 720 \
      --force

  2. Recommended lymphoma WGS run:

    LYMPHOMA_SAMPLES="SAMPLE_01,SAMPLE_02,SAMPLE_03"

    bash precompute_ffperase_picard_metrics_by_sample_v3.sh \
      --root /path/to/variant-calling-output \
      --ref /path/to/reference/genome.fa \
      --samples "$LYMPHOMA_SAMPLES" \
      --jobs 1 \
      --java-mem 12g \
      --timeout-min 720 \
      --keep-going

Then run FFPErase post-pipeline with precomputed metrics
--------------------------------------------------------
  ROOT=/path/to/variant-calling-output
  REF=/path/to/reference/genome.fa

  bash run_ffperase_post_pipeline.sh \
    --root "$ROOT" \
    --outroot "$ROOT/ffperase_classification" \
    --ref "$REF" \
    --picard-metrics-root "$ROOT/ffperase_picard_metrics" \
    --jobs 1 \
    --max-forks 2 \
    --timeout-min 720 \
    --uppercase-reference true \
    --ffperase-model-name ARTIFACT \
    --use-resume false \
    --fail-on-error false
EOF_HELP
}

log() { printf '[%(%F %T)T] %s\n' -1 "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }
truthy() { [[ "${1:-}" =~ ^([Tt][Rr][Uu][Ee]|1|[Yy][Ee][Ss]|[Yy])$ ]]; }
falsy() { [[ "${1:-}" =~ ^([Ff][Aa][Ll][Ss][Ee]|0|[Nn][Oo]|[Nn])$ ]]; }

# Predictable locale for Java/Picard/GATK.
export LC_ALL="${LC_ALL:-C.UTF-8}"
export LANG="${LANG:-C.UTF-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$2"; shift 2 ;;
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --bam-suffix) BAM_SUFFIX="$2"; shift 2 ;;
    --outroot) OUTROOT="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --dbsnp|--db-snp) DBSNP="$2"; shift 2 ;;
    --samples) SAMPLES_CSV="$2"; shift 2 ;;
    --samples-file|--sample-file) SAMPLES_FILE_IN="$2"; shift 2 ;;
    --samples-from-dir|--sample-dir|--samples-dir) SAMPLES_FROM_DIR="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --java-mem) JAVA_MEM="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-baseq) MIN_BASEQ="$2"; shift 2 ;;
    --samtools-threads) SAMTOOLS_THREADS="$2"; shift 2 ;;
    --backend) PICARD_BACKEND="$2"; shift 2 ;;
    --gatk-bin) GATK_BIN="$2"; shift 2 ;;
    --picard-jar) PICARD_JAR="$2"; shift 2 ;;
    --project-dir) PROJECT_DIR="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --conda-prefix) CONDA_ENV="$2"; shift 2 ;;
    --create-conda-env)
      if truthy "$2"; then CREATE_CONDA_ENV=true; elif falsy "$2"; then CREATE_CONDA_ENV=false; else die "--create-conda-env must be true or false"; fi
      shift 2 ;;
    --no-create-conda-env) CREATE_CONDA_ENV=false; shift ;;
    --no-conda) USE_CONDA=false; CREATE_CONDA_ENV=false; shift ;;
    --conda-frontend) CONDA_FRONTEND_FORCED="$2"; shift 2 ;;
    --timeout-min) TIMEOUT_MIN="$2"; shift 2 ;;
    --tmp-root) TMP_ROOT="$2"; shift 2 ;;
    --uppercase-reference)
      if truthy "$2"; then UPPERCASE_REFERENCE=true; elif falsy "$2"; then UPPERCASE_REFERENCE=false; else die "--uppercase-reference must be true or false"; fi
      shift 2 ;;
    --force) FORCE=true; shift ;;
    --skip-missing) SKIP_MISSING=true; shift ;;
    --keep-going) KEEP_GOING=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$ROOT" ]] || die "--root is required"
[[ -n "$REF" ]] || die "--ref is required"
[[ "$JOBS" =~ ^[0-9]+$ ]] || die "--jobs must be an integer"
[[ "$JOBS" -ge 1 ]] || die "--jobs must be >= 1"
[[ "$SAMTOOLS_THREADS" =~ ^[0-9]+$ ]] || die "--samtools-threads must be an integer"
[[ "$SAMTOOLS_THREADS" -ge 1 ]] || die "--samtools-threads must be >= 1"
[[ "$TIMEOUT_MIN" =~ ^[0-9]+$ ]] || die "--timeout-min must be an integer"
case "$PICARD_BACKEND" in auto|gatk|jar|picard) ;; *) die "--backend must be auto, gatk, jar, or picard" ;; esac
case "$CONDA_FRONTEND_FORCED" in ""|conda|mamba) ;; *) die "--conda-frontend must be conda or mamba" ;; esac

ROOT="$(readlink -m "$ROOT")"
REF="$(readlink -m "$REF")"
BAM_DIR="${BAM_DIR:-$ROOT/preprocessed_bam}"
BAM_DIR="$(readlink -m "$BAM_DIR")"
OUTROOT="${OUTROOT:-$ROOT/ffperase_picard_metrics}"
OUTROOT="$(readlink -m "$OUTROOT")"
[[ -n "$PROJECT_DIR" ]] && PROJECT_DIR="$(readlink -m "$PROJECT_DIR")"
[[ -n "$TMP_ROOT" ]] && TMP_ROOT="$(readlink -m "$TMP_ROOT")"

[[ -d "$ROOT" ]] || die "Root directory not found: $ROOT"
[[ -d "$BAM_DIR" ]] || die "BAM directory not found: $BAM_DIR"
[[ -s "$REF" ]] || die "Reference FASTA not found/empty: $REF"
mkdir -p "$OUTROOT" "$OUTROOT/logs"

###############################################################################
# Conda environment creation / activation.
###############################################################################

find_conda_sh() {
  local base=""
  if command -v conda >/dev/null 2>&1; then
    base="$(conda info --base 2>/dev/null || true)"
    [[ -n "$base" && -f "$base/etc/profile.d/conda.sh" ]] && { echo "$base/etc/profile.d/conda.sh"; return 0; }
  fi
  for f in \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh"
  do
    [[ -f "$f" ]] && { echo "$f"; return 0; }
  done
  return 1
}

conda_env_exists() {
  local env="$1"
  if [[ "$env" == */* ]]; then
    [[ -d "$env/conda-meta" || -x "$env/bin/gatk" || -x "$env/bin/samtools" ]]
  else
    conda env list | awk 'NF {print $1"\t"$NF}' | awk -v e="$env" '($1==e){found=1} ($2 ~ "/"e"$"){found=1} END{exit found?0:1}'
  fi
}

create_conda_env_if_needed() {
  local env="$1"
  local frontend="conda"

  if [[ -n "$CONDA_FRONTEND_FORCED" ]]; then
    frontend="$CONDA_FRONTEND_FORCED"
  elif command -v mamba >/dev/null 2>&1; then
    frontend="mamba"
  fi

  command -v "$frontend" >/dev/null 2>&1 || die "Requested conda frontend not found: $frontend"

  if conda_env_exists "$env"; then
    log "Conda env already exists: $env"
    return 0
  fi

  [[ "$CREATE_CONDA_ENV" == true ]] || die "Conda env '$env' does not exist and --no-create-conda-env was used"

  log "Creating conda env: $env using $frontend"
  local -a target_args
  if [[ "$env" == */* ]]; then
    mkdir -p "$(dirname "$(readlink -m "$env")")"
    target_args=(-p "$(readlink -m "$env")")
  else
    target_args=(-n "$env")
  fi

  "$frontend" create -y "${target_args[@]}" \
    -c conda-forge -c bioconda \
    python=3.10 \
    gatk4 \
    samtools \
    htslib \
    openjdk=17
}

activate_conda_env_if_requested() {
  [[ "$USE_CONDA" == true ]] || { log "Conda disabled by --no-conda; using current PATH"; return 0; }

  local conda_sh
  conda_sh="$(find_conda_sh || true)"
  [[ -n "$conda_sh" ]] || die "conda.sh not found. Install conda/mamba or rerun with --no-conda and explicit --gatk-bin."

  # shellcheck disable=SC1090
  source "$conda_sh"

  # mamba shell function may live in a separate file.
  local conda_base
  conda_base="$(conda info --base 2>/dev/null || true)"
  if [[ -n "$conda_base" && -f "$conda_base/etc/profile.d/mamba.sh" ]]; then
    # shellcheck disable=SC1090
    source "$conda_base/etc/profile.d/mamba.sh" || true
  fi

  create_conda_env_if_needed "$CONDA_ENV"

  log "Activating conda env: $CONDA_ENV"
  if [[ "$CONDA_ENV" == */* ]]; then
    conda activate "$(readlink -m "$CONDA_ENV")"
  else
    conda activate "$CONDA_ENV"
  fi
}

activate_conda_env_if_requested

###############################################################################
# Tool resolution.
###############################################################################

if [[ -n "$PROJECT_DIR" && -d "$PROJECT_DIR/.conda" && -z "$GATK_BIN" ]]; then
  found_gatk="$(find "$PROJECT_DIR/.conda" -path '*/bin/gatk' -type f -executable 2>/dev/null | sort | head -n 1 || true)"
  [[ -n "$found_gatk" ]] && GATK_BIN="$found_gatk"
fi

if [[ -n "$GATK_BIN" ]]; then
  GATK_BIN="$(readlink -m "$GATK_BIN")"
  [[ -x "$GATK_BIN" ]] || die "--gatk-bin is not executable: $GATK_BIN"
elif command -v gatk >/dev/null 2>&1; then
  GATK_BIN="$(command -v gatk)"
fi

if [[ -n "$PICARD_JAR" ]]; then
  PICARD_JAR="$(readlink -m "$PICARD_JAR")"
  [[ -s "$PICARD_JAR" ]] || die "--picard-jar not found/empty: $PICARD_JAR"
else
  candidate=""
  candidate="$(find "${HOME}/.nextflow/assets" -path '*/papaemmelab/nf-ffperase/*/assets/picard.jar' -type f 2>/dev/null | sort | head -n 1 || true)"
  if [[ -z "$candidate" ]]; then
    candidate="$(find "${HOME}/.nextflow/assets" -path '*/nf-ffperase*/assets/picard.jar' -type f 2>/dev/null | sort | head -n 1 || true)"
  fi
  [[ -n "$candidate" ]] && PICARD_JAR="$candidate"
fi

if command -v picard >/dev/null 2>&1; then
  PICARD_BIN="$(command -v picard)"
fi

case "$PICARD_BACKEND" in
  auto)
    if [[ -n "$GATK_BIN" ]]; then
      PICARD_BACKEND="gatk"
    elif [[ -n "$PICARD_JAR" ]]; then
      command -v java >/dev/null 2>&1 || die "java not found but picard.jar was detected"
      PICARD_BACKEND="jar"
    elif [[ -n "$PICARD_BIN" ]]; then
      PICARD_BACKEND="picard"
      log "WARNING: falling back to picard wrapper. If it aborts, rerun with default backend/gatk."
    else
      die "No GATK, picard.jar, or picard wrapper found."
    fi
    ;;
  gatk)
    [[ -n "$GATK_BIN" ]] || die "--backend gatk requested, but gatk was not found even after conda setup"
    ;;
  jar)
    [[ -n "$PICARD_JAR" ]] || die "--backend jar requested, but picard.jar was not found. Use --picard-jar."
    command -v java >/dev/null 2>&1 || die "java not found for --backend jar"
    ;;
  picard)
    [[ -n "$PICARD_BIN" ]] || die "--backend picard requested, but picard wrapper was not found"
    log "WARNING: --backend picard is fragile on your server. Prefer default --backend gatk."
    ;;
esac

command -v samtools >/dev/null 2>&1 || die "samtools not found after conda/tool setup"
command -v python3 >/dev/null 2>&1 || die "python3 not found after conda/tool setup"
if [[ "$TIMEOUT_MIN" -gt 0 ]]; then
  command -v timeout >/dev/null 2>&1 || die "timeout command not found but --timeout-min was set"
fi

###############################################################################
# Reference sidecars and optional uppercase reference.
###############################################################################

if [[ ! -s "${REF}.fai" ]]; then
  log "Indexing original reference FASTA with samtools faidx: $REF"
  samtools faidx "$REF"
fi

REF_FOR_PICARD="$REF"
if [[ "$UPPERCASE_REFERENCE" == true ]]; then
  UPPER_DIR="$OUTROOT/reference_uppercase"
  mkdir -p "$UPPER_DIR"
  UPPER_REF="$UPPER_DIR/$(basename "$REF").uppercase.fa"
  if [[ ! -s "$UPPER_REF" ]]; then
    log "Building uppercase reference: $UPPER_REF"
    python3 - "$REF" "$UPPER_REF" <<'PY_UPPER'
import sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, 'rt') as inp, open(dst, 'wt') as out:
    for line in inp:
        if line.startswith('>'):
            out.write(line)
        else:
            out.write(line.upper())
PY_UPPER
  fi
  if [[ ! -s "${UPPER_REF}.fai" ]]; then
    log "Indexing uppercase reference with samtools faidx: $UPPER_REF"
    samtools faidx "$UPPER_REF"
  fi
  REF_FOR_PICARD="$UPPER_REF"
fi

ref_dict_for() {
  local fasta="$1"
  case "$fasta" in
    *.fa) echo "${fasta%.fa}.dict" ;;
    *.fasta) echo "${fasta%.fasta}.dict" ;;
    *) echo "${fasta}.dict" ;;
  esac
}

create_dict_if_missing() {
  local fasta="$1"
  local dict
  dict="$(ref_dict_for "$fasta")"
  if [[ -s "$dict" ]]; then
    return 0
  fi
  log "Creating reference sequence dictionary: $dict"
  if [[ "$PICARD_BACKEND" == "gatk" ]]; then
    "$GATK_BIN" CreateSequenceDictionary -R "$fasta" -O "$dict"
  elif [[ "$PICARD_BACKEND" == "jar" ]]; then
    java -Xmx2g -jar "$PICARD_JAR" CreateSequenceDictionary R="$fasta" O="$dict"
  else
    "$PICARD_BIN" CreateSequenceDictionary R="$fasta" O="$dict"
  fi
}
create_dict_if_missing "$REF_FOR_PICARD"

if [[ -n "$DBSNP" ]]; then
  DBSNP="$(readlink -m "$DBSNP")"
  [[ -s "$DBSNP" ]] || die "--dbsnp not found/empty: $DBSNP"
  if [[ "$DBSNP" == *.gz && ! -s "${DBSNP}.tbi" && ! -s "${DBSNP}.idx" ]]; then
    if [[ -n "$GATK_BIN" ]]; then
      log "Indexing dbSNP with GATK IndexFeatureFile: $DBSNP"
      "$GATK_BIN" IndexFeatureFile -I "$DBSNP"
    else
      log "WARNING: dbSNP has no .tbi/.idx and GATK is unavailable. Picard may fail."
    fi
  fi
fi

log "Script     : $SCRIPT_NAME"
log "Root       : $ROOT"
log "BAM dir    : $BAM_DIR"
log "BAM suffix : $BAM_SUFFIX"
log "Outroot    : $OUTROOT"
log "Reference  : $REF"
log "Picard ref : $REF_FOR_PICARD"
log "Uppercase  : $UPPERCASE_REFERENCE"
log "dbSNP      : ${DBSNP:-none}"
log "Jobs       : $JOBS"
log "Java mem   : $JAVA_MEM per job"
log "MIN_MAPQ   : $MIN_MAPQ"
log "MIN_BASEQ  : $MIN_BASEQ"
log "Timeout    : ${TIMEOUT_MIN} min per sample Picard/GATK job; 0=disabled"
log "Conda      : use=$USE_CONDA env=${CONDA_ENV:-none} create=$CREATE_CONDA_ENV"
log "Backend    : $PICARD_BACKEND"
log "GATK       : ${GATK_BIN:-none}"
log "Picard jar : ${PICARD_JAR:-none}"
log "Picard bin : ${PICARD_BIN:-none}"
log "Force      : $FORCE"
log "Keep going : $KEEP_GOING"

###############################################################################
# Create sample/BAM manifest.
###############################################################################

MANIFEST="$OUTROOT/picard_metrics_manifest.tsv"

python3 - \
  --bam-dir "$BAM_DIR" \
  --bam-suffix "$BAM_SUFFIX" \
  --samples "$SAMPLES_CSV" \
  --samples-file "$SAMPLES_FILE_IN" \
  --samples-from-dir "$SAMPLES_FROM_DIR" \
  --skip-missing "$SKIP_MISSING" \
  --out "$MANIFEST" <<'PY_MANIFEST'
import argparse
import csv
import re
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument('--bam-dir', required=True)
ap.add_argument('--bam-suffix', required=True)
ap.add_argument('--samples', default='')
ap.add_argument('--samples-file', default='')
ap.add_argument('--samples-from-dir', default='')
ap.add_argument('--skip-missing', default='false')
ap.add_argument('--out', required=True)
args = ap.parse_args()

bam_dir = Path(args.bam_dir)
bam_suffix = args.bam_suffix
skip_missing = str(args.skip_missing).lower() == 'true'

suffixes = [
    bam_suffix,
    '.final_variants.vcf.gz', '.final_variants.vcf', '.final_variants',
    '.preprocessed.bam', '.sorted.bam', '.dedup.bam', '.bam',
    '.vcf.gz', '.vcf', '.bam.bai', '.bai', '.tbi',
]

def strip_suffix(x):
    s = str(x).strip()
    s = Path(s).name
    changed = True
    while changed:
        changed = False
        for suf in suffixes:
            if suf and s.endswith(suf):
                s = s[:-len(suf)]
                changed = True
    return s


def add_token(tokens, token):
    token = str(token).strip().strip('"').strip("'")
    if token and not token.startswith('#'):
        tokens.append(token)


tokens = []

if args.samples:
    for part in re.split(r'[,;\n\r\t ]+', args.samples):
        add_token(tokens, part)

if args.samples_file:
    p = Path(args.samples_file)
    if not p.exists():
        raise SystemExit(f'ERROR: --samples-file not found: {p}')
    with p.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            first = re.split(r'[\t, ]+', line)[0]
            if first.lower() in {'sample', 'sample_id', 'id', 'bam', 'path'}:
                continue
            add_token(tokens, first)

if args.samples_from_dir:
    d = Path(args.samples_from_dir)
    if not d.exists() or not d.is_dir():
        raise SystemExit(f'ERROR: --samples-from-dir not found or not a directory: {d}')
    for p in sorted(d.iterdir()):
        add_token(tokens, p.name)

if not tokens:
    for bam in sorted(bam_dir.glob(f'*{bam_suffix}')):
        add_token(tokens, strip_suffix(bam.name))

seen = set()
unique = []
for t in tokens:
    s = strip_suffix(t)
    if not s or s in seen:
        continue
    seen.add(s)
    unique.append(t)

rows = []
missing = []
for token in unique:
    p = Path(token)
    if p.exists() and p.is_file() and p.name.endswith('.bam'):
        sample = strip_suffix(p.name)
        bam = p.resolve()
    else:
        sample = strip_suffix(token)
        candidates = [
            bam_dir / f'{sample}{bam_suffix}',
            bam_dir / f'{sample}.preprocessed.bam',
            bam_dir / f'{sample}.bam',
            bam_dir / sample,
        ]
        bam = None
        for c in candidates:
            if c.exists() and c.is_file():
                bam = c.resolve()
                break
        if bam is None:
            missing.append((sample, str(candidates[0])))
            continue
    rows.append((sample, str(bam)))

if missing and not skip_missing:
    msg = ['ERROR: missing BAMs for selected samples:']
    msg.extend([f'  {s}\t{expected}' for s, expected in missing])
    raise SystemExit('\n'.join(msg))

with open(args.out, 'w', newline='') as out:
    w = csv.writer(out, delimiter='\t', lineterminator='\n')
    w.writerow(['sample', 'bam'])
    w.writerows(rows)

if missing:
    miss_path = str(Path(args.out).with_suffix('.missing.tsv'))
    with open(miss_path, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'expected_bam'])
        w.writerows(missing)

if not rows:
    raise SystemExit('ERROR: no valid samples/BAMs selected')
PY_MANIFEST

[[ -s "$MANIFEST" ]] || die "Manifest not created: $MANIFEST"
log "Manifest   : $MANIFEST"
log "Samples    : $(($(wc -l < "$MANIFEST") - 1))"

if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$MANIFEST" >&2 || cat "$MANIFEST" >&2
else
  cat "$MANIFEST" >&2
fi

if [[ "$DRY_RUN" == true ]]; then
  log "DRY RUN: exiting before Picard/GATK."
  exit 0
fi

###############################################################################
# Embedded Picard table cleaner.
###############################################################################

CONSOLIDATOR="$OUTROOT/clean_picard_metrics_for_ffperase.py"
cat > "$CONSOLIDATOR" <<'PY_CLEAN'
#!/usr/bin/env python3
import argparse
import csv
import glob
from pathlib import Path


def extract_table(path):
    with open(path, 'r', newline='') as fh:
        lines = fh.readlines()
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('SAMPLE_ALIAS\t'):
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError(f'Could not find SAMPLE_ALIAS header in {path}')
    header = lines[header_idx].rstrip('\n\r').split('\t')
    rows = []
    for raw in lines[header_idx + 1:]:
        line = raw.rstrip('\n\r')
        if not line:
            break
        if line.startswith('#') or line.startswith('##'):
            break
        fields = line.split('\t')
        if len(fields) != len(header):
            continue
        if not fields[0] or fields[0] == 'SAMPLE_ALIAS':
            continue
        rows.append(fields)
    return header, rows


def write_merged(pattern, outpath):
    files = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError(f'No files matched {pattern}')
    final_header = None
    final_rows = []
    for f in files:
        header, rows = extract_table(f)
        if final_header is None:
            final_header = header
        elif header != final_header:
            raise RuntimeError(f'Header mismatch in {f}')
        final_rows.extend(rows)
    if not final_rows:
        raise RuntimeError(f'No metric rows found for {pattern}')
    with open(outpath, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(final_header)
        writer.writerows(final_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tmp-picard', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    tmp = Path(args.tmp_picard)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_merged(str(tmp / '*.pre_adapter_detail_metrics'), outdir / 'pre_adapter_metrics.tsv')
    write_merged(str(tmp / '*.bait_bias_detail_metrics'), outdir / 'bait_bias_metrics.tsv')


if __name__ == '__main__':
    main()
PY_CLEAN
chmod +x "$CONSOLIDATOR"

###############################################################################
# Worker functions.
###############################################################################

tail_if_exists() {
  local label="$1"
  local file="$2"
  if [[ -s "$file" ]]; then
    log "---- last 80 lines of $label: $file ----"
    tail -n 80 "$file" >&2 || true
    log "---- end $label ----"
  fi
}

run_picard_cmd() {
  local sample="$1"
  local bam="$2"
  local prefix="$3"
  local stdout_log="$4"
  local stderr_log="$5"
  local tmpdir
  tmpdir="$(dirname "$prefix")"

  local -a cmd
  if [[ "$PICARD_BACKEND" == "gatk" ]]; then
    cmd=("$GATK_BIN" --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${tmpdir}" CollectSequencingArtifactMetrics
      -I "$bam"
      -O "$prefix"
      -R "$REF_FOR_PICARD"
      --MINIMUM_MAPPING_QUALITY "$MIN_MAPQ"
      --MINIMUM_QUALITY_SCORE "$MIN_BASEQ"
      --VALIDATION_STRINGENCY LENIENT)
    [[ -n "$DBSNP" ]] && cmd+=(--DB_SNP "$DBSNP")
  elif [[ "$PICARD_BACKEND" == "jar" ]]; then
    cmd=(java -Xmx"$JAVA_MEM" -Djava.io.tmpdir="$tmpdir" -jar "$PICARD_JAR" CollectSequencingArtifactMetrics
      I="$bam" O="$prefix" R="$REF_FOR_PICARD"
      MINIMUM_MAPPING_QUALITY="$MIN_MAPQ"
      MINIMUM_QUALITY_SCORE="$MIN_BASEQ"
      VALIDATION_STRINGENCY=LENIENT)
    [[ -n "$DBSNP" ]] && cmd+=(DB_SNP="$DBSNP")
  else
    cmd=("$PICARD_BIN" CollectSequencingArtifactMetrics
      I="$bam" O="$prefix" R="$REF_FOR_PICARD"
      MINIMUM_MAPPING_QUALITY="$MIN_MAPQ"
      MINIMUM_QUALITY_SCORE="$MIN_BASEQ"
      VALIDATION_STRINGENCY=LENIENT)
    [[ -n "$DBSNP" ]] && cmd+=(DB_SNP="$DBSNP")
  fi

  log "Picard/GATK command for $sample: ${cmd[*]}"
  export TMPDIR="$tmpdir"

  local status=0
  if [[ "$TIMEOUT_MIN" -gt 0 ]]; then
    timeout --kill-after=60s "${TIMEOUT_MIN}m" "${cmd[@]}" >"$stdout_log" 2>"$stderr_log" || status=$?
  else
    "${cmd[@]}" >"$stdout_log" 2>"$stderr_log" || status=$?
  fi

  if [[ "$status" -ne 0 ]]; then
    log "Picard/GATK failed for sample=$sample status=$status"
    tail_if_exists "stdout" "$stdout_log"
    tail_if_exists "stderr" "$stderr_log"
    return "$status"
  fi
  return 0
}

run_one_sample() {
  local sample="$1"
  local bam="$2"

  local sdir="$OUTROOT/$sample"
  local tmp_base
  if [[ -n "$TMP_ROOT" ]]; then
    tmp_base="$TMP_ROOT/$sample"
  else
    tmp_base="$sdir/tmp"
  fi
  local tmp="$sdir/tmpPicard"
  local prefix="$tmp/picard_${sample}"

  mkdir -p "$sdir/logs" "$tmp_base" "$tmp"
  log "START sample=$sample bam=$bam"

  [[ -s "$bam" ]] || die "$sample: BAM missing/empty: $bam"
  samtools quickcheck "$bam" >"$sdir/logs/${sample}.quickcheck.log" 2>&1 || die "$sample: samtools quickcheck failed; see $sdir/logs/${sample}.quickcheck.log"

  local bai1="${bam}.bai"
  local bai2="${bam%.bam}.bai"
  if [[ ! -s "$bai1" && ! -s "$bai2" ]]; then
    log "Indexing BAM for sample=$sample"
    samtools index -@ "$SAMTOOLS_THREADS" "$bam"
  fi

  if [[ "$FORCE" == true ]]; then
    rm -f "$sdir/pre_adapter_metrics.tsv" "$sdir/bait_bias_metrics.tsv" "$sdir/.picard_metrics.done"
    rm -rf "$tmp"
    mkdir -p "$tmp"
  fi

  if [[ -s "$sdir/pre_adapter_metrics.tsv" && -s "$sdir/bait_bias_metrics.tsv" && "$FORCE" == false ]]; then
    log "SKIP sample=$sample: FFPErase-ready Picard metrics already exist"
    return 0
  fi

  run_picard_cmd \
    "$sample" "$bam" "$prefix" \
    "$sdir/logs/${sample}.picard.stdout.log" \
    "$sdir/logs/${sample}.picard.stderr.log"

  [[ -s "${prefix}.pre_adapter_detail_metrics" ]] || {
    tail_if_exists "Picard stderr" "$sdir/logs/${sample}.picard.stderr.log"
    die "$sample: missing ${prefix}.pre_adapter_detail_metrics"
  }
  [[ -s "${prefix}.bait_bias_detail_metrics" ]] || {
    tail_if_exists "Picard stderr" "$sdir/logs/${sample}.picard.stderr.log"
    die "$sample: missing ${prefix}.bait_bias_detail_metrics"
  }

  log "Cleaning Picard metrics for FFPErase: sample=$sample"
  "$CONSOLIDATOR" \
    --tmp-picard "$tmp" \
    --outdir "$sdir" \
    >"$sdir/logs/${sample}.clean_picard.stdout.log" \
    2>"$sdir/logs/${sample}.clean_picard.stderr.log"

  [[ -s "$sdir/pre_adapter_metrics.tsv" ]] || die "$sample: missing $sdir/pre_adapter_metrics.tsv"
  [[ -s "$sdir/bait_bias_metrics.tsv" ]] || die "$sample: missing $sdir/bait_bias_metrics.tsv"

  date '+%F %T' > "$sdir/.picard_metrics.done"
  log "DONE sample=$sample metrics=$sdir"
}

export -f log die truthy falsy tail_if_exists run_picard_cmd run_one_sample
export ROOT BAM_DIR BAM_SUFFIX OUTROOT REF REF_FOR_PICARD DBSNP JAVA_MEM MIN_MAPQ MIN_BASEQ SAMTOOLS_THREADS
export PICARD_JAR PICARD_BACKEND GATK_BIN PICARD_BIN FORCE TIMEOUT_MIN TMP_ROOT CONSOLIDATOR

###############################################################################
# Controlled launcher with Ctrl+C cleanup.
###############################################################################

terminate_children() {
  local pids
  pids="$(jobs -pr || true)"
  if [[ -n "$pids" ]]; then
    log "Stopping child process groups: $pids"
    while IFS= read -r pid; do
      [[ -n "$pid" ]] || continue
      kill -TERM -- "-$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
    done <<< "$pids"
    sleep 5
    while IFS= read -r pid; do
      [[ -n "$pid" ]] || continue
      kill -KILL -- "-$pid" 2>/dev/null || kill -KILL "$pid" 2>/dev/null || true
    done <<< "$pids"
  fi
}

on_signal() {
  local sig="$1"
  log "Received $sig. Cleaning child jobs."
  trap - INT TERM HUP
  terminate_children
  exit 143
}

trap 'on_signal INT' INT
trap 'on_signal TERM' TERM
trap 'on_signal HUP' HUP

failures=0
launched=0

launch_sample() {
  local sample="$1"
  local bam="$2"
  if command -v setsid >/dev/null 2>&1; then
    setsid bash -c 'run_one_sample "$@"' _ "$sample" "$bam" &
  else
    bash -c 'run_one_sample "$@"' _ "$sample" "$bam" &
  fi
  local pid=$!
  launched=$((launched + 1))
  log "Launched sample=$sample pid=$pid"
}

wait_for_slot() {
  while [[ "$(jobs -pr | wc -l)" -ge "$JOBS" ]]; do
    if ! wait -n; then
      failures=$((failures + 1))
      log "A sample job failed. failures=$failures"
      if [[ "$KEEP_GOING" != true ]]; then
        terminate_children
        exit 1
      fi
    fi
  done
}

while IFS=$'\t' read -r sample bam; do
  [[ "$sample" == "sample" ]] && continue
  [[ -n "${sample:-}" && -n "${bam:-}" ]] || continue
  wait_for_slot
  launch_sample "$sample" "$bam"
done < "$MANIFEST"

while [[ "$(jobs -pr | wc -l)" -gt 0 ]]; do
  if ! wait -n; then
    failures=$((failures + 1))
    log "A sample job failed. failures=$failures"
    if [[ "$KEEP_GOING" != true ]]; then
      terminate_children
      exit 1
    fi
  fi
done

if [[ "$failures" -gt 0 ]]; then
  die "Completed with $failures failed sample job(s). Check per-sample logs under $OUTROOT/<sample>/logs"
fi

log "All requested Picard metrics are ready: $OUTROOT"
find "$OUTROOT" -maxdepth 2 -type f \( -name 'pre_adapter_metrics.tsv' -o -name 'bait_bias_metrics.tsv' \) | sort

###__FFPERASE_EMBEDDED_PICARD_END__###

###__FFPERASE_EMBEDDED_NF_WRAPPER_BEGIN__###
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_VERSION="2026-06-01_parser-v1_model-cache_uppercase-ref_plotfix_precomputed-picard_features-cache_threads_v9_sanitize-snvs"
USE_UPPERCASE_REF=true
SANITIZE_SNVS=true
REF_FFPERASE=""
FFPERASE_MODEL_NAME="ARTIFACT"

ROOT=""
REF=""
OUTROOT=""
MANIFEST=""
SAMPLES_CSV=""
PICARD_METRICS_ROOT=""
ALLOW_MISSING_PICARD_METRICS=false
FEATURES_ROOT=""
FEATURES_CACHE_ROOT=""
REUSE_FEATURES=true
CLASSIFY_ONLY=false
PREPROCESS_ONLY=false
ALLOW_MISSING_FEATURES=false

JOBS=1
MAX_FORKS=2
THREADS=0
TIMEOUT_MIN=360
MIN_MAPQ=20
MIN_BASEQ=20
MIN_DEPTH=1
SPLIT_PILEUP=1000
SPLIT_READS=7500000
VARIANT_WINDOW_BP=500

RUN_SNVS=true
RUN_INDELS=true
FAIL_ON_ERROR=false
FORCE=false
KEEP_WORK=false
USE_RESUME=false
KILL_STALE=false
MODELS_DIR=""
DOWNLOAD_MODELS=true
PREFER_LOCAL_MODEL=true

REPOSITORY="papaemmelab/nf-ffperase"
REVISION="b0dd56cbd0a939896a966b9ce30c4d719b158170"
CONTAINER="auto"
ENGINE="auto"     # auto|apptainer|singularity|docker

CHILD_GROUPS=()
LOG=""

usage() {
  cat <<'EOF'
run_ffperase_post_pipeline.sh

Purpose
-------
Post-process final per-sample VCFs with papaemmelab/nf-ffperase. This is the NF-based wrapper, not the direct annotate_w_pileup wrapper.
It is agnostic to disease/project type. It expects a previous variant pipeline
root containing at least:

  <root>/final_vcf/<sample>.final_variants.vcf.gz
  <root>/preprocessed_bam/<sample>.preprocessed.bam

The script automatically splits each final VCF into SNVs and indels, creates a
small variant-window BED for FFPErase, runs nf-ffperase per sample/type, and
collects classified_df_snvs.tsv / classified_df_indels.tsv into one output tree.

Key fixes in this version
-------------------------
  1. Correct sample detection from files such as SAMPLE_01.final_variants.vcf.gz.
     It strips .final_variants, so it looks for SAMPLE_01.preprocessed.bam, not
     SAMPLE_01.final_variants.preprocessed.bam.

  2. Forces NXF_SYNTAX_PARSER=v1 for nested nf-ffperase. This is required
     because nf-ffperase main still uses legacy Groovy syntax such as
     params."${key}", which breaks under Nextflow 26 strict parser.

  3. Fixes the Hugging Face read-only cache problem inside Apptainer/Singularity
     by passing a local --model file when possible and by forcing HOME,
     XDG_CACHE_HOME, HF_HOME, HUGGINGFACE_HUB_CACHE and TRANSFORMERS_CACHE to
     writable folders under <outroot>/cache/.

  4. Avoids stale Nextflow lock reuse by running every nested nf-ffperase call
     from an isolated launch directory and NOT using -resume by default.

  5. Ctrl+C / TERM / HUP cleanup kills the entire process group for every nested
     Nextflow/Picard/Java task started by this script.

  6. Builds an uppercase copy of the reference for FFPErase by default. This
     avoids nf-ffperase annotate_variants.py failures such as KeyError: gGg
     when the input reference contains soft-masked/lowercase bases. The BAM
     alignments are not changed; only the reference sequence case used by
     FFPErase context extraction is normalized.

  7. Uses --modelName ARTIFACT by default. nf-ffperase PLOT_REPORT expects
     columns ARTIFACT_raw_predicts and ARTIFACT_predicts. Using sample-specific
     modelName values such as FFPErase_snvs_SAMPLE_01 makes classification finish
     but PLOT_REPORT fail with KeyError: ARTIFACT_* not in index.

  8. If nf-ffperase exits non-zero after CLASSIFY_RANDOM_FOREST but the
     classified_df_*.tsv file exists, this wrapper still copies it and marks
     the status as DONE_CLASSIFIED_PLOT_FAILED instead of losing the result.

  9. Deduplicates manifest/sample rows. Some pipeline roots contain both
     <sample>.final_variants.vcf.gz and <sample>.final_variants.vcf. Older
     versions added both to the manifest, so the same sample was classified
     twice. This version prefers the compressed .vcf.gz and runs each sample
     only once per invocation.

  10. Sanitizes SNV VCFs before nf-ffperase. FFPErase/Picard tables only have
      A/C/G/T contexts. Variants with REF/ALT=N or +/-1 reference context=N
      crash annotate_variants.py with KeyError: N. These records are skipped
      and reported instead of aborting the run.

  10. Supports precomputed Picard CollectSequencingArtifactMetrics with
      --picard-metrics-root. For each sample it normalizes/copies the metrics
      into <outroot>/<sample>/input/picard_metrics and passes --picardMetrics
      to nf-ffperase. This avoids the huge genome-wide PICARD scatter step.

  11. Supports FFPErase feature caching. nf-ffperase does not expose a
      --pileupMetrics input, so the supported way to avoid re-running pileup is
      to cache <sample>/<type>/features.tsv after preprocessing and later run
      nf-ffperase with --step classify --features. Use --preprocess-only true
      to generate/cache features, then re-run with --classify-only true or the
      default --reuse-features true.

  12. Adds --threads as a practical alias for nested nf-ffperase maxForks.
      This is the useful parallelism for PILEUP chunks. Picard itself is not
      meaningfully multithreaded per BAM, so use precomputed Picard metrics
      with --picard-metrics-root to avoid internal Picard scatter.

Required flags
--------------
  --root PATH
      Previous variant-calling output root.

  --ref FASTA
      Reference FASTA used for the BAM/VCF.

  --uppercase-reference true|false
      Build/use an uppercase reference copy for FFPErase. Default: true.
      This fixes context-key errors such as KeyError: gGg caused by
      soft-masked/lowercase hg38 sequence.

  --sanitize-snvs true|false
      Remove SNV records that cannot be annotated by FFPErase because REF, ALT
      or the +/-1 reference context is not A/C/G/T. Default: true. Skipped
      variants are written to <outroot>/<sample>/input/*.skipped_for_ffperase.tsv.

Common flags
------------
  --outroot PATH
      Output root for FFPErase post-processing.
      Default: <root>/ffperase_post

  --samples CSV
      Optional comma-separated sample list to keep.
      Example: --samples SAMPLE_01,SAMPLE_02,SAMPLE_03

  --jobs N
      Number of samples/types to run in parallel. Default: 1.
      Recommended for FFPErase on local workstation/server: 1.

  --max-forks N
      Max concurrent nf-ffperase internal tasks per nested run. Default: 2.
      Recommended: 2 for local runs to avoid Picard/IO overload.

  --threads N
      Convenience alias for --max-forks. It controls how many nested
      nf-ffperase tasks can run at once inside each sample/type run. This is
      the useful knob for PILEUP chunk parallelism. It does not make a single
      Picard CollectSequencingArtifactMetrics job multithreaded.

  --timeout-min N
      Timeout per sample/type nf-ffperase run. Default: 360.

  --force
      Delete partial FFPErase sample/type output and rerun. Do not use this
      after a successful run unless you intentionally want to recompute.
      Duplicate manifest/sample rows are still skipped within the same run.

  --ffperase-model-name NAME
      Model name passed to nf-ffperase. Default: ARTIFACT. Do not change unless
      you intentionally want to skip/ignore nf-ffperase PLOT_REPORT, because
      the upstream plot script expects ARTIFACT_raw_predicts and
      ARTIFACT_predicts.

  --keep-work
      Keep nested Nextflow work directories after successful runs.
      By default, successful work directories are removed to reduce disk usage.

  --fail-on-error true|false
      If true, stop/return non-zero on the first failed sample/type.
      Default: false, meaning continue and summarize failures.

  --use-resume true|false
      Use -resume inside each nested nf-ffperase run. Default: false.
      Keep false unless you know the nested launch directory has no stale lock.

Container / engine flags
------------------------
  --container PATH_OR_URI
      Default: auto. Auto uses:
        $APPTAINER_CACHEDIR/nf-ffperase_v1.0.0.sif if present,
        otherwise docker://papaemmelab/nf-ffperase:v1.0.0

  --engine auto|apptainer|singularity|docker
      Default: auto. For your server, apptainer is expected.

  --repository REPO
      Nextflow repo. Default: papaemmelab/nf-ffperase

  --revision REVISION
      nf-ffperase revision. Default: b0dd56cbd0a939896a966b9ce30c4d719b158170

  --models-dir PATH
      Directory with model.snvs.joblib and model.indels.joblib, or where the
      script will try to download/cache them. Default: <outroot>/models

  --download-models true|false
      Try to download missing model files from Hugging Face with curl/wget.
      Default: true. If false and model files are missing, nf-ffperase will try
      its internal DOWNLOAD_MODEL step, using the writable cache/home fix.

  --prefer-local-model true|false
      If a model file exists locally, pass --model to nf-ffperase and bypass
      DOWNLOAD_MODEL. Default: true.

  --picard-metrics-root PATH
      Optional root containing precomputed Picard CollectSequencingArtifactMetrics
      per sample. Expected layouts accepted:

        PATH/<sample>/pre_adapter_metrics.tsv
        PATH/<sample>/bait_bias_metrics.tsv
        PATH/<sample>/tmpPicard/*pre_adapter_detail_metrics
        PATH/<sample>/tmpPicard/*bait_bias_detail_metrics

      The wrapper copies these small files into the per-sample FFPErase input
      folder and passes --picardMetrics to nf-ffperase. This prevents
      nf-ffperase from launching many preprocessWorkflow:PICARD jobs in
      genome-wide mode.

  --allow-missing-picard-metrics true|false
      If --picard-metrics-root is set and a sample lacks valid Picard metrics,
      false means skip/fail that sample/type instead of silently falling back to
      expensive internal Picard. Default: false.

Variant/type flags
------------------
  --snvs-only
      Run SNVs only.

  --indels-only
      Run indels only.

  --skip-snvs
      Do not run SNVs.

  --skip-indels
      Do not run indels.

  --variant-window-bp N
      Window around each variant to build FFPErase BED. Default: 500.

FFPErase parameter flags
------------------------
  --min-mapq N       Default: 20
  --min-baseq N      Default: 20
  --min-depth N      Default: 1
  --split-pileup N   Default: 1000
  --split-reads N    Default: 7500000

Stale process handling
----------------------
  --kill-stale true|false
      If true, kill old processes containing this outroot path before starting.
      Default: false. Use only when you are sure no wanted FFPErase run is active.

Outputs
-------
  <outroot>/input_manifest.normalized.tsv
  <outroot>/<sample>/snvs/classify/classified_df_snvs.tsv
  <outroot>/<sample>/indels/classify/classified_df_indels.tsv
  <outroot>/classified/<sample>.classified_df_snvs.tsv
  <outroot>/classified/<sample>.classified_df_indels.tsv
  <outroot>/features_cache/<sample>/snvs/features.tsv
  <outroot>/features_cache/<sample>/indels/features.tsv
  <outroot>/classified/all_samples.ffperase_snvs.tsv
  <outroot>/classified/all_samples.ffperase_indels.tsv
  <outroot>/status/*.status.tsv
  <outroot>/ffperase_post_pipeline.log

Example
-------
  unset NXF_PLUGINS_DIR
  unset NXF_WORK
  export NXF_HOME="${NXF_HOME:-$HOME/.nextflow}"
  export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$HOME/.apptainer/cache}"
  mkdir -p "$NXF_HOME/plugins" "$APPTAINER_CACHEDIR"

  bash run_ffperase_post_pipeline.sh \
    --root /path/to/variant-calling-output \
    --ref /path/to/reference/genome.fa \
    --picard-metrics-root /path/to/variant-calling-output/ffperase_picard_metrics \
    --reuse-features true \
    --jobs 1 \
    --max-forks 2 \
    --timeout-min 360 \
    --force
EOF
}

as_bool() {
  local v="${1:-}"
  v="$(printf '%s' "$v" | tr '[:upper:]' '[:lower:]')"
  case "$v" in
    true|1|yes|y|on) echo true ;;
    false|0|no|n|off) echo false ;;
    *) echo "ERROR: expected true/false, got: $1" >&2; exit 2 ;;
  esac
}

log() {
  local msg="$*"
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  if [[ -n "${LOG:-}" ]]; then
    echo "[$ts] $msg" | tee -a "$LOG" >&2
  else
    echo "[$ts] $msg" >&2
  fi
}

abs_path() {
  python3 - "$1" <<'PY'
import os, sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PY
}

kill_group() {
  local pgid="$1"
  [[ -n "$pgid" ]] || return 0
  if kill -0 "$pgid" 2>/dev/null || kill -0 "-$pgid" 2>/dev/null; then
    log "Sending TERM to process group $pgid"
    kill -TERM "-$pgid" 2>/dev/null || kill -TERM "$pgid" 2>/dev/null || true
    sleep 8
    if kill -0 "$pgid" 2>/dev/null || kill -0 "-$pgid" 2>/dev/null; then
      log "Sending KILL to process group $pgid"
      kill -KILL "-$pgid" 2>/dev/null || kill -KILL "$pgid" 2>/dev/null || true
    fi
  fi
}

cleanup_all() {
  local st=$?
  trap - EXIT INT TERM HUP
  if [[ ${#CHILD_GROUPS[@]} -gt 0 ]]; then
    log "Cleanup requested; killing ${#CHILD_GROUPS[@]} active child process group(s)."
    local p
    for p in "${CHILD_GROUPS[@]}"; do
      kill_group "$p"
    done
  fi
  exit "$st"
}

trap cleanup_all EXIT
trap 'log "Interrupted"; exit 130' INT
trap 'log "Terminated"; exit 143' TERM HUP

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --outroot|--outdir) OUTROOT="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --samples) SAMPLES_CSV="$2"; shift 2 ;;
    --picard-metrics-root|--picardMetrics-root|--picard-metrics) PICARD_METRICS_ROOT="$2"; shift 2 ;;
    --allow-missing-picard-metrics) ALLOW_MISSING_PICARD_METRICS="$(as_bool "$2")"; shift 2 ;;
    --features-root|--ffperase-features-root) FEATURES_ROOT="$2"; shift 2 ;;
    --features-cache-root) FEATURES_CACHE_ROOT="$2"; shift 2 ;;
    --reuse-features) REUSE_FEATURES="$(as_bool "$2")"; shift 2 ;;
    --preprocess-only) PREPROCESS_ONLY="$(as_bool "$2")"; shift 2 ;;
    --classify-only) CLASSIFY_ONLY="$(as_bool "$2")"; shift 2 ;;
    --allow-missing-features) ALLOW_MISSING_FEATURES="$(as_bool "$2")"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --max-forks) MAX_FORKS="$2"; THREADS="$2"; shift 2 ;;
    --threads) THREADS="$2"; MAX_FORKS="$2"; shift 2 ;;
    --timeout-min) TIMEOUT_MIN="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-baseq|--min-bq) MIN_BASEQ="$2"; shift 2 ;;
    --min-depth) MIN_DEPTH="$2"; shift 2 ;;
    --split-pileup) SPLIT_PILEUP="$2"; shift 2 ;;
    --split-reads) SPLIT_READS="$2"; shift 2 ;;
    --variant-window-bp) VARIANT_WINDOW_BP="$2"; shift 2 ;;
    --repository|--repo) REPOSITORY="$2"; shift 2 ;;
    --revision|-r) REVISION="$2"; shift 2 ;;
    --models-dir) MODELS_DIR="$2"; shift 2 ;;
    --download-models) DOWNLOAD_MODELS="$(as_bool "$2")"; shift 2 ;;
    --prefer-local-model) PREFER_LOCAL_MODEL="$(as_bool "$2")"; shift 2 ;;
    --container) CONTAINER="$2"; shift 2 ;;
    --engine) ENGINE="$2"; shift 2 ;;
    --uppercase-reference|--uppercase-ref) USE_UPPERCASE_REF="$(as_bool "$2")"; shift 2 ;;
    --sanitize-snvs) SANITIZE_SNVS="$(as_bool "$2")"; shift 2 ;;
    --ffperase-model-name|--model-name) FFPERASE_MODEL_NAME="$2"; shift 2 ;;
    --snvs-only) RUN_SNVS=true; RUN_INDELS=false; shift ;;
    --indels-only) RUN_SNVS=false; RUN_INDELS=true; shift ;;
    --skip-snvs) RUN_SNVS=false; shift ;;
    --skip-indels) RUN_INDELS=false; shift ;;
    --force|--overwrite) FORCE=true; shift ;;
    --keep-work) KEEP_WORK=true; shift ;;
    --fail-on-error) FAIL_ON_ERROR="$(as_bool "$2")"; shift 2 ;;
    --use-resume) USE_RESUME="$(as_bool "$2")"; shift 2 ;;
    --kill-stale) KILL_STALE="$(as_bool "$2")"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$ROOT" ]] || { echo "ERROR: --root is required" >&2; usage >&2; exit 2; }
[[ -n "$REF" ]] || { echo "ERROR: --ref is required" >&2; usage >&2; exit 2; }

ROOT="$(abs_path "$ROOT")"
REF="$(abs_path "$REF")"
if [[ -z "$OUTROOT" ]]; then
  OUTROOT="$ROOT/ffperase_classification"
else
  OUTROOT="$(abs_path "$OUTROOT")"
fi
if [[ -n "$MANIFEST" ]]; then
  MANIFEST="$(abs_path "$MANIFEST")"
fi
if [[ -z "$MODELS_DIR" ]]; then
  MODELS_DIR="$OUTROOT/models"
else
  MODELS_DIR="$(abs_path "$MODELS_DIR")"
fi
if [[ -z "$PICARD_METRICS_ROOT" && -d "$ROOT/ffperase_picard_metrics" ]]; then
  PICARD_METRICS_ROOT="$ROOT/ffperase_picard_metrics"
fi
if [[ -n "$PICARD_METRICS_ROOT" ]]; then
  PICARD_METRICS_ROOT="$(abs_path "$PICARD_METRICS_ROOT")"
fi
if [[ -n "$FEATURES_ROOT" ]]; then
  FEATURES_ROOT="$(abs_path "$FEATURES_ROOT")"
fi
if [[ -z "$FEATURES_CACHE_ROOT" ]]; then
  FEATURES_CACHE_ROOT="$OUTROOT/features_cache"
else
  FEATURES_CACHE_ROOT="$(abs_path "$FEATURES_CACHE_ROOT")"
fi
if [[ -z "$FEATURES_ROOT" ]]; then
  FEATURES_ROOT="$FEATURES_CACHE_ROOT"
fi

if [[ "$PREPROCESS_ONLY" == "true" && "$CLASSIFY_ONLY" == "true" ]]; then
  echo "ERROR: --preprocess-only true and --classify-only true are mutually exclusive" >&2
  exit 2
fi

mkdir -p "$OUTROOT" "$OUTROOT"/{classified,status,logs,work,launch,cache,models} "$MODELS_DIR" "$FEATURES_CACHE_ROOT"
LOG="$OUTROOT/ffperase_post_pipeline.log"
: > "$LOG"

log "Script      : run_ffperase_post_pipeline.sh $SCRIPT_VERSION"
log "Root        : $ROOT"
log "Outroot     : $OUTROOT"
log "Reference   : $REF"
log "UppercaseRef: $USE_UPPERCASE_REF"
log "SanitizeSNVs: $SANITIZE_SNVS"
log "Model name  : $FFPERASE_MODEL_NAME"
log "Jobs        : $JOBS"
log "Max forks   : $MAX_FORKS"
log "Threads     : ${THREADS:-$MAX_FORKS} (alias for max-forks / PILEUP chunk parallelism)"
log "Timeout min : $TIMEOUT_MIN"
log "SNVs        : $RUN_SNVS"
log "Indels      : $RUN_INDELS"
log "Models dir  : $MODELS_DIR"
log "Picard root : ${PICARD_METRICS_ROOT:-none}"
log "Picard allow-missing: $ALLOW_MISSING_PICARD_METRICS"
log "Features root: ${FEATURES_ROOT:-none}"
log "Features cache: $FEATURES_CACHE_ROOT"
log "Reuse features: $REUSE_FEATURES"
log "Preprocess only: $PREPROCESS_ONLY"
log "Classify only : $CLASSIFY_ONLY"
log "Parser      : nested nf-ffperase forced with NXF_SYNTAX_PARSER=v1"

[[ -d "$ROOT" ]] || { log "ERROR: root does not exist: $ROOT"; exit 2; }
[[ -s "$REF" ]] || { log "ERROR: reference FASTA does not exist: $REF"; exit 2; }
if [[ -n "$PICARD_METRICS_ROOT" && ! -d "$PICARD_METRICS_ROOT" ]]; then
  log "ERROR: --picard-metrics-root does not exist: $PICARD_METRICS_ROOT"
  exit 2
fi
if [[ -n "$FEATURES_ROOT" && ! -d "$FEATURES_ROOT" ]]; then
  if [[ "$CLASSIFY_ONLY" == "true" ]]; then
    log "ERROR: --features-root does not exist: $FEATURES_ROOT"
    exit 2
  else
    log "WARNING: --features-root does not exist yet: $FEATURES_ROOT"
  fi
fi

prepare_ffperase_reference() {
  local src="$1"
  if [[ "$USE_UPPERCASE_REF" != "true" ]]; then
    REF_FFPERASE="$src"
    return 0
  fi

  local ref_outdir="$OUTROOT/reference_uppercase"
  mkdir -p "$ref_outdir"
  local src_base
  src_base="$(basename "$src")"
  local dst="$ref_outdir/${src_base}.uppercase.fa"
  local tmp="$dst.tmp.$$"

  if [[ ! -s "$dst" || "$FORCE" == "true" ]]; then
    log "Building uppercase FFPErase reference: $dst"
    # Preserve FASTA headers exactly; uppercase sequence lines only. This avoids
    # nf-ffperase annotate_variants.py KeyError on mixed-case trinucleotide
    # contexts such as gGg when hg38 contains soft-masked bases.
    awk 'BEGIN{FS=""} /^>/{print; next} {print toupper($0)}' "$src" > "$tmp"
    mv -f "$tmp" "$dst"
  else
    log "Using existing uppercase FFPErase reference: $dst"
  fi

  [[ -s "$dst.fai" && "$dst" -ot "$dst.fai" ]] || samtools faidx "$dst"

  local dict="${dst%.*}.dict"
  if [[ ! -s "$dict" || "$dst" -nt "$dict" ]]; then
    if samtools dict "$dst" -o "$dict" >/dev/null 2>&1; then
      true
    elif command -v gatk >/dev/null 2>&1; then
      gatk CreateSequenceDictionary -R "$dst" -O "$dict" >/dev/null 2>&1 || true
    else
      log "WARNING: could not create reference dict for $dst; nf-ffperase may still run if not required."
    fi
  fi

  REF_FFPERASE="$dst"
}

for exe in python3 nextflow samtools bcftools bgzip tabix; do
  command -v "$exe" >/dev/null 2>&1 || { log "ERROR: required command not found: $exe"; exit 2; }
done

[[ -s "$REF.fai" ]] || samtools faidx "$REF"
prepare_ffperase_reference "$REF"
log "FFPErase ref: $REF_FFPERASE"

if [[ "$ENGINE" == "auto" ]]; then
  if command -v apptainer >/dev/null 2>&1; then
    ENGINE="apptainer"
  elif command -v singularity >/dev/null 2>&1; then
    ENGINE="singularity"
  elif command -v docker >/dev/null 2>&1; then
    ENGINE="docker"
  else
    log "ERROR: no apptainer, singularity or docker executable found."
    exit 2
  fi
fi
case "$ENGINE" in
  apptainer|singularity|docker) ;;
  *) log "ERROR: --engine must be auto, apptainer, singularity, or docker"; exit 2 ;;
esac

if [[ "$CONTAINER" == "auto" ]]; then
  cache_base="${APPTAINER_CACHEDIR:-${HOME}/.apptainer/cache}"
  local_sif="$cache_base/nf-ffperase_v1.0.0.sif"
  if [[ -s "$local_sif" ]]; then
    CONTAINER="$local_sif"
  else
    CONTAINER="docker://papaemmelab/nf-ffperase:v1.0.0"
  fi
fi

# Protect against the old/private image typo that caused Docker Hub access errors.
if [[ "$CONTAINER" == *"nf-ffperase.0.0"* ]]; then
  log "ERROR: bad container image detected: $CONTAINER"
  log "Use ${APPTAINER_CACHEDIR:-${HOME}/.apptainer/cache}/nf-ffperase_v1.0.0.sif or docker://papaemmelab/nf-ffperase:v1.0.0"
  exit 2
fi

log "Engine      : $ENGINE"
log "Repository  : $REPOSITORY"
log "Revision    : $REVISION"
log "Container   : $CONTAINER"

if [[ "$KILL_STALE" == "true" ]]; then
  log "kill-stale=true; killing old processes containing outroot path: $OUTROOT"
  mapfile -t stale_pids < <(pgrep -af "nf-ffperase|run_ffperase_post_pipeline|FFPErase|picard.jar" | grep -F "$OUTROOT" | awk '{print $1}' | sort -u || true)
  for sp in "${stale_pids[@]:-}"; do
    [[ -n "$sp" ]] || continue
    [[ "$sp" == "$$" ]] && continue
    log "Killing stale PID $sp"
    kill -TERM "$sp" 2>/dev/null || true
  done
  sleep 5
fi

MANIFEST_OUT="$OUTROOT/input_manifest.normalized.tsv"
MANIFEST_SKIP="$OUTROOT/input_manifest.skipped.log"

python3 - "$ROOT" "$MANIFEST" "$MANIFEST_OUT" "$MANIFEST_SKIP" "$SAMPLES_CSV" <<'PY_MANIFEST'
import csv
import os
import sys
from pathlib import Path

root = Path(sys.argv[1])
manifest_in = sys.argv[2]
out = Path(sys.argv[3])
skiplog = Path(sys.argv[4])
samples_csv = sys.argv[5]
keep = None
if samples_csv.strip():
    keep = {x.strip() for x in samples_csv.replace(';', ',').split(',') if x.strip()}

def strip_vcf_suffix(name: str) -> str:
    for suf in ('.vcf.gz', '.vcf'):
        if name.endswith(suf):
            name = name[:-len(suf)]
            break
    for suf in ('.final_variants.snvs', '.final_variants.indels', '.final_variants', '.variants', '.final'):
        if name.endswith(suf):
            name = name[:-len(suf)]
            break
    return name

def find_bam(sample: str):
    bamdir = root / 'preprocessed_bam'
    candidates = [
        bamdir / f'{sample}.preprocessed.bam',
        bamdir / f'{sample}.bam',
        bamdir / f'{sample}.sorted.bam',
    ]
    for c in candidates:
        if c.exists() and c.stat().st_size > 0:
            return c.resolve()
    if bamdir.exists():
        hits = sorted([p for p in bamdir.glob(f'{sample}*.bam') if not str(p).endswith('.bai')])
        if hits:
            return hits[0].resolve()
    return None

def find_bai(bam: Path):
    cands = [Path(str(bam) + '.bai'), bam.with_suffix('.bai')]
    for c in cands:
        if c.exists() and c.stat().st_size > 0:
            return c.resolve()
    return ''

def find_thresholds(sample: str):
    cands = [
        root / 'coverage' / f'{sample}.adaptive_thresholds.tsv',
        root / 'coverage' / f'{sample}.thresholds.tsv',
        root / 'reports' / 'per_sample' / f'{sample}.adaptive_thresholds.tsv',
    ]
    for c in cands:
        if c.exists() and c.stat().st_size > 0:
            return c.resolve()
    return ''

rows = []
skips = []

if manifest_in:
    with open(manifest_in, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            vcf = Path(r.get('vcf') or r.get('final_vcf') or '')
            bam = Path(r.get('bam') or '') if (r.get('bam') or '') else None
            sample = r.get('sample') or (strip_vcf_suffix(vcf.name) if str(vcf) else '')
            sample = strip_vcf_suffix(sample)
            if keep is not None and sample not in keep:
                continue
            if not vcf.exists():
                skips.append(f'SKIP manifest row: sample={sample}; missing VCF {vcf}')
                continue
            if bam is None or not bam.exists():
                bam = find_bam(sample)
            if bam is None:
                skips.append(f'SKIP manifest row: sample={sample}; missing BAM for VCF {vcf}')
                continue
            bai = find_bai(bam)
            rows.append((sample, str(vcf.resolve()), str(bam.resolve()), str(bai), str(find_thresholds(sample))))
else:
    final_dir = root / 'final_vcf'
    vcfs = []
    if final_dir.exists():
        vcfs.extend(sorted(final_dir.glob('*.final_variants.vcf.gz')))
        vcfs.extend(sorted(final_dir.glob('*.final_variants.vcf')))
    if not vcfs and final_dir.exists():
        vcfs.extend(sorted(p for p in final_dir.glob('*.vcf.gz') if '.snvs.' not in p.name and '.indels.' not in p.name))
        vcfs.extend(sorted(p for p in final_dir.glob('*.vcf') if '.snvs.' not in p.name and '.indels.' not in p.name))
    for vcf in vcfs:
        sample = strip_vcf_suffix(vcf.name)
        if keep is not None and sample not in keep:
            continue
        bam = find_bam(sample)
        if bam is None:
            skips.append(f'SKIP manifest pairing: sample={sample}; missing BAM {root / "preprocessed_bam" / (sample + ".preprocessed.bam")}')
            continue
        bai = find_bai(bam)
        rows.append((sample, str(vcf.resolve()), str(bam.resolve()), str(bai), str(find_thresholds(sample))))

# Deduplicate by sample. Prefer compressed final VCFs over uncompressed VCFs,
# and prefer paths containing .final_variants. over generic VCF names. This
# prevents roots with both sample.final_variants.vcf.gz and
# sample.final_variants.vcf from running the same sample twice.
def row_priority(row):
    sample, vcf, bam, bai, thresholds = row
    name = Path(vcf).name
    score = 0
    if name.endswith('.vcf.gz'):
        score += 100
    if '.final_variants.' in name:
        score += 20
    if name.endswith('.final_variants.vcf.gz'):
        score += 10
    return score

best = {}
for row in rows:
    sample = row[0]
    if sample not in best:
        best[sample] = row
    else:
        old = best[sample]
        if row_priority(row) > row_priority(old):
            skips.append(f'DEDUP manifest row: sample={sample}; keeping {row[1]}; dropping duplicate {old[1]}')
            best[sample] = row
        else:
            skips.append(f'DEDUP manifest row: sample={sample}; keeping {old[1]}; dropping duplicate {row[1]}')

rows = [best[k] for k in sorted(best.keys())]

with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='	')
    w.writerow(['sample', 'vcf', 'bam', 'bai', 'thresholds'])
    w.writerows(rows)

with open(skiplog, 'w') as fh:
    for s in skips:
        fh.write(s + '\n')

if not rows:
    print(f'ERROR: no valid sample/VCF/BAM pairs found under {root}', file=sys.stderr)
    if skips:
        print('\n'.join(skips[:50]), file=sys.stderr)
    sys.exit(3)

print(f'Wrote {len(rows)} deduplicated manifest rows to {out}')
if skips:
    print(f'Wrote {len(skips)} skipped/deduplicated rows to {skiplog}')

PY_MANIFEST

while IFS= read -r line; do
  [[ -n "$line" ]] && log "$line"
done < "$MANIFEST_SKIP"
log "Manifest    : $MANIFEST_OUT"
log "Manifest rows: $(($(wc -l < "$MANIFEST_OUT") - 1))"

if [[ -z "${NXF_HOME:-}" ]]; then
  export NXF_HOME="$HOME/.nextflow"
fi
unset NXF_PLUGINS_DIR || true
unset NXF_WORK || true
mkdir -p "$NXF_HOME/plugins"

log "Nextflow    : $(nextflow -version 2>/dev/null | tr '\n' ' ' | sed 's/[[:space:]][[:space:]]*/ /g')"
log "NXF_VER     : ${NXF_VER:-unset}"
log "NXF_HOME    : ${NXF_HOME:-unset}"
log "NXF_SYNTAX  : v1 for nested nf-ffperase commands only"

estimate_coverage() {
  local thresholds="$1"
  local default_cov="3"
  if [[ -n "$thresholds" && -s "$thresholds" ]]; then
    python3 - "$thresholds" "$default_cov" <<'PY_COV'
import csv, sys, math
path, default = sys.argv[1], sys.argv[2]
try:
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        row = next(reader, None)
        if row:
            for key in ('mean_depth', 'coverage', 'mean_coverage', 'depth', 'th_mean_depth'):
                if key in row and row[key] not in ('', '.', None):
                    try:
                        v = float(row[key])
                        if math.isfinite(v) and v > 0:
                            print(max(1, int(round(v))))
                            raise SystemExit
                    except Exception:
                        pass
except Exception:
    pass
print(default)
PY_COV
  else
    echo "$default_cov"
  fi
}

estimate_median_insert() {
  local sample="$1"
  local bam="$2"
  local cache_file="$OUTROOT/cache/${sample}.median_insert.txt"
  if [[ -s "$cache_file" ]]; then
    cat "$cache_file"
    return 0
  fi
  local val=""
  val="$(samtools stats "$bam" 2>/dev/null | awk -F'\t' '$1=="SN" && $2=="insert size average:" {printf "%.0f\n", $3; exit}' || true)"
  if [[ -z "$val" || "$val" == "0" ]]; then
    val="250"
  fi
  echo "$val" > "$cache_file"
  echo "$val"
}

make_variant_window_bed() {
  local vcfgz="$1"
  local bed="$2"
  local window="$3"
  python3 - "$vcfgz" "$bed" "$window" <<'PY_BED'
import gzip
import sys
from pathlib import Path
vcf, bed, window = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
rows = []
open_func = gzip.open if vcf.endswith('.gz') else open
with open_func(vcf, 'rt', errors='replace') as fh:
    for line in fh:
        if not line or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        chrom = f[0]
        pos = int(f[1])
        ref = f[3]
        start = max(0, pos - 1 - window)
        end = pos - 1 + max(1, len(ref)) + window
        rows.append((chrom, start, end))
rows.sort(key=lambda x: (x[0], x[1], x[2]))
merged = []
for chrom, start, end in rows:
    if not merged or merged[-1][0] != chrom or start > merged[-1][2]:
        merged.append([chrom, start, end])
    else:
        if end > merged[-1][2]:
            merged[-1][2] = end
with open(bed, 'w') as out:
    for r in merged:
        out.write(f'{r[0]}\t{r[1]}\t{r[2]}\n')
if not merged:
    raise SystemExit('no intervals')
PY_BED
}


sanitize_snv_vcf_for_ffperase() {
  local in_vcfgz="$1"
  local out_vcfgz="$2"
  local skipped_tsv="$3"
  local ref_fa="$4"
  local tmp_vcf="${out_vcfgz}.sanitized.$$.vcf"
  local stats_tsv="${skipped_tsv%.tsv}.stats.tsv"

  [[ -s "$ref_fa.fai" ]] || samtools faidx "$ref_fa"

  python3 - "$in_vcfgz" "$tmp_vcf" "$skipped_tsv" "$stats_tsv" "$ref_fa" <<'PY_SANITIZE_SNVS'
import gzip
import sys
from pathlib import Path

in_vcf, out_vcf, skipped_tsv, stats_tsv, ref_fa = sys.argv[1:6]
VALID = set("ACGT")

fai = {}
with open(ref_fa + ".fai", "rt") as fh:
    for line in fh:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        name, length, offset, line_bases, line_width = f[:5]
        fai[name] = (int(length), int(offset), int(line_bases), int(line_width))

aliases = {}
for name in list(fai):
    aliases[name] = name
    if name.startswith("chr"):
        aliases[name[3:]] = name
    else:
        aliases["chr" + name] = name

rfh = open(ref_fa, "rb")

def fetch_base(chrom, pos1):
    key = aliases.get(chrom)
    if key is None:
        return "N", "contig_missing"
    length, offset, line_bases, line_width = fai[key]
    if pos1 < 1 or pos1 > length:
        return "N", "context_out_of_bounds"
    zero = pos1 - 1
    byte_offset = offset + (zero // line_bases) * line_width + (zero % line_bases)
    rfh.seek(byte_offset)
    b = rfh.read(1).decode("ascii", errors="ignore").upper()
    if b not in VALID:
        return b if b else "N", "non_acgt_context"
    return b, "ok"

def opener(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, "rt", errors="replace")

kept = skipped = total = 0
reasons = {}
Path(skipped_tsv).parent.mkdir(parents=True, exist_ok=True)
with opener(in_vcf) as inp, open(out_vcf, "wt") as out, open(skipped_tsv, "wt") as sk:
    sk.write("CHROM\tPOS\tID\tREF\tALT\tREASON\tCONTEXT_5\tCONTEXT_REF\tCONTEXT_3\n")
    for line in inp:
        if line.startswith("#"):
            out.write(line)
            continue
        if not line.strip():
            continue
        total += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            chrom = pos = vid = ref = alt = "."
            c5 = cref = c3 = "."
            reason = "malformed_vcf_record"
        else:
            chrom, pos_s, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            ref_u = ref.upper()
            alt_u = alt.upper()
            c5 = cref = c3 = "."
            if "," in alt_u or len(ref_u) != 1 or len(alt_u) != 1:
                reason = "not_biallelic_single_base_snv"
            elif ref_u not in VALID or alt_u not in VALID:
                reason = "non_acgt_ref_or_alt"
            else:
                try:
                    pos_i = int(pos_s)
                except Exception:
                    pos_i = -1
                if pos_i < 1:
                    reason = "invalid_position"
                else:
                    c5, r5 = fetch_base(chrom, pos_i - 1)
                    cref, r0 = fetch_base(chrom, pos_i)
                    c3, r3 = fetch_base(chrom, pos_i + 1)
                    if r5 != "ok" or r0 != "ok" or r3 != "ok":
                        reason = "non_acgt_trinucleotide_context"
                    else:
                        f[3] = ref_u
                        f[4] = alt_u
                        out.write("\t".join(f) + "\n")
                        kept += 1
                        continue
            pos = pos_s if len(f) > 1 else "."
        skipped += 1
        reasons[reason] = reasons.get(reason, 0) + 1
        sk.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, pos, vid, ref, alt, reason, c5, cref, c3))

rfh.close()
with open(stats_tsv, "wt") as st:
    st.write("metric\tvalue\n")
    st.write(f"total_input_snvs\t{total}\n")
    st.write(f"kept_ffperase_safe_snvs\t{kept}\n")
    st.write(f"skipped_snvs\t{skipped}\n")
    for reason, count in sorted(reasons.items()):
        st.write(f"skipped_{reason}\t{count}\n")
PY_SANITIZE_SNVS

  bgzip -f -c "$tmp_vcf" > "$out_vcfgz"
  rm -f "$tmp_vcf"
}

write_nf_config() {
  local cfg="$1"
  local sample_cache_root="$2"
  local bind_root="$3"
  local ref_dir="$4"
  local container_for_cfg="$5"
  local engine_for_cfg="$6"

  local home_cache="$sample_cache_root/home"
  local xdg_cache="$sample_cache_root/xdg"
  local hf_home="$sample_cache_root/huggingface"
  local hf_hub="$hf_home/hub"
  local transformers="$hf_home/transformers"
  local mpl="$sample_cache_root/matplotlib"
  mkdir -p "$home_cache" "$xdg_cache" "$hf_home" "$hf_hub" "$transformers" "$mpl"

  local nxf_home_bind="${NXF_HOME:-$HOME/.nextflow}"
  mkdir -p "$nxf_home_bind" || true

  # Bind everything the nested nf-ffperase task needs. --home + --env are
  # intentional: the older FFPErase image/HuggingFace stack ignores HF_HOME in
  # some cases and falls back to ~/.cache. Making HOME writable avoids this.
  local bind_opts="--bind ${sample_cache_root}:${sample_cache_root} --bind ${bind_root}:${bind_root} --bind ${ref_dir}:${ref_dir} --bind ${nxf_home_bind}:${nxf_home_bind} --home ${home_cache}:${home_cache} --env HOME=${home_cache} --env XDG_CACHE_HOME=${xdg_cache} --env HF_HOME=${hf_home} --env HUGGINGFACE_HUB_CACHE=${hf_hub} --env TRANSFORMERS_CACHE=${transformers} --env MPLCONFIGDIR=${mpl} --env LC_ALL=C.UTF-8 --env LANG=C.UTF-8"

  # beforeScript runs inside the task context and is a second layer of defense
  # for commands that read HOME/HF_HOME directly.
  local before_script="export HOME='${home_cache}'; export XDG_CACHE_HOME='${xdg_cache}'; export HF_HOME='${hf_home}'; export HUGGINGFACE_HUB_CACHE='${hf_hub}'; export TRANSFORMERS_CACHE='${transformers}'; export MPLCONFIGDIR='${mpl}'; export LC_ALL='C.UTF-8'; export LANG='C.UTF-8'"

  if [[ "$engine_for_cfg" == "docker" ]]; then
    local docker_container="$container_for_cfg"
    docker_container="${docker_container#docker://}"
    cat > "$cfg" <<EOF_CFG
process {
  executor = 'local'
  cpus = 1
  maxForks = ${MAX_FORKS}
  container = '${docker_container}'
  beforeScript = '''${before_script}'''
}

env {
  HOME = '${home_cache}'
  XDG_CACHE_HOME = '${xdg_cache}'
  HF_HOME = '${hf_home}'
  HUGGINGFACE_HUB_CACHE = '${hf_hub}'
  TRANSFORMERS_CACHE = '${transformers}'
  MPLCONFIGDIR = '${mpl}'
  LC_ALL = 'C.UTF-8'
  LANG = 'C.UTF-8'
}

docker {
  enabled = true
  runOptions = '-v ${sample_cache_root}:${sample_cache_root} -v ${bind_root}:${bind_root} -v ${ref_dir}:${ref_dir} -v ${nxf_home_bind}:${nxf_home_bind} -e HOME=${home_cache} -e XDG_CACHE_HOME=${xdg_cache} -e HF_HOME=${hf_home} -e HUGGINGFACE_HUB_CACHE=${hf_hub} -e TRANSFORMERS_CACHE=${transformers} -e MPLCONFIGDIR=${mpl}'
}

apptainer.enabled = false
singularity.enabled = false
EOF_CFG
  elif [[ "$engine_for_cfg" == "singularity" ]]; then
    cat > "$cfg" <<EOF_CFG
process {
  executor = 'local'
  cpus = 1
  maxForks = ${MAX_FORKS}
  container = '${container_for_cfg}'
  beforeScript = '''${before_script}'''
}

env {
  HOME = '${home_cache}'
  XDG_CACHE_HOME = '${xdg_cache}'
  HF_HOME = '${hf_home}'
  HUGGINGFACE_HUB_CACHE = '${hf_hub}'
  TRANSFORMERS_CACHE = '${transformers}'
  MPLCONFIGDIR = '${mpl}'
  LC_ALL = 'C.UTF-8'
  LANG = 'C.UTF-8'
}

singularity {
  enabled = true
  autoMounts = true
  runOptions = '${bind_opts}'
}

apptainer.enabled = false
docker.enabled = false
EOF_CFG
  else
    cat > "$cfg" <<EOF_CFG
process {
  executor = 'local'
  cpus = 1
  maxForks = ${MAX_FORKS}
  container = '${container_for_cfg}'
  beforeScript = '''${before_script}'''
}

env {
  HOME = '${home_cache}'
  XDG_CACHE_HOME = '${xdg_cache}'
  HF_HOME = '${hf_home}'
  HUGGINGFACE_HUB_CACHE = '${hf_hub}'
  TRANSFORMERS_CACHE = '${transformers}'
  MPLCONFIGDIR = '${mpl}'
  LC_ALL = 'C.UTF-8'
  LANG = 'C.UTF-8'
}

apptainer {
  enabled = true
  autoMounts = true
  runOptions = '${bind_opts}'
}

singularity.enabled = false
docker.enabled = false
EOF_CFG
  fi
}

ensure_ffperase_model() {
  local mtype="$1"
  local model_file="$MODELS_DIR/model.${mtype}.joblib"

  if [[ -s "$model_file" ]]; then
    printf '%s\n' "$model_file"
    return 0
  fi

  if [[ "$DOWNLOAD_MODELS" != "true" ]]; then
    return 1
  fi

  mkdir -p "$MODELS_DIR"
  local url="https://huggingface.co/papaemmelab/ffperase/resolve/main/model.${mtype}.joblib"
  local tmp="${model_file}.tmp.$$"
  log "Model missing for ${mtype}; trying to download to $model_file"

  if command -v curl >/dev/null 2>&1; then
    if curl -L --fail --retry 3 --connect-timeout 20 -o "$tmp" "$url"; then
      mv -f "$tmp" "$model_file"
      printf '%s\n' "$model_file"
      return 0
    fi
  elif command -v wget >/dev/null 2>&1; then
    if wget -O "$tmp" "$url"; then
      mv -f "$tmp" "$model_file"
      printf '%s\n' "$model_file"
      return 0
    fi
  fi

  rm -f "$tmp"
  log "WARNING: could not download model.${mtype}.joblib. nf-ffperase will try its internal DOWNLOAD_MODEL step."
  return 1
}

wait_with_timeout() {
  local pid="$1"
  local timeout_min="$2"
  local start now elapsed
  start="$(date +%s)"
  while kill -0 "$pid" 2>/dev/null; do
    sleep 5
    now="$(date +%s)"
    elapsed=$(( (now - start) / 60 ))
    if (( timeout_min > 0 && elapsed >= timeout_min )); then
      log "Timeout after ${timeout_min} min for process group $pid"
      kill_group "$pid"
      wait "$pid" 2>/dev/null || true
      return 124
    fi
  done
  wait "$pid"
}

find_picard_metrics_source_dir() {
  local sample="$1"
  [[ -n "$PICARD_METRICS_ROOT" ]] || return 1

  local candidates=(
    "$PICARD_METRICS_ROOT/$sample"
    "$PICARD_METRICS_ROOT/${sample}.final_variants"
    "$PICARD_METRICS_ROOT"
  )

  local d
  for d in "${candidates[@]}"; do
    [[ -d "$d" ]] || continue
    local pre_final="$d/pre_adapter_metrics.tsv"
    local bait_final="$d/bait_bias_metrics.tsv"
    local pre_raw=""
    local bait_raw=""
    pre_raw="$(find "$d" -maxdepth 2 -type f -name '*pre_adapter_detail_metrics*' 2>/dev/null | head -n 1 || true)"
    bait_raw="$(find "$d" -maxdepth 2 -type f -name '*bait_bias_detail_metrics*' 2>/dev/null | head -n 1 || true)"

    if [[ -s "$pre_final" && -s "$bait_final" ]]; then
      printf '%s\n' "$d"
      return 0
    fi
    if [[ -n "$pre_raw" && -s "$pre_raw" && -n "$bait_raw" && -s "$bait_raw" ]]; then
      printf '%s\n' "$d"
      return 0
    fi
  done

  return 1
}

prepare_picard_metrics_for_sample() {
  local sample="$1"
  local input_dir="$2"
  [[ -n "$PICARD_METRICS_ROOT" ]] || return 1

  local src=""
  if ! src="$(find_picard_metrics_source_dir "$sample")"; then
    log "WARNING: no valid precomputed Picard metrics found for sample=$sample under $PICARD_METRICS_ROOT"
    return 2
  fi

  local dst="$input_dir/picard_metrics"
  rm -rf "$dst"
  mkdir -p "$dst" "$dst/tmpPicard"

  # Copy final FFPErase-ready TSVs when present.
  [[ -s "$src/pre_adapter_metrics.tsv" ]] && cp -L -f "$src/pre_adapter_metrics.tsv" "$dst/pre_adapter_metrics.tsv"
  [[ -s "$src/bait_bias_metrics.tsv" ]] && cp -L -f "$src/bait_bias_metrics.tsv" "$dst/bait_bias_metrics.tsv"

  # Also copy raw Picard detail metrics into tmpPicard. nf-ffperase
  # collect_picard.py aggregates raw detail metrics only when tmpPicard exists.
  find "$src" -maxdepth 2 -type f \
    \( -name '*pre_adapter_detail_metrics*' -o -name '*bait_bias_detail_metrics*' \) \
    -exec cp -L -f {} "$dst/tmpPicard/" \; 2>/dev/null || true

  local have_final=false
  local have_raw=false
  if [[ -s "$dst/pre_adapter_metrics.tsv" && -s "$dst/bait_bias_metrics.tsv" ]]; then
    have_final=true
  fi
  if find "$dst/tmpPicard" -maxdepth 1 -type f -name '*pre_adapter_detail_metrics*' | grep -q . \
     && find "$dst/tmpPicard" -maxdepth 1 -type f -name '*bait_bias_detail_metrics*' | grep -q .; then
    have_raw=true
  fi

  if [[ "$have_final" != "true" && "$have_raw" != "true" ]]; then
    log "WARNING: Picard metrics source was found but normalized folder is incomplete for sample=$sample: $src"
    return 2
  fi

  log "Using precomputed Picard metrics: sample=$sample src=$src input=$dst final_tsv=$have_final raw_detail=$have_raw"
  printf '%s\n' "$dst"
  return 0
}


find_features_source_file() {
  local sample="$1"
  local mtype="$2"
  [[ -n "$FEATURES_ROOT" ]] || return 1

  local candidates=(
    "$FEATURES_ROOT/$sample/$mtype/features.tsv"
    "$FEATURES_ROOT/$sample/$mtype/preprocess/features.tsv"
    "$FEATURES_ROOT/$sample/${mtype}/features.tsv"
    "$FEATURES_ROOT/$sample/${mtype}/preprocess/features.tsv"
    "$FEATURES_ROOT/$sample.${mtype}.features.tsv"
    "$FEATURES_ROOT/${sample}.classified_df_${mtype}.features.tsv"
    "$FEATURES_ROOT/${sample}/features.${mtype}.tsv"
    "$FEATURES_ROOT/${sample}_${mtype}.features.tsv"
    "$OUTROOT/features_cache/$sample/$mtype/features.tsv"
    "$OUTROOT/$sample/$mtype/preprocess/features.tsv"
  )

  local f
  for f in "${candidates[@]}"; do
    if [[ -s "$f" ]]; then
      printf '%s\n' "$f"
      return 0
    fi
  done

  return 1
}

prepare_features_for_sample_type() {
  local sample="$1"
  local mtype="$2"
  local input_dir="$3"
  local src=""

  if ! src="$(find_features_source_file "$sample" "$mtype")"; then
    return 1
  fi

  local dst="$input_dir/${sample}.${mtype}.features.tsv"
  cp -L -f "$src" "$dst"
  [[ -s "$dst" ]] || return 1
  log "Using precomputed FFPErase features: sample=$sample type=$mtype src=$src input=$dst"
  printf '%s\n' "$dst"
  return 0
}

cache_features_after_preprocess() {
  local sample="$1"
  local mtype="$2"
  local type_out="$3"
  local feature_src="$type_out/preprocess/features.tsv"
  local cache_dir="$FEATURES_CACHE_ROOT/$sample/$mtype"
  local cache_file="$cache_dir/features.tsv"

  if [[ -s "$feature_src" ]]; then
    mkdir -p "$cache_dir"
    cp -f "$feature_src" "$cache_file"
    log "Cached FFPErase features: sample=$sample type=$mtype cache=$cache_file"
    return 0
  fi

  return 1
}

count_feature_rows() {
  local f="$1"
  if [[ -s "$f" ]]; then
    awk 'NR>1 {n++} END {print n+0}' "$f"
  else
    echo 0
  fi
}

run_one() {
  local sample="$1"
  local final_vcf="$2"
  local bam="$3"
  local bai="$4"
  local thresholds="$5"
  local mtype="$6"

  local sample_dir="$OUTROOT/$sample"
  local input_dir="$sample_dir/input"
  local type_out="$sample_dir/$mtype"
  local workdir="$OUTROOT/work/$sample/$mtype"
  local launch_dir="$OUTROOT/launch/$sample/$mtype"
  local sample_cache="$OUTROOT/cache/$sample/$mtype"
  local console="$type_out/${sample}.nf_ffperase_${mtype}.console.log"
  local status_file="$OUTROOT/status/${sample}.${mtype}.status.tsv"
  local ref_dir
  ref_dir="$(dirname "$REF_FFPERASE")"

  mkdir -p "$input_dir" "$type_out" "$workdir" "$launch_dir" "$sample_cache" "$OUTROOT/classified" "$FEATURES_CACHE_ROOT"

  if [[ "$PREPROCESS_ONLY" == "true" && "$FORCE" != "true" ]]; then
    local existing_feature=""
    if existing_feature="$(find_features_source_file "$sample" "$mtype")"; then
      local cache_dir="$FEATURES_CACHE_ROOT/$sample/$mtype"
      local cache_file="$cache_dir/features.tsv"
      mkdir -p "$cache_dir"
      if [[ "$existing_feature" != "$cache_file" ]]; then
        cp -L -f "$existing_feature" "$cache_file"
      fi
      log "SKIP existing FFPErase features: sample=$sample type=$mtype features=$cache_file"
      printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "SKIP_EXISTING_FEATURES" "." "$cache_file" > "$OUTROOT/status/${sample}.${mtype}.status.tsv"
      return 0
    fi
  fi

  local classified_canonical="$OUTROOT/classified/${sample}.classified_df_${mtype}.tsv"
  local classified_nested="$type_out/classify/classified_df_${mtype}.tsv"

  if [[ "$FORCE" == "true" ]]; then
    if [[ "$CLASSIFY_ONLY" == "true" ]]; then
      # In classify-only mode, keep preprocess logs/features; only clear classify/work.
      rm -rf "$type_out/classify" "$workdir" "$launch_dir"
      mkdir -p "$type_out" "$type_out/classify" "$workdir" "$launch_dir"
    else
      rm -rf "$type_out" "$workdir" "$launch_dir"
      mkdir -p "$type_out" "$workdir" "$launch_dir"
    fi
  elif [[ "$PREPROCESS_ONLY" != "true" && -s "$classified_nested" ]]; then
    log "SKIP existing classified output: sample=$sample type=$mtype"
    cp -f "$classified_nested" "$classified_canonical"
    printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "SKIP_EXISTING" "." "nested classified output exists" > "$status_file"
    return 0
  elif [[ "$PREPROCESS_ONLY" != "true" && -s "$classified_canonical" ]]; then
    log "SKIP existing collected classified output: sample=$sample type=$mtype"
    printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "SKIP_EXISTING_COLLECTED" "." "collected classified output exists" > "$status_file"
    return 0
  elif [[ "$CLASSIFY_ONLY" == "true" ]]; then
    # Do not delete preprocess output/logs while checking cached features.
    rm -rf "$type_out/classify" "$workdir" "$launch_dir"
    mkdir -p "$type_out" "$type_out/classify" "$workdir" "$launch_dir"
  else
    # Partial failed output should not keep stale .nextflow locks or partial classify dirs.
    rm -rf "$type_out" "$workdir" "$launch_dir"
    mkdir -p "$type_out" "$workdir" "$launch_dir"
  fi

  local bam_in="$input_dir/${sample}.bam"
  ln -sf "$bam" "$bam_in"
  if [[ -n "$bai" && -s "$bai" ]]; then
    ln -sf "$bai" "$bam_in.bai"
  else
    log "Indexing BAM for sample=$sample because .bai was missing"
    samtools index "$bam_in"
  fi

  local split_vcf="$input_dir/${sample}.final_variants.${mtype}.vcf.gz"
  local split_tmp="$split_vcf.tmp.vcf.gz"
  if [[ "$mtype" == "snvs" ]]; then
    if [[ "$SANITIZE_SNVS" == "true" ]]; then
      local raw_snvs="$split_vcf.raw.tmp.vcf.gz"
      local skipped_snvs="$input_dir/${sample}.snvs.skipped_for_ffperase.tsv"
      bcftools view -v snps -Oz -o "$raw_snvs" "$final_vcf"
      tabix -f -p vcf "$raw_snvs" || true
      sanitize_snv_vcf_for_ffperase "$raw_snvs" "$split_tmp" "$skipped_snvs" "$REF_FFPERASE"
      rm -f "$raw_snvs" "$raw_snvs.tbi"
      local skipped_count=0 kept_count=0
      skipped_count="$(awk 'NR>1{n++} END{print n+0}' "$skipped_snvs" 2>/dev/null || echo 0)"
      kept_count="$(bcftools view -H "$split_tmp" 2>/dev/null | wc -l | awk '{print $1}')"
      log "SNV sanitize: sample=$sample kept=$kept_count skipped=$skipped_count skipped_report=$skipped_snvs"
    else
      bcftools view -v snps -Oz -o "$split_tmp" "$final_vcf"
    fi
  else
    bcftools view -v indels -Oz -o "$split_tmp" "$final_vcf"
  fi
  mv -f "$split_tmp" "$split_vcf"
  tabix -f -p vcf "$split_vcf"

  local nvar
  nvar="$(bcftools view -H "$split_vcf" | wc -l | awk '{print $1}')"
  if [[ "$nvar" == "0" ]]; then
    log "SKIP no variants: sample=$sample type=$mtype"
    printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "SKIP_NO_VARIANTS" "0" "no variants after split" > "$status_file"
    return 0
  fi

  local bed="$input_dir/${sample}.${mtype}.ffperase.regions.bed"
  make_variant_window_bed "$split_vcf" "$bed" "$VARIANT_WINDOW_BP"

  local coverage median_insert
  coverage="$(estimate_coverage "$thresholds")"
  median_insert="$(estimate_median_insert "$sample" "$bam")"

  local cfg="$sample_dir/nf_ffperase_LOCAL.config"
  write_nf_config "$cfg" "$sample_cache" "$OUTROOT" "$ref_dir" "$CONTAINER" "$ENGINE"

  local cache_home="$sample_cache/home"
  local xdg_cache="$sample_cache/xdg"
  local hf_home="$sample_cache/huggingface"
  local hf_hub="$hf_home/hub"
  local transformers="$hf_home/transformers"
  local mpl="$sample_cache/matplotlib"
  mkdir -p "$cache_home" "$xdg_cache" "$hf_home" "$hf_hub" "$transformers" "$mpl"

  export HOME="$cache_home"
  export XDG_CACHE_HOME="$xdg_cache"
  export HF_HOME="$hf_home"
  export HUGGINGFACE_HUB_CACHE="$hf_hub"
  export TRANSFORMERS_CACHE="$transformers"
  export MPLCONFIGDIR="$mpl"
  export APPTAINERENV_XDG_CACHE_HOME="$xdg_cache"
  export APPTAINERENV_HF_HOME="$hf_home"
  export APPTAINERENV_HUGGINGFACE_HUB_CACHE="$hf_hub"
  export APPTAINERENV_TRANSFORMERS_CACHE="$transformers"
  export APPTAINERENV_MPLCONFIGDIR="$mpl"
  export APPTAINERENV_LC_ALL="C.UTF-8"
  export APPTAINERENV_LANG="C.UTF-8"
  export SINGULARITYENV_XDG_CACHE_HOME="$xdg_cache"
  export SINGULARITYENV_HF_HOME="$hf_home"
  export SINGULARITYENV_HUGGINGFACE_HUB_CACHE="$hf_hub"
  export SINGULARITYENV_TRANSFORMERS_CACHE="$transformers"
  export SINGULARITYENV_MPLCONFIGDIR="$mpl"
  export SINGULARITYENV_LC_ALL="C.UTF-8"
  export SINGULARITYENV_LANG="C.UTF-8"

  local resume_args=()
  if [[ "$USE_RESUME" == "true" ]]; then
    resume_args=(-resume)
  fi

  local model_args=()
  local local_model=""
  local run_step="full"
  local features_in=""
  local feature_rows=""

  # There is no nf-ffperase --pileupMetrics parameter. The reusable output of
  # preprocessing is features.tsv. Therefore, if features exist, run classify
  # only and avoid PILEUP/PICARD completely.
  if [[ "$CLASSIFY_ONLY" == "true" || "$REUSE_FEATURES" == "true" ]]; then
    if features_in="$(prepare_features_for_sample_type "$sample" "$mtype" "$input_dir")"; then
      run_step="classify"
      feature_rows="$(count_feature_rows "$features_in")"
      if [[ "$feature_rows" != "0" ]]; then
        nvar="$feature_rows"
      fi
    elif [[ "$CLASSIFY_ONLY" == "true" ]]; then
      local msg="missing precomputed FFPErase features for sample=$sample type=$mtype under ${FEATURES_ROOT}"
      printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "FAILED_MISSING_FEATURES" "$nvar" "$msg" > "$status_file"
      if [[ "$ALLOW_MISSING_FEATURES" == "true" ]]; then
        log "SKIP sample=$sample type=$mtype because $msg"
        return 0
      fi
      log "FAIL sample=$sample type=$mtype because $msg"
      if [[ "$FAIL_ON_ERROR" == "true" ]]; then
        return 2
      fi
      return 0
    fi
  fi

  if [[ "$PREPROCESS_ONLY" == "true" ]]; then
    run_step="preprocess"
  fi

  if [[ "$run_step" != "preprocess" ]]; then
    if [[ "$PREFER_LOCAL_MODEL" == "true" ]]; then
      if local_model="$(ensure_ffperase_model "$mtype")"; then
        model_args=(--model "$local_model")
        log "Using local FFPErase model: $local_model"
      fi
    fi
  fi

  local picard_metrics_dir=""
  local picard_args=()
  if [[ "$run_step" != "classify" && -n "$PICARD_METRICS_ROOT" ]]; then
    local pm_st=0
    set +e
    picard_metrics_dir="$(prepare_picard_metrics_for_sample "$sample" "$input_dir")"
    pm_st=$?
    set -e
    if [[ "$pm_st" -ne 0 || -z "$picard_metrics_dir" ]]; then
      local msg="missing or invalid precomputed Picard metrics under ${PICARD_METRICS_ROOT}"
      printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "FAILED_MISSING_PICARD_METRICS" "$nvar" "$msg" > "$status_file"
      if [[ "$ALLOW_MISSING_PICARD_METRICS" == "true" ]]; then
        log "WARNING: $msg; falling back to nf-ffperase internal Picard for sample=$sample type=$mtype"
      else
        log "SKIP/FAIL sample=$sample type=$mtype because $msg. Set --allow-missing-picard-metrics true to allow slow fallback."
        if [[ "$FAIL_ON_ERROR" == "true" ]]; then
          return 2
        fi
        return 0
      fi
    else
      picard_args=(--picardMetrics "$picard_metrics_dir")
    fi
  fi

  local cmd=()
  if [[ "$run_step" == "classify" ]]; then
    cmd=(nextflow run "$REPOSITORY" -r "$REVISION" -c "$cfg" "${resume_args[@]}" -work-dir "$workdir"
      --step classify
      --features "$features_in"
      --outdir "$type_out"
      --outdirClassify "$type_out/classify"
      --outdirTrain "$type_out/train"
      --mutationType "$mtype"
      --modelName "$FFPERASE_MODEL_NAME"
      "${model_args[@]}")
    log "Running FFPERASE classify-only: sample=$sample type=$mtype variants=$nvar features=$features_in"
  elif [[ "$run_step" == "preprocess" ]]; then
    cmd=(nextflow run "$REPOSITORY" -r "$REVISION" -c "$cfg" "${resume_args[@]}" -work-dir "$workdir"
      --step preprocess
      --vcf "$split_vcf"
      --bam "$bam_in"
      --reference "$REF_FFPERASE"
      --bed "$bed"
      "${picard_args[@]}"
      --outdir "$type_out"
      --outdirPreprocess "$type_out/preprocess"
      --outdirClassify "$type_out/classify"
      --outdirTrain "$type_out/train"
      --coverage "$coverage"
      --medianInsert "$median_insert"
      --mutationType "$mtype"
      --modelName "$FFPERASE_MODEL_NAME"
      --minMapq "$MIN_MAPQ"
      --minBaseq "$MIN_BASEQ"
      --minDepth "$MIN_DEPTH"
      --splitPileup "$SPLIT_PILEUP"
      --splitReads "$SPLIT_READS")
    log "Running FFPERASE preprocess-only: sample=$sample type=$mtype variants=$nvar coverage=$coverage medianInsert=$median_insert picardMetrics=${picard_metrics_dir:-none}"
  else
    cmd=(nextflow run "$REPOSITORY" -r "$REVISION" -c "$cfg" "${resume_args[@]}" -work-dir "$workdir"
      --step full
      --vcf "$split_vcf"
      --bam "$bam_in"
      --reference "$REF_FFPERASE"
      --bed "$bed"
      "${picard_args[@]}"
      --outdir "$type_out"
      --outdirPreprocess "$type_out/preprocess"
      --outdirClassify "$type_out/classify"
      --outdirTrain "$type_out/train"
      --coverage "$coverage"
      --medianInsert "$median_insert"
      --mutationType "$mtype"
      --modelName "$FFPERASE_MODEL_NAME"
      "${model_args[@]}"
      --minMapq "$MIN_MAPQ"
      --minBaseq "$MIN_BASEQ"
      --minDepth "$MIN_DEPTH"
      --splitPileup "$SPLIT_PILEUP"
      --splitReads "$SPLIT_READS")
    log "Running FFPERASE full: sample=$sample type=$mtype variants=$nvar coverage=$coverage medianInsert=$median_insert picardMetrics=${picard_metrics_dir:-none}"
  fi

  log "Log: $console"
  log "Launch dir: $launch_dir"

  local nf_pid=""
  set +e
  setsid env \
    NXF_SYNTAX_PARSER=v1 \
    NXF_HOME="${NXF_HOME:-${HOME}/.nextflow}" \
    HOME="$cache_home" \
    XDG_CACHE_HOME="$xdg_cache" \
    HF_HOME="$hf_home" \
    HUGGINGFACE_HUB_CACHE="$hf_hub" \
    TRANSFORMERS_CACHE="$transformers" \
    MPLCONFIGDIR="$mpl" \
    APPTAINERENV_XDG_CACHE_HOME="$xdg_cache" \
    APPTAINERENV_HF_HOME="$hf_home" \
    APPTAINERENV_HUGGINGFACE_HUB_CACHE="$hf_hub" \
    APPTAINERENV_TRANSFORMERS_CACHE="$transformers" \
    APPTAINERENV_MPLCONFIGDIR="$mpl" \
    APPTAINERENV_LC_ALL="C.UTF-8" \
    APPTAINERENV_LANG="C.UTF-8" \
    SINGULARITYENV_XDG_CACHE_HOME="$xdg_cache" \
    SINGULARITYENV_HF_HOME="$hf_home" \
    SINGULARITYENV_HUGGINGFACE_HUB_CACHE="$hf_hub" \
    SINGULARITYENV_TRANSFORMERS_CACHE="$transformers" \
    SINGULARITYENV_MPLCONFIGDIR="$mpl" \
    SINGULARITYENV_LC_ALL="C.UTF-8" \
    SINGULARITYENV_LANG="C.UTF-8" \
    bash -c 'cd "$1" || exit 99; shift; exec "$@"' _ "$launch_dir" "${cmd[@]}" > "$console" 2>&1 &
  nf_pid=$!
  CHILD_GROUPS+=("$nf_pid")
  wait_with_timeout "$nf_pid" "$TIMEOUT_MIN"
  local st=$?
  set -e

  # Remove this process group from the active cleanup list.
  local new_groups=()
  local g
  for g in "${CHILD_GROUPS[@]}"; do
    [[ "$g" != "$nf_pid" ]] && new_groups+=("$g")
  done
  CHILD_GROUPS=("${new_groups[@]}")

  local classified="$type_out/classify/classified_df_${mtype}.tsv"

  if [[ "$run_step" == "preprocess" ]]; then
    if [[ "$st" -ne 0 ]]; then
      if cache_features_after_preprocess "$sample" "$mtype" "$type_out"; then
        local cache_file="$FEATURES_CACHE_ROOT/$sample/$mtype/features.tsv"
        log "FFPERASE preprocess returned non-zero, but features.tsv exists: sample=$sample type=$mtype status=$st"
        printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "DONE_FEATURES_RECOVERED" "$nvar" "$cache_file; recovered despite exit_status_${st}; see ${console}" > "$status_file"
        return 0
      fi
      log "FFPERASE preprocess failed: sample=$sample type=$mtype status=$st"
      printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "FAILED_PREPROCESS" "$nvar" "exit_status_${st}; see ${console}" > "$status_file"
      if [[ "$FAIL_ON_ERROR" == "true" ]]; then
        return "$st"
      fi
      return 0
    fi

    if cache_features_after_preprocess "$sample" "$mtype" "$type_out"; then
      local cache_file="$FEATURES_CACHE_ROOT/$sample/$mtype/features.tsv"
      printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "DONE_FEATURES_ONLY" "$nvar" "$cache_file" > "$status_file"
      log "Done preprocess-only: sample=$sample type=$mtype features=$cache_file"
      if [[ "$KEEP_WORK" != "true" ]]; then
        rm -rf "$workdir" "$launch_dir"
      fi
      return 0
    fi

    log "FFPERASE preprocess ended but features.tsv missing: sample=$sample type=$mtype"
    printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "FAILED_NO_FEATURES_TSV" "$nvar" "missing ${type_out}/preprocess/features.tsv" > "$status_file"
    if [[ "$FAIL_ON_ERROR" == "true" ]]; then
      return 1
    fi
    return 0
  fi

  # For full mode, keep a reusable copy of features.tsv whenever it exists.
  if [[ "$run_step" == "full" ]]; then
    cache_features_after_preprocess "$sample" "$mtype" "$type_out" >/dev/null 2>&1 || true
  fi

  if [[ "$st" -ne 0 ]]; then
    if [[ -s "$classified" ]]; then
      log "FFPERASE returned non-zero after classification, but classified TSV exists: sample=$sample type=$mtype status=$st"
      log "Treating as classified output recovered. Most likely PLOT_REPORT failed; see ${console}"
      cp -f "$classified" "$OUTROOT/classified/${sample}.classified_df_${mtype}.tsv"
      printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "DONE_CLASSIFIED_PLOT_FAILED" "$nvar" "$classified; wrapper recovered despite exit_status_${st}; see ${console}" > "$status_file"
      if [[ "$KEEP_WORK" != "true" ]]; then
        rm -rf "$workdir" "$launch_dir"
      fi
      return 0
    fi

    log "FFPERASE failed: sample=$sample type=$mtype status=$st"
    printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "FAILED" "$nvar" "exit_status_${st}; see ${console}" > "$status_file"
    if [[ "$FAIL_ON_ERROR" == "true" ]]; then
      return "$st"
    fi
    return 0
  fi

  if [[ ! -s "$classified" ]]; then
    log "FFPERASE ended but classified output missing: sample=$sample type=$mtype"
    printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "FAILED_NO_CLASSIFIED_TSV" "$nvar" "missing ${classified}" > "$status_file"
    if [[ "$FAIL_ON_ERROR" == "true" ]]; then
      return 1
    fi
    return 0
  fi

  cp -f "$classified" "$OUTROOT/classified/${sample}.classified_df_${mtype}.tsv"
  printf 'sample\ttype\tstatus\tvariants\tmessage\n%s\t%s\t%s\t%s\t%s\n' "$sample" "$mtype" "DONE" "$nvar" "$classified" > "$status_file"
  log "Done: sample=$sample type=$mtype classified=$classified"

  if [[ "$KEEP_WORK" != "true" ]]; then
    rm -rf "$workdir" "$launch_dir"
  fi

  return 0
}

run_row() {
  local sample="$1" final_vcf="$2" bam="$3" bai="$4" thresholds="$5"
  if [[ "$RUN_SNVS" == "true" ]]; then
    run_one "$sample" "$final_vcf" "$bam" "$bai" "$thresholds" "snvs"
  fi
  if [[ "$RUN_INDELS" == "true" ]]; then
    run_one "$sample" "$final_vcf" "$bam" "$bai" "$thresholds" "indels"
  fi
}

active_jobs=0
pids=()
labels=()
declare -A SEEN_MANIFEST_SAMPLES=()

wait_for_one() {
  local pid="${pids[0]}"
  local label="${labels[0]}"
  set +e
  wait "$pid"
  local st=$?
  set -e
  if [[ "$st" -ne 0 ]]; then
    log "Worker failed: $label status=$st"
    if [[ "$FAIL_ON_ERROR" == "true" ]]; then
      exit "$st"
    fi
  fi
  pids=("${pids[@]:1}")
  labels=("${labels[@]:1}")
  active_jobs=$((active_jobs - 1))
}

while IFS=$'\t' read -r sample final_vcf bam bai thresholds; do
  [[ "$sample" == "sample" ]] && continue
  [[ -n "$sample" ]] || continue
  if [[ -n "${SEEN_MANIFEST_SAMPLES[$sample]:-}" ]]; then
    log "SKIP duplicate manifest sample within this invocation: sample=$sample vcf=$final_vcf"
    continue
  fi
  SEEN_MANIFEST_SAMPLES[$sample]=1
  if (( JOBS <= 1 )); then
    run_row "$sample" "$final_vcf" "$bam" "$bai" "$thresholds"
  else
    run_row "$sample" "$final_vcf" "$bam" "$bai" "$thresholds" &
    pids+=("$!")
    labels+=("$sample")
    active_jobs=$((active_jobs + 1))
    while (( active_jobs >= JOBS )); do
      wait_for_one
    done
  fi
done < "$MANIFEST_OUT"

while (( active_jobs > 0 )); do
  wait_for_one
done

MERGED_LIST="$(python3 - "$OUTROOT" <<'PY_MERGE'
import csv
import glob
import os
import sys
from pathlib import Path
outroot = Path(sys.argv[1])
classified_dir = outroot / 'classified'
for mtype in ('snvs', 'indels'):
    files = sorted(classified_dir.glob(f'*.classified_df_{mtype}.tsv'))
    if not files:
        continue
    out = classified_dir / f'all_samples.ffperase_{mtype}.tsv'
    wrote_header = False
    with open(out, 'w', newline='') as oh:
        writer = None
        for f in files:
            sample = f.name.split('.classified_df_')[0]
            with open(f, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                if reader.fieldnames is None:
                    continue
                fieldnames = ['SAMPLE', 'MUTATION_TYPE'] + reader.fieldnames
                if not wrote_header:
                    writer = csv.DictWriter(oh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
                    writer.writeheader()
                    wrote_header = True
                for row in reader:
                    row2 = {'SAMPLE': sample, 'MUTATION_TYPE': mtype}
                    row2.update(row)
                    writer.writerow(row2)
    print(out)
PY_MERGE
)"
while IFS= read -r merged; do
  [[ -n "$merged" ]] && log "Merged table: $merged"
done <<< "$MERGED_LIST"

log "Done. Main output files:"
find "$OUTROOT" -path '*/classify/classified_df_*.tsv' -type f -print | sort | tee -a "$LOG"
find "$OUTROOT/classified" -maxdepth 1 -type f -name 'all_samples.ffperase_*.tsv' -print | sort | tee -a "$LOG"

failed_count="$(find "$OUTROOT/status" -type f -name '*.status.tsv' -exec awk 'NR==2 && $3 ~ /^FAILED/ {c++} END{print c+0}' {} + 2>/dev/null || echo 0)"
if [[ "${failed_count:-0}" != "0" ]]; then
  log "WARNING: failed sample/type count=$failed_count. See $OUTROOT/status and per-run console logs."
  if [[ "$FAIL_ON_ERROR" == "true" ]]; then
    exit 1
  fi
fi

exit 0

###__FFPERASE_EMBEDDED_NF_WRAPPER_END__###
