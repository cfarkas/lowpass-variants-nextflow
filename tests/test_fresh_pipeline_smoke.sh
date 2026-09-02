#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${TEST_DIR}/.." && pwd -P)"
PIPELINE="${PACKAGE_DIR}/main.nf"

fail() {
  printf 'not ok - %s\n' "$*" >&2
  exit 1
}

skip() {
  printf 'ok - SKIP fresh/BQSR pipeline smoke: %s\n' "$*"
  exit 0
}

if [[ -n "${LOWPASS_TEST_TOOL_ENV:-}" ]]; then
  [[ -d "${LOWPASS_TEST_TOOL_ENV}/bin" ]] ||
    fail "LOWPASS_TEST_TOOL_ENV has no bin directory: ${LOWPASS_TEST_TOOL_ENV}"
  export PATH="${LOWPASS_TEST_TOOL_ENV}/bin:${PATH}"
fi

for executable in gatk samtools bcftools bgzip tabix python; do
  command -v "$executable" >/dev/null 2>&1 || skip "required executable not found: ${executable}"
done
command -v timeout >/dev/null 2>&1 || skip "required executable not found: timeout"

if [[ -n "${NEXTFLOW_BIN:-}" ]]; then
  NEXTFLOW="${NEXTFLOW_BIN}"
else
  NEXTFLOW="$(command -v nextflow 2>/dev/null || true)"
fi
[[ -n "$NEXTFLOW" && -x "$NEXTFLOW" ]] || skip "Nextflow was not found"
[[ -f "$PIPELINE" ]] || fail "pipeline is missing: ${PIPELINE}"

TEST_TMP="$(mktemp -d)"
trap 'rm -rf -- "$TEST_TMP"' EXIT

REFERENCE="${TEST_TMP}/tiny.fa"
SAM="${TEST_TMP}/SAMPLE.sam"
UNSORTED_BAM="${TEST_TMP}/SAMPLE.unsorted.bam"
INPUT_BAM="${TEST_TMP}/SAMPLE.bam"
KNOWN_SITES="${TEST_TMP}/known-sites.vcf"
OUTDIR="${TEST_TMP}/results"
WORKDIR="${TEST_TMP}/work"
RUNDIR="${TEST_TMP}/run"
RUN_LOG="${TEST_TMP}/nextflow.log"

reference_sequence=""
for _ in {1..50}; do
  reference_sequence+="ACGT"
done
printf '>chr1\n%s\n' "$reference_sequence" > "$REFERENCE"
samtools faidx "$REFERENCE"
gatk CreateSequenceDictionary \
  -R "$REFERENCE" \
  -O "${REFERENCE%.fa}.dict" \
  > "${TEST_TMP}/create-dictionary.log" 2>&1

printf '%s\n' \
  '##fileformat=VCFv4.2' \
  '##contig=<ID=chr1,length=200>' \
  $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' \
  $'chr1\t25\t.\tA\tC\t.\tPASS\t.' \
  > "$KNOWN_SITES"
gatk IndexFeatureFile -I "$KNOWN_SITES" > "${TEST_TMP}/index-known-sites.log" 2>&1
[[ -s "${KNOWN_SITES}.idx" ]] || fail "GATK did not index the known-sites VCF"

read_one="${reference_sequence:0:50}"
read_two="${reference_sequence:50:50}"
quality="$(printf 'I%.0s' {1..50})"
printf '%s\n' \
  $'@HD\tVN:1.6\tSO:coordinate' \
  $'@SQ\tSN:chr1\tLN:200' \
  $'@RG\tID:rg1\tSM:SAMPLE\tLB:lib1\tPL:ILLUMINA\tPU:unit1' \
  "$(printf 'read1\t0\tchr1\t1\t60\t50M\t*\t0\t0\t%s\t%s\tRG:Z:rg1' "$read_one" "$quality")" \
  "$(printf 'read2\t0\tchr1\t51\t60\t50M\t*\t0\t0\t%s\t%s\tRG:Z:rg1' "$read_two" "$quality")" \
  > "$SAM"
samtools view -b -o "$UNSORTED_BAM" "$SAM"
samtools sort -o "$INPUT_BAM" "$UNSORTED_BAM"
samtools quickcheck "$INPUT_BAM"
[[ ! -e "${INPUT_BAM}.bai" && ! -e "${INPUT_BAM%.bam}.bai" ]] ||
  fail "fixture BAM unexpectedly had an input index"

mkdir -p "$RUNDIR" "$WORKDIR"
status=0
set +e
(
  cd "$RUNDIR"
  timeout "${FRESH_SMOKE_TIMEOUT:-300}" \
    env \
      NXF_SYNTAX_PARSER=v2 \
      NXF_OFFLINE=true \
      NXF_DISABLE_CHECK_LATEST=true \
      "$NEXTFLOW" run "$PIPELINE" \
        -ansi-log false \
        -work-dir "$WORKDIR" \
        --fresh true \
        --ffpe false \
        --input "$INPUT_BAM" \
        --outdir "$OUTDIR" \
        --ref "$REFERENCE" \
        --known_sites "$KNOWN_SITES" \
        --skip_mutect2 true \
        --skip_freebayes true \
        --skip_bcftools true \
        --auto_thresholds false \
        --threads 1 \
        --prepare_cpus 1 \
        --mutect2_cpus 1 \
        --freebayes_cpus 1 \
        --bcftools_cpus 1 \
        --normalize_cpus 1 \
        --finalize_cpus 1 \
        --max_samples_parallel 1 \
        --max_prepare_parallel 1 \
        --max_mutect2_parallel 1 \
        --max_freebayes_parallel 1 \
        --max_bcftools_parallel 1 \
        --max_normalize_parallel 1 \
        --max_finalize_parallel 1
) > "$RUN_LOG" 2>&1
status=$?
set -e

if [[ "$status" -ne 0 ]]; then
  sed -n '1,300p' "$RUN_LOG" >&2
  fail "fresh/BQSR Nextflow smoke failed with status ${status}"
fi

RECAL_TABLE="${OUTDIR}/preprocessed_bam/SAMPLE.recal_data.table"
FINAL_VCF="${OUTDIR}/final_vcf/SAMPLE.final_variants.vcf.gz"
FINAL_TBI="${FINAL_VCF}.tbi"
THRESHOLDS="${OUTDIR}/coverage/SAMPLE.adaptive_thresholds.tsv"

[[ -s "$RECAL_TABLE" ]] || fail "missing BQSR recalibration table"
[[ -s "$FINAL_VCF" ]] || fail "missing final normalized VCF"
[[ -s "$FINAL_TBI" ]] || fail "missing final VCF tabix index"
[[ -s "$THRESHOLDS" ]] || fail "missing adaptive-threshold table"
bcftools index --stats "$FINAL_VCF" >/dev/null 2>&1 || fail "final VCF index is invalid"

if bcftools view -h "$FINAL_VCF" | grep -q 'FFPE_'; then
  fail "fresh final VCF contains FFPE-specific INFO declarations"
fi
if head -n 1 "$THRESHOLDS" | grep -qi 'ffpe'; then
  fail "fresh threshold output contains FFPE-specific columns"
fi
if find "$OUTDIR" -maxdepth 3 -iname '*ffperase*' -print -quit | grep -q .; then
  fail "fresh run unexpectedly created FFPErase output"
fi
[[ ! -e "${INPUT_BAM}.bai" && ! -e "${INPUT_BAM%.bam}.bai" ]] ||
  fail "workflow mutated the input area by creating a BAM index"

printf 'ok - real fresh pipeline runs BQSR, skips callers, emits indexed normalized VCF, and bypasses FFPErase\n'
