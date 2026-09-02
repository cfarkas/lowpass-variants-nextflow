#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${TEST_DIR}/.." && pwd -P)"
DATA_DIR="${PACKAGE_DIR}/examples/tiny"

fail() {
  printf 'not ok - %s\n' "$*" >&2
  exit 1
}

skip() {
  printf 'ok - SKIP bundled example data: %s\n' "$*"
  exit 0
}

for tool in samtools bcftools; do
  command -v "$tool" >/dev/null 2>&1 || skip "required executable not found: $tool"
done

for path in \
  tiny.fa tiny.fa.fai tiny.dict SAMPLE.sam SAMPLE.bam SAMPLE.bam.bai \
  known-sites.vcf known-sites.vcf.gz known-sites.vcf.gz.tbi
do
  [[ -s "${DATA_DIR}/${path}" ]] || fail "missing bundled fixture: ${path}"
done

reference_length="$(awk '!/^>/ {gsub(/[[:space:]]/, ""); n += length} END {print n+0}' "${DATA_DIR}/tiny.fa")"
[[ "$reference_length" == 200 ]] || fail "tiny reference length is ${reference_length}, expected 200"
grep -Fqx $'chr1\t200\t6\t200\t201' "${DATA_DIR}/tiny.fa.fai" ||
  fail "unexpected FASTA index"
grep -Fq $'@SQ\tSN:chr1\tLN:200' "${DATA_DIR}/tiny.dict" ||
  fail "reference dictionary lacks chr1 length 200"

samtools quickcheck "${DATA_DIR}/SAMPLE.bam" || fail "bundled BAM failed quickcheck"
samtools idxstats "${DATA_DIR}/SAMPLE.bam" |
  grep -Fqx $'chr1\t200\t2\t0' || fail "unexpected BAM index statistics"
samtools view -H --no-PG "${DATA_DIR}/SAMPLE.bam" |
  grep -Fq $'@RG\tID:rg1\tSM:SAMPLE' || fail "BAM read group/sample is missing"

bcftools index --stats "${DATA_DIR}/known-sites.vcf.gz" |
  grep -Fqx $'chr1\t200\t1' || fail "unexpected known-sites index statistics"
known_site="$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${DATA_DIR}/known-sites.vcf.gz")"
[[ "$known_site" == $'chr1\t25\tA\tC' ]] || fail "unexpected known-sites record"

TEST_TMP="$(mktemp -d)"
trap 'rm -rf -- "$TEST_TMP"' EXIT
FAKE_LOG="${TEST_TMP}/nextflow.calls"
touch "$FAKE_LOG"

env \
  NEXTFLOW_BIN="${TEST_DIR}/fixtures/nextflow" \
  FAKE_NF_LOG="$FAKE_LOG" \
  FAKE_NF_MODE=both \
  NEXTFLOW_PROFILE= \
  OUTDIR="${TEST_TMP}/fresh-results" \
  bash "${PACKAGE_DIR}/examples/run_minimal_fresh.sh" >/dev/null
grep -Fq -- '--fresh' "$FAKE_LOG" || fail "fresh example did not forward --fresh"
grep -Fq -- "${DATA_DIR}/SAMPLE.bam" "$FAKE_LOG" || fail "fresh example did not use bundled BAM"

: > "$FAKE_LOG"
env \
  NEXTFLOW_BIN="${TEST_DIR}/fixtures/nextflow" \
  FAKE_NF_LOG="$FAKE_LOG" \
  FAKE_NF_MODE=both \
  NEXTFLOW_PROFILE= \
  OUTDIR="${TEST_TMP}/ffpe-results" \
  bash "${PACKAGE_DIR}/examples/run_minimal_ffpe.sh" >/dev/null
grep -Fq -- '--ffpe' "$FAKE_LOG" || fail "FFPE example did not forward --ffpe"
grep -Fq -- '--skip_mutect2' "$FAKE_LOG" || fail "FFPE smoke does not skip callers"

printf 'ok - bundled synthetic data and minimal fresh/FFPE launchers are valid\n'
