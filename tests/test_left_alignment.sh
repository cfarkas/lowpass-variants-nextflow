#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${TEST_DIR}/.." && pwd -P)"

fail() {
  printf 'not ok - %s\n' "$*" >&2
  exit 1
}

for executable in bcftools samtools bgzip; do
  command -v "$executable" >/dev/null 2>&1 || fail "required executable not found: ${executable}"
done

TEST_TMP="$(mktemp -d)"
trap 'rm -rf -- "$TEST_TMP"' EXIT

REFERENCE="${TEST_TMP}/tiny.fa"
INPUT_VCF="${TEST_TMP}/right_shifted.vcf"
OUTPUT_VCF="${TEST_TMP}/left_aligned.vcf.gz"
NORMALIZE_LOG="${TEST_TMP}/normalize.log"

printf '>chr1\nCAAAAAG\n' > "$REFERENCE"
samtools faidx "$REFERENCE"

printf '%s\n' \
  '##fileformat=VCFv4.2' \
  '##contig=<ID=chr1,length=7>' \
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
  $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE' \
  $'chr1\t5\t.\tAA\tA\t60\tPASS\t.\tGT\t0/1' \
  > "$INPUT_VCF"

# shellcheck source=../bin/nf_common.sh
source "${PACKAGE_DIR}/bin/nf_common.sh"
normalize_vcf "$INPUT_VCF" "$OUTPUT_VCF" "$REFERENCE" "$NORMALIZE_LOG" SAMPLE TEST_EMPTY

[[ -s "$OUTPUT_VCF" ]] || fail "normalization did not create a compressed VCF"
observed="$(bcftools query -f $'%CHROM\t%POS\t%REF\t%ALT\n' "$OUTPUT_VCF")"
expected=$'chr1\t1\tCA\tC'
[[ "$observed" == "$expected" ]] || fail \
  "right-shifted deletion was not left-aligned; expected '${expected}', observed '${observed}'"

record_count="$(bcftools index -n "$OUTPUT_VCF" 2>/dev/null)" || fail "VCF index is missing or invalid"
[[ "$record_count" == "1" ]] || fail "expected one indexed variant, found ${record_count}"

printf 'ok - bcftools norm left-aligns a homopolymer deletion and creates a valid index\n'
