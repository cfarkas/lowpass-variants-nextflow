#!/usr/bin/env bash
set -Eeuo pipefail

TEST_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PACKAGE_DIR="$(cd -- "${TEST_DIR}/.." && pwd -P)"
EXAMPLE="${PACKAGE_DIR}/examples/run_minimal_fresh.sh"

fail() {
  printf 'not ok - %s\n' "$*" >&2
  exit 1
}

skip() {
  printf 'ok - SKIP bundled fresh quickstart: %s\n' "$*"
  exit 0
}

if [[ -n "${LOWPASS_TEST_TOOL_ENV:-}" ]]; then
  [[ -d "${LOWPASS_TEST_TOOL_ENV}/bin" ]] ||
    fail "LOWPASS_TEST_TOOL_ENV has no bin directory: ${LOWPASS_TEST_TOOL_ENV}"
  export PATH="${LOWPASS_TEST_TOOL_ENV}/bin:${PATH}"
fi

for executable in gatk samtools bcftools bgzip tabix freebayes python timeout; do
  command -v "$executable" >/dev/null 2>&1 ||
    skip "required executable not found: ${executable}"
done

if [[ -n "${NEXTFLOW_BIN:-}" ]]; then
  [[ -x "$NEXTFLOW_BIN" ]] || skip "NEXTFLOW_BIN is not executable: ${NEXTFLOW_BIN}"
else
  command -v nextflow >/dev/null 2>&1 || skip "Nextflow was not found"
fi

[[ -x "$EXAMPLE" ]] || fail "fresh quickstart is missing or not executable: ${EXAMPLE}"

TEST_TMP="$(mktemp -d)"
trap 'rm -rf -- "$TEST_TMP"' EXIT
OUTDIR="${TEST_TMP}/results"
RUN_LOG="${TEST_TMP}/fresh-quickstart.log"

status=0
set +e
(
  cd "$TEST_TMP"
  timeout "${FRESH_QUICKSTART_TIMEOUT:-600}" \
    env \
      NEXTFLOW_PROFILE= \
      OUTDIR="$OUTDIR" \
      NXF_SYNTAX_PARSER=v2 \
      NXF_OFFLINE=true \
      NXF_DISABLE_CHECK_LATEST=true \
      bash "$EXAMPLE"
) > "$RUN_LOG" 2>&1
status=$?
set -e

if [[ "$status" -ne 0 ]]; then
  sed -n '1,360p' "$RUN_LOG" >&2
  fail "bundled fresh quickstart failed with status ${status}"
fi

FINAL_VCF="${OUTDIR}/final_vcf/SAMPLE.final_variants.vcf.gz"
RETAINED_TSV="${OUTDIR}/reports/per_sample/SAMPLE.final_variants.retained.tsv"
CALLER_GATE="${OUTDIR}/logs/SAMPLE.callers.ready.tsv"
RECAL_TABLE="${OUTDIR}/preprocessed_bam/SAMPLE.recal_data.table"
FREEBAYES_VCF="${OUTDIR}/freebayes/SAMPLE.freebayes.vcf.gz"
BCFTOOLS_VCF="${OUTDIR}/bcftools/SAMPLE.bcftools.vcf.gz"

for required in \
  "$FINAL_VCF" "${FINAL_VCF}.tbi" "$RETAINED_TSV" "$CALLER_GATE" \
  "$RECAL_TABLE" "$FREEBAYES_VCF" "${FREEBAYES_VCF}.tbi" \
  "$BCFTOOLS_VCF" "${BCFTOOLS_VCF}.tbi"
do
  [[ -s "$required" ]] || fail "missing quickstart output: ${required}"
done

bcftools index --stats "$FINAL_VCF" >/dev/null 2>&1 ||
  fail "final quickstart VCF index is invalid"

mapfile -t final_records < <(
  bcftools query \
    -f $'%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DECISION\t%INFO/VOTE_COUNT\t%INFO/SOURCE\t%INFO/PIPELINE_MODE\n' \
    "$FINAL_VCF"
)
[[ "${#final_records[@]}" -eq 1 ]] ||
  fail "expected exactly one final VCF record, found ${#final_records[@]}"
expected_record=$'chr1\t25\tA\tC\tPASS\tHIGH_CONFIDENCE\t2\tFREEBAYES,BCFTOOLS\tfresh'
[[ "${final_records[0]}" == "$expected_record" ]] ||
  fail "unexpected final VCF record: ${final_records[0]}"

awk -F '\t' '
  NR == 2 && $1 == "SAMPLE" && $3 == "skipped" &&
  $5 == "completed" && $7 == "completed" && $8 == "caller_stages_finished" {
    ok = 1
  }
  END { exit(ok ? 0 : 1) }
' "$CALLER_GATE" ||
  fail "caller gate does not report skipped Mutect2 and completed FreeBayes/BCFtools"

awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) column[$i] = i
    next
  }
  {
    rows++
    if ($(column["chrom"]) == "chr1" && $(column["pos"]) == "25" &&
        $(column["ref"]) == "A" && $(column["alt"]) == "C" &&
        $(column["decision"]) == "HIGH_CONFIDENCE" &&
        $(column["vote_count"]) == "2") {
      matched++
    }
  }
  END { exit(rows == 1 && matched == 1 ? 0 : 1) }
' "$RETAINED_TSV" ||
  fail "retained report is not exactly the expected one-variant high-confidence row"

if bcftools view -h "$FINAL_VCF" | grep -qi 'FFPE_'; then
  fail "fresh quickstart VCF contains FFPE-specific INFO declarations"
fi
if head -n 1 "$RETAINED_TSV" | grep -Eqi 'ffpe|is_CT_or_GA'; then
  fail "fresh retained report contains FFPE-specific columns"
fi
if find "$OUTDIR" -iname '*ffpe*' -print -quit | grep -q .; then
  fail "fresh quickstart unexpectedly created FFPE output"
fi

printf 'ok - bundled fresh quickstart emits one chr1:25 A>C high-confidence two-caller record\n'
