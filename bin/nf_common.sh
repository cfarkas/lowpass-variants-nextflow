#!/usr/bin/env bash
set -Eeuo pipefail

sanitize_id() {
  local s="$1"
  s="${s// /_}"
  printf '%s' "$s" | sed 's/[^A-Za-z0-9_.-]/_/g'
}

link_or_copy() {
  local src="$1"
  local dest="$2"
  rm -f "$dest"
  if ln "$src" "$dest" 2>/dev/null; then
    return 0
  fi
  if ln -s "$src" "$dest" 2>/dev/null; then
    return 0
  fi
  cp -f "$src" "$dest"
}

bam_sm() {
  local bam="$1"
  samtools view -H "$bam" \
    | awk -F'\t' '$1=="@RG"{for(i=1;i<=NF;i++){if($i ~ /^SM:/){sub(/^SM:/,"",$i); print $i}}}' \
    | sort -u \
    | head -1
}

bam_is_coordinate_sorted() {
  local bam="$1"
  samtools view -H "$bam" \
    | awk '$1=="@HD"{for(i=1;i<=NF;i++){if($i=="SO:coordinate") found=1}} END{exit(found?0:1)}'
}

write_empty_vcf() {
  local sample="$1"
  local out_vcf_gz="$2"
  local source_label="${3:-EMPTY}"
  mkdir -p "$(dirname "$out_vcf_gz")"
  {
    echo '##fileformat=VCFv4.2'
    echo "##source=${source_label}"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$sample"
  } | bgzip -c > "$out_vcf_gz"
  tabix -f -p vcf "$out_vcf_gz"
  [[ -s "${out_vcf_gz}.tbi" ]]
}

index_vcf() {
  local vcf="$1"
  [[ -s "$vcf" ]] || return 1
  rm -f "${vcf}.csi"
  if ! bcftools index -f -t "$vcf" 2>/dev/null \
    && ! tabix -f -p vcf "$vcf" 2>/dev/null; then
    echo "ERROR: could not create the required TBI index for $vcf" >&2
    echo "ERROR: this workflow currently requires TBI-compatible contig lengths/positions" >&2
    return 1
  fi
  [[ -s "${vcf}.tbi" ]]
}

normalize_vcf() {
  local input_vcf="$1"
  local output_vcf="$2"
  local ref="$3"
  local log="$4"
  local sample="${5:-SAMPLE}"
  local source_label="${6:-EMPTY_NORMALIZE_INPUT}"

  mkdir -p "$(dirname "$output_vcf")" "$(dirname "$log")"

  if [[ ! -s "$input_vcf" ]]; then
    write_empty_vcf "$sample" "$output_vcf" "$source_label"
    return 0
  fi

  if ! bcftools norm \
    -m -any \
    -f "$ref" \
    -Oz \
    -o "$output_vcf" \
    "$input_vcf" \
    > "$log" 2>&1; then
    echo "ERROR: bcftools norm failed for $input_vcf; refusing to emit an unnormalized VCF" >&2
    cat "$log" >&2 || true
    return 1
  fi
  [[ -s "$output_vcf" ]] || {
    echo "ERROR: bcftools norm produced no output for $input_vcf" >&2
    return 1
  }

  index_vcf "$output_vcf"
}

ensure_bam_index() {
  local bam="$1"
  if [[ -s "${bam}.bai" ]]; then
    return 0
  fi
  if [[ -s "${bam%.bam}.bai" ]]; then
    return 0
  fi
  if [[ -s "${bam}.csi" ]]; then
    return 0
  fi
  samtools index -@ "${THREADS:-1}" "$bam"
}
