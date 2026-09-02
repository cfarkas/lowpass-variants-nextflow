#!/usr/bin/env bash
set -Eeuo pipefail

DATA_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
DEST="${1:-}"

if [[ -z "$DEST" ]]; then
  printf 'Usage: %s DESTINATION\n' "$0" >&2
  exit 2
fi

for tool in samtools bcftools bgzip tabix; do
  command -v "$tool" >/dev/null 2>&1 || {
    printf 'Missing required command: %s\n' "$tool" >&2
    exit 127
  }
done

mkdir -p "$DEST"
DEST="$(cd -- "$DEST" && pwd -P)"
cp "$DATA_DIR/tiny.fa" "$DATA_DIR/SAMPLE.sam" "$DATA_DIR/known-sites.vcf" "$DEST/"

samtools faidx "$DEST/tiny.fa"
samtools dict -u tiny.fa -o "$DEST/tiny.dict" "$DEST/tiny.fa"
samtools view --no-PG -b -o "$DEST/SAMPLE.unsorted.bam" "$DEST/SAMPLE.sam"
samtools sort --no-PG -o "$DEST/SAMPLE.bam" "$DEST/SAMPLE.unsorted.bam"
samtools index "$DEST/SAMPLE.bam"
bgzip -c "$DEST/known-sites.vcf" > "$DEST/known-sites.vcf.gz"
tabix -f -p vcf "$DEST/known-sites.vcf.gz"

rm -f "$DEST/SAMPLE.unsorted.bam"
samtools quickcheck "$DEST/SAMPLE.bam"
bcftools index --stats "$DEST/known-sites.vcf.gz" >/dev/null
printf 'Generated and validated tiny data under %s\n' "$DEST"
