#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Unified low-pass final per-sample VCF workflow
 * ---------------------------------------------------------------
 * Visible pipeline stages:
 *   RESOLVE_BAMS -> PREPARE_REFERENCE + PREPARE_SCOPE -> PREPARE_BAM (BQSR)
 *   -> COVERAGE -> CALIBRATE_THRESHOLDS
 *   -> MUTECT2_CALL -> MUTECT2_ORIENTATION_MODEL -> MUTECT2_FILTER -> NORMALIZE_MUTECT2
 *   -> FREEBAYES_CALL -> NORMALIZE_FREEBAYES
 *   -> BCFTOOLS_CALL -> NORMALIZE_BCFTOOLS
 *   -> MERGE_CALLERS_PER_SAMPLE -> NORMALIZE_FINAL_VCF
 *   -> FFPERASE (only with --ffpe)
 *
 * Run behavior:
 *   - Exactly one sample mode is required: --ffpe or --fresh.
 *   - Genome-wide by default; optional --genes or --bed/--target_bed.
 *   - One independent final VCF per sample only.
 *   - No cohort-level union VCF and no all_samples.final_variants.* TSV.
 *   - Mutect2, FreeBayes, and BCFtools evidence are kept as SOURCE annotations.
 *   - FFPE C>T/G>A risk and FFPErase are enabled only in --ffpe mode.
 *   - BQSR is enabled by default and requires --known_sites; --skip_bqsr is explicit.
 *   - Every caller VCF and the final VCF are split, normalized, and
 *     left-aligned against --ref with bcftools norm.
 */



def truthyParam(v) {
    if( v == null ) return false
    if( v instanceof Boolean ) return v
    return ['true', '1', 'yes', 'y'].contains(v.toString().toLowerCase())
}

def asIntParam(v, fallback) {
    if( v == null ) return fallback as int
    def s = v.toString().trim()
    if( s == '' || s == 'null' ) return fallback as int
    return s as int
}

def threadsParam() {
    return asIntParam(params.threads, 12)
}

def prepareCpus() {
    return asIntParam(params.prepare_cpus, 12)
}

def mutect2Cpus() {
    return asIntParam(params.mutect2_cpus, 12)
}

def freebayesCpus() {
    return asIntParam(params.freebayes_cpus, 1)
}

def bcftoolsCpus() {
    return asIntParam(params.bcftools_cpus, 8)
}

def normalizeCpus() {
    return asIntParam(params.normalize_cpus, 4)
}

def finalizeCpus() {
    return asIntParam(params.finalize_cpus, 4)
}

def baseParallel() {
    return asIntParam(params.max_samples_parallel ?: params.max_parallel_samples, 4)
}

def prepareParallel() {
    return asIntParam(params.max_prepare_parallel, baseParallel())
}

def mutect2Parallel() {
    return asIntParam(params.max_mutect2_parallel, baseParallel())
}

def freebayesParallel() {
    return asIntParam(params.max_freebayes_parallel, baseParallel())
}

def bcftoolsParallel() {
    return asIntParam(params.max_bcftools_parallel, baseParallel())
}

def normalizeParallel() {
    return asIntParam(params.max_normalize_parallel, 4)
}

def finalizeParallel() {
    return asIntParam(params.max_finalize_parallel, baseParallel())
}

def workflowHelp(fullHelp=false, showHidden=false, assetsDir='assets') {
    def helpName = (fullHelp || showHidden) ? 'help_full.txt' : 'help.txt'
    def helpFile = new java.io.File(assetsDir, helpName)
    if( helpFile.exists() ) {
        return helpFile.text
    }
    return "ERROR: help file not found: ${helpFile}\n"
}

process RESOLVE_BAMS {
    tag 'resolve_bams'
    cache false
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 1
    memory '2 GB'
    time '2h'
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: 'selected_samples.tsv', overwrite: params.overwrite

    output:
    path 'selected_samples.tsv', emit: samplesheet

    script:
    def recursiveFlag = params.recursive ? '--recursive' : ''
    def metadataPath = params.database_tsv ?: params.metadata_tsv ?: ''
    def metadataArg = metadataPath ? "--database-tsv \"${metadataPath}\"" : ''
    """
    set -Eeuo pipefail
    mkdir -p "${params.outdir}/reports"
    python "${projectDir}/bin/resolve_bams.py" \
      --input "${params.input}" \
      --output selected_samples.tsv \
      --samples "${params.samples}" \
      --cohort "${params.cohort}" \
      ${recursiveFlag} \
      ${metadataArg}

    echo "Selected samples:"
    cat selected_samples.tsv
    """
}

process PREPARE_REFERENCE {
    tag 'reference_sidecars'
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 2
    memory '8 GB'
    time '4h'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'reference.ready.txt', overwrite: params.overwrite

    input:
    path source_reference

    output:
    tuple path('reference.fa'), path('reference.fa.fai'), path('reference.dict'), path('reference.ready.txt'), emit: bundle

    script:
    """
    set -Eeuo pipefail
    [[ -s "${source_reference}" ]] || { echo "ERROR: missing reference FASTA: ${source_reference}" >&2; exit 1; }
    cp -L "${source_reference}" reference.fa
    samtools faidx reference.fa
    gatk CreateSequenceDictionary -R reference.fa -O reference.dict

    {
      echo "reference=reference.fa"
      echo "fai=reference.fa.fai"
      echo "dict=reference.dict"
      date
    } > reference.ready.txt
    """
}

process PREPARE_SCOPE {
    tag 'prepare_scope'
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 2
    memory '8 GB'
    time '4h'
    publishDir "${params.outdir}/targets", mode: 'copy', pattern: 'targets.*.bed', overwrite: params.overwrite
    publishDir "${params.outdir}/targets", mode: 'copy', pattern: 'scope.env', overwrite: params.overwrite
    publishDir "${params.outdir}/targets", mode: 'copy', pattern: 'gencode*.gtf.gz', overwrite: params.overwrite

    input:
    tuple val(bed_name), val(gtf_name), path(scope_resource_files)

    output:
    tuple val('scope'), val('unused'), path('scope.env'), path('targets.annotated.bed'), path('targets.merged.bed'), emit: scope_tuple
    path 'gencode*.gtf.gz', optional: true, emit: downloaded_gtf

    script:
    def bedPath = bed_name ?: ''
    def gtfPath = gtf_name ?: ''
    """
    set -Eeuo pipefail
    python "${projectDir}/bin/prepare_scope.py" \
      --genes "${params.genes}" \
      --bed "${bedPath}" \
      --gtf "${gtfPath}" \
      --gencode-version "${params.gencode_version}" \
      --promoter-upstream "${params.promoter_upstream}" \
      --promoter-downstream "${params.promoter_downstream}" \
      --outdir .

    cat scope.env
    """
}

process PREPARE_BAM {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.prepare_cpus
    memory params.memory
    time params.quick_time
    maxForks params.max_prepare_parallel
    publishDir "${params.outdir}/preprocessed_bam", mode: 'copy', pattern: '*.preprocessed.bam*', overwrite: params.overwrite
    publishDir "${params.outdir}/preprocessed_bam", mode: 'copy', pattern: '*.tumor_sm.txt', overwrite: params.overwrite
    publishDir "${params.outdir}/preprocessed_bam", mode: 'copy', pattern: '*.recal_data.table', overwrite: params.overwrite
    publishDir "${params.outdir}/preprocessed_bam", mode: 'copy', pattern: '*.markduplicates.metrics.txt', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.prepare_bam.*', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.bqsr.*.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(bam_path), val(scope_key), val(scope_unused), path(scope_env), path(target_annot), path(target_merged), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), val(known_sites_names), path(known_site_files)

    output:
    tuple val(sample), path("${sample}.preprocessed.bam"), path("${sample}.preprocessed.bam.bai"), path("${sample}.tumor_sm.txt"), emit: prepared
    path "${sample}.recal_data.table", optional: true, emit: recalibration
    path "${sample}.prepare_bam.*", optional: true, emit: logs
    path "${sample}.markduplicates.metrics.txt", optional: true, emit: duplicate_metrics
    path "${sample}.bqsr.*.log", optional: true, emit: bqsr_logs

    script:
    def forceRg = truthyParam(params.force_rg_sm_basename) ? 'true' : 'false'
    def skipBqsr = truthyParam(params.skip_bqsr) ? 'true' : 'false'
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"

    sample="${sample}"
    INPUT_BAM="${bam_path}"
    FINAL_BAM="${sample}.preprocessed.bam"
    WORK_DIR="prepare_${sample}"
    THREADS=${task.cpus}
    SKIP_BQSR="${skipBqsr}"
    KNOWN_SITES_RAW="${known_sites_names}"
    REF="${reference_fasta}"
    mkdir -p "\$WORK_DIR"

    echo "[PREPARE_BAM] sample=\$sample input=\$INPUT_BAM cpus=${task.cpus} duplicate_mode=${params.duplicate_mode} skip_bqsr=\$SKIP_BQSR"

    current="\$INPUT_BAM"
    need_sort=false

    if ! samtools quickcheck -v "\$INPUT_BAM" > "${sample}.prepare_bam.quickcheck.log" 2>&1; then
      echo "WARNING: input BAM failed quickcheck. Trying coordinate sort repair." >&2
      need_sort=true
    fi
    if ! bam_is_coordinate_sorted "\$INPUT_BAM"; then
      need_sort=true
    fi

    if [[ "\$need_sort" == "true" ]]; then
      current="\$WORK_DIR/${sample}.sorted.input.bam"
      samtools sort -@ "\$THREADS" -T "\$WORK_DIR/${sample}.sort" -o "\$current" "\$INPUT_BAM"
      samtools index -@ "\$THREADS" "\$current"
    fi

    sm="\$(bam_sm "\$current" || true)"
    if [[ "${forceRg}" == "true" || -z "\$sm" ]]; then
      rg_bam="\$WORK_DIR/${sample}.rg.bam"
      gatk AddOrReplaceReadGroups \
        -I "\$current" \
        -O "\$rg_bam" \
        -RGID "${sample}" \
        -RGLB "${sample}" \
        -RGPL ILLUMINA \
        -RGPU "${sample}" \
        -RGSM "${sample}" \
        --TMP_DIR "\$WORK_DIR" \
        --CREATE_INDEX false
      current="\$WORK_DIR/${sample}.rg.sorted.bam"
      samtools sort -@ "\$THREADS" -T "\$WORK_DIR/${sample}.rg.sort" -o "\$current" "\$rg_bam"
      samtools index -@ "\$THREADS" "\$current"
      sm="${sample}"
    fi

    [[ -n "\$sm" ]] || sm="${sample}"
    echo "\$sm" > "${sample}.tumor_sm.txt"

    if [[ "${params.duplicate_mode}" != "skip" ]]; then
      remove_dups=false
      if [[ "${params.duplicate_mode}" == "remove" ]]; then
        remove_dups=true
      fi
      duplicate_bam="\$WORK_DIR/${sample}.duplicates.bam"
      gatk MarkDuplicates \
        -I "\$current" \
        -O "\$duplicate_bam" \
        -M "${sample}.markduplicates.metrics.txt" \
        --REMOVE_DUPLICATES "\$remove_dups" \
        --CREATE_INDEX false \
        --TMP_DIR "\$WORK_DIR"
      current="\$duplicate_bam"
    fi

    # Never create an index beside a user-owned input BAM. When the original
    # BAM has no usable sidecar, stage a private copy in the task work area
    # before indexing it.
    if [[ "\$current" == "\$INPUT_BAM" \
          && ! -s "\$current.bai" \
          && ! -s "\${current%.bam}.bai" \
          && ! -s "\$current.csi" ]]; then
      staged_bam="\$WORK_DIR/${sample}.indexable.bam"
      cp -f "\$current" "\$staged_bam"
      current="\$staged_bam"
    fi
    ensure_bam_index "\$current"

    if [[ "\$SKIP_BQSR" == "true" ]]; then
      echo "[PREPARE_BAM] BQSR explicitly skipped for sample=\$sample" | tee "${sample}.bqsr.skipped.log"
      # The work directory is removed below, so the final BAM must not be a
      # symlink into it.
      cp -f "\$current" "\$FINAL_BAM"
      samtools index -@ "\$THREADS" "\$FINAL_BAM"
    else
      [[ -n "\$KNOWN_SITES_RAW" ]] || {
        echo "ERROR: BQSR requires --known_sites unless --skip_bqsr is used" >&2
        exit 1
      }

      IFS=',' read -r -a known_sites_array <<< "\$KNOWN_SITES_RAW"
      bqsr_known_args=()
      for known_site in "\${known_sites_array[@]}"; do
        known_site="\${known_site#"\${known_site%%[![:space:]]*}"}"
        known_site="\${known_site%"\${known_site##*[![:space:]]}"}"
        [[ -n "\$known_site" ]] || continue
        [[ -s "\$known_site" ]] || {
          echo "ERROR: known-sites VCF not found or empty: \$known_site" >&2
          exit 1
        }
        if [[ "\$known_site" == *.gz ]]; then
          [[ -s "\$known_site.tbi" || -s "\$known_site.csi" || -s "\$known_site.idx" ]] || {
            echo "ERROR: indexed known-sites VCF required (.tbi/.csi/.idx): \$known_site" >&2
            exit 1
          }
        else
          [[ -s "\$known_site.idx" ]] || {
            echo "ERROR: indexed known-sites VCF required (.idx): \$known_site" >&2
            exit 1
          }
        fi
        bqsr_known_args+=(--known-sites "\$known_site")
      done
      [[ "\${#bqsr_known_args[@]}" -gt 0 ]] || {
        echo "ERROR: --known_sites did not contain a usable VCF" >&2
        exit 1
      }

      gatk BaseRecalibrator \
        -R "\$REF" \
        -I "\$current" \
        "\${bqsr_known_args[@]}" \
        -O "${sample}.recal_data.table" \
        > "${sample}.bqsr.base_recalibrator.log" 2>&1

      gatk ApplyBQSR \
        -R "\$REF" \
        -I "\$current" \
        --bqsr-recal-file "${sample}.recal_data.table" \
        -O "\$FINAL_BAM" \
        --create-output-bam-index true \
        > "${sample}.bqsr.apply.log" 2>&1

      if [[ -s "\${FINAL_BAM%.bam}.bai" && ! -s "\$FINAL_BAM.bai" ]]; then
        mv "\${FINAL_BAM%.bam}.bai" "\$FINAL_BAM.bai"
      fi
      [[ -s "\$FINAL_BAM.bai" ]] || samtools index -@ "\$THREADS" "\$FINAL_BAM"
    fi

    samtools quickcheck "\$FINAL_BAM"
    [[ -s "\$FINAL_BAM.bai" ]] || {
      echo "ERROR: missing BAM index after preparation: \$FINAL_BAM.bai" >&2
      exit 1
    }
    rm -rf "\$WORK_DIR"
    """
}

process COVERAGE {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 1
    memory '4 GB'
    time params.quick_time
    maxForks params.max_prepare_parallel
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: '*', overwrite: params.overwrite

    input:
    tuple val(sample), path(bam), path(bai), path(sm_file), val(scope_key), val(scope_unused), path(scope_env), path(target_annot), path(target_merged)

    output:
    tuple val(sample), path("${sample}.genome_coverage.tsv"), path("${sample}.target_region_coverage.tsv"), path("${sample}.coverage_mode.txt"), emit: coverage_bundle

    script:
    """
    set -Eeuo pipefail
    source "${scope_env}"
    echo "[COVERAGE] sample=${sample} target_mode=\$TARGET_MODE"
    echo "\$TARGET_MODE" > "${sample}.coverage_mode.txt"

    # samtools coverage gives mean depth per contig and is suitable for genome-wide low-pass heuristics.
    samtools coverage "${bam}" > "${sample}.genome_coverage.tsv"

    if [[ "\$TARGET_MODE" != "genome_wide" && -s "${target_annot}" ]]; then
      bedtools coverage -a "${target_annot}" -b "${bam}" -mean > "${sample}.target_region_coverage.tsv"
    else
      : > "${sample}.target_region_coverage.tsv"
    fi
    """
}

process CALIBRATE_THRESHOLDS {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 1
    memory '2 GB'
    time '2h'
    maxForks params.max_prepare_parallel
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: '*.adaptive_thresholds.tsv', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.adaptive_thresholds.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(genome_coverage), path(target_coverage), path(mode_file)

    output:
    tuple val(sample), path("${sample}.adaptive_thresholds.tsv"), emit: thresholds
    path "${sample}.adaptive_thresholds.log", optional: true, emit: logs

    script:
    def sampleMode = truthyParam(params.ffpe) ? 'ffpe' : 'fresh'
    """
    set -Eeuo pipefail
    TARGET_MODE="\$(cat "${mode_file}" 2>/dev/null || echo genome_wide)"
    echo "[CALIBRATE_THRESHOLDS] sample=${sample} sample_mode=${sampleMode} target_mode=\$TARGET_MODE auto_thresholds=${params.auto_thresholds}" | tee "${sample}.adaptive_thresholds.log"
    CALIBRATE_ARGS=(
      python "${projectDir}/bin/calibrate_thresholds.py"
      --sample "${sample}"
      --sample-mode "${sampleMode}"
      --target-mode "\$TARGET_MODE"
      --genome-coverage "${genome_coverage}"
      --target-coverage "${target_coverage}"
      --auto-thresholds "${params.auto_thresholds}"
      --manual-min-dp "${params.min_dp}"
      --manual-min-alt-reads "${params.min_alt_reads}"
      --manual-min-af "${params.min_af}"
      --manual-vote-threshold "${params.vote_threshold}"
      --output "${sample}.adaptive_thresholds.tsv"
    )
    if [[ "${sampleMode}" == "ffpe" ]]; then
      CALIBRATE_ARGS+=(
        --manual-ffpe-ct-vaf-keep "${params.ffpe_ct_vaf_keep}"
        --manual-ffpe-ct-min-alt "${params.ffpe_ct_min_alt}"
      )
    fi
    "\${CALIBRATE_ARGS[@]}" | tee -a "${sample}.adaptive_thresholds.log"
    """
}

process MUTECT2_CALL {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.mutect2_cpus
    memory params.mutect2_memory
    time params.caller_time
    maxForks params.max_mutect2_parallel
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*.unfiltered.vcf.gz*', overwrite: params.overwrite
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*.f1r2.tar.gz', overwrite: params.overwrite
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*.mutect2_call_status.txt', overwrite: params.overwrite
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*.mutect2.*.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(bam), path(bai), path(sm_file), val(scope_key), val(scope_unused), path(scope_env), path(target_annot), path(target_merged), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), val(germline_name), val(pon_name), path(mutect_resource_files)

    output:
    tuple val(sample), path("${sample}.unfiltered.vcf.gz"), path("${sample}.unfiltered.vcf.gz.tbi"), path("${sample}.f1r2.tar.gz"), path("${sample}.unfiltered.vcf.gz.stats"), path("${sample}.mutect2_call_status.txt"), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), emit: raw
    path "${sample}.mutect2.*.log", optional: true, emit: logs

    script:
    def skip = truthyParam(params.skip_mutect2) ? 'true' : 'false'
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"
    source "${scope_env}"

    sample="${sample}"
    TUMOR_SM="\$(cat "${sm_file}")"
    raw="${sample}.unfiltered.vcf.gz"
    f1r2="${sample}.f1r2.tar.gz"
    stats="${sample}.unfiltered.vcf.gz.stats"
    call_status="${sample}.mutect2_call_status.txt"

    echo "[MUTECT2_CALL] sample=\$sample cpus=${task.cpus} tumor_sm=\$TUMOR_SM target_mode=\$TARGET_MODE"

    if [[ "${skip}" == "true" ]]; then
      write_empty_vcf "\$sample" "\$raw" "MUTECT2_SKIPPED"
      : > "\$f1r2"
      : > "\$stats"
      echo "skipped" > "\$call_status"
      exit 0
    fi

    MUTECT_ARGS=(
      gatk --java-options "-Xmx12g" Mutect2
      -R "${reference_fasta}"
      -I "${bam}"
      -tumor "\$TUMOR_SM"
      --f1r2-tar-gz "\$f1r2"
      --native-pair-hmm-threads ${task.cpus}
      --max-reads-per-alignment-start 50
      --annotation AlleleFraction
      --annotation AS_FisherStrand
      --annotation TandemRepeat
      -O "\$raw"
    )
    if [[ "\$TARGET_MODE" != "genome_wide" && -s "${target_merged}" ]]; then
      MUTECT_ARGS+=(-L "${target_merged}")
    fi
    if [[ -n "${germline_name}" ]]; then
      MUTECT_ARGS+=(--germline-resource "${germline_name}")
    fi
    if [[ -n "${pon_name}" ]]; then
      MUTECT_ARGS+=(--panel-of-normals "${pon_name}")
    fi

    if ! "\${MUTECT_ARGS[@]}" > "${sample}.mutect2.stdout.log" 2> "${sample}.mutect2.stderr.log"; then
      echo "ERROR: Mutect2 failed for ${sample}; see ${sample}.mutect2.stderr.log" >&2
      exit 1
    fi

    [[ -s "\$raw" ]] || { echo "ERROR: Mutect2 produced no VCF for ${sample}" >&2; exit 1; }
    [[ -e "\$f1r2" ]] || { echo "ERROR: Mutect2 produced no F1R2 archive for ${sample}" >&2; exit 1; }
    [[ -e "\$stats" ]] || { echo "ERROR: Mutect2 produced no stats file for ${sample}" >&2; exit 1; }
    index_vcf "\$raw"
    echo "completed" > "\$call_status"
    """
}

process MUTECT2_ORIENTATION_MODEL {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 1
    memory '8 GB'
    time params.quick_time
    maxForks params.max_mutect2_parallel
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*orientation_artifact_prior.tar.gz', overwrite: params.overwrite
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*.learn_orientation.*.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(raw_vcf), path(raw_tbi), path(f1r2), path(stats), path(call_status_file), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.unfiltered.vcf.gz"), path("${sample}.unfiltered.vcf.gz.tbi"), path("${sample}.f1r2.tar.gz"), path("${sample}.unfiltered.vcf.gz.stats"), path("${sample}.mutect2_call_status.txt"), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), path("${sample}.orientation_artifact_prior.tar.gz"), emit: artifact
    path "${sample}.learn_orientation.*.log", optional: true, emit: logs

    script:
    """
    set -Eeuo pipefail
    artifact="${sample}.orientation_artifact_prior.tar.gz"
    call_status="\$(cat "${call_status_file}")"
    echo "[MUTECT2_ORIENTATION_MODEL] sample=${sample}"
    if [[ "\$call_status" == "skipped" ]]; then
      : > "\$artifact"
      exit 0
    fi
    [[ -s "${f1r2}" ]] || { echo "ERROR: missing F1R2 evidence for ${sample}" >&2; exit 1; }
    gatk LearnReadOrientationModel \
      -I "${f1r2}" \
      -O "\$artifact" \
      > "${sample}.learn_orientation.stdout.log" \
      2> "${sample}.learn_orientation.stderr.log"
    [[ -s "\$artifact" ]] || { echo "ERROR: orientation model is empty for ${sample}" >&2; exit 1; }
    """
}

process MUTECT2_FILTER {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 2
    memory '12 GB'
    time params.quick_time
    maxForks params.max_mutect2_parallel
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*mutect2.candidate*', overwrite: params.overwrite
    publishDir "${params.outdir}/mutect2", mode: 'copy', pattern: '*.filter_mutect.*.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(raw_vcf), path(raw_tbi), path(f1r2), path(stats), path(call_status_file), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), path(artifact)

    output:
    tuple val(sample), path("${sample}.mutect2.candidate.vcf.gz"), path("${sample}.mutect2.candidate.vcf.gz.tbi"), path("${sample}.mutect2_status.txt"), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), emit: candidate
    path "${sample}.filter_mutect.*.log", optional: true, emit: logs

    script:
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"

    candidate="${sample}.mutect2.candidate.vcf.gz"
    status_file="${sample}.mutect2_status.txt"
    call_status="\$(cat "${call_status_file}")"
    echo "[MUTECT2_FILTER] sample=${sample}"

    [[ -s "${raw_vcf}" ]] || { echo "ERROR: missing raw Mutect2 VCF for ${sample}" >&2; exit 1; }

    if [[ "\$call_status" == "skipped" ]]; then
      cp -f "${raw_vcf}" "\$candidate"
      index_vcf "\$candidate"
      echo "skipped" > "\$status_file"
      exit 0
    fi

    [[ -s "${artifact}" ]] || { echo "ERROR: missing orientation model for ${sample}" >&2; exit 1; }

    FILTER_ARGS=(
      gatk --java-options "-Xmx8g" FilterMutectCalls
      -R "${reference_fasta}"
      -V "${raw_vcf}"
      --orientation-bias-artifact-priors "${artifact}"
      -O "\$candidate"
    )
    if [[ -s "${stats}" ]]; then
      FILTER_ARGS+=(--stats "${stats}")
    fi

    "\${FILTER_ARGS[@]}" > "${sample}.filter_mutect.stdout.log" 2> "${sample}.filter_mutect.stderr.log"
    [[ -s "\$candidate" ]] || { echo "ERROR: FilterMutectCalls produced no VCF for ${sample}" >&2; exit 1; }
    index_vcf "\$candidate"
    echo "filtered" > "\$status_file"
    """
}

process NORMALIZE_MUTECT2 {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.normalize_cpus
    memory params.normalize_memory
    time params.quick_time
    maxForks params.max_normalize_parallel
    publishDir "${params.outdir}/normalized_vcf", mode: 'copy', pattern: '*mutect2.norm.vcf.gz*', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.mutect2.norm.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(candidate_vcf), path(candidate_tbi), path(status_file), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.mutect2.norm.vcf.gz"), path("${sample}.mutect2.norm.vcf.gz.tbi"), path("${sample}.mutect2_status.txt"), emit: norm
    path "${sample}.mutect2.norm.log", optional: true, emit: logs

    script:
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"
    echo "[NORMALIZE_MUTECT2] sample=${sample}"
    normalize_vcf "${candidate_vcf}" "${sample}.mutect2.norm.vcf.gz" "${reference_fasta}" "${sample}.mutect2.norm.log" "${sample}" "MUTECT2_NORMALIZE_EMPTY"
    """
}

process FREEBAYES_CALL {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.freebayes_cpus
    memory params.freebayes_memory
    time params.caller_time
    maxForks params.max_freebayes_parallel
    publishDir "${params.outdir}/freebayes", mode: 'copy', pattern: '*freebayes.vcf.gz*', overwrite: params.overwrite
    publishDir "${params.outdir}/freebayes", mode: 'copy', pattern: '*.freebayes.*.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(bam), path(bai), path(sm_file), val(scope_key), val(scope_unused), path(scope_env), path(target_annot), path(target_merged), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.freebayes.vcf.gz"), path("${sample}.freebayes.vcf.gz.tbi"), path("${sample}.freebayes_status.txt"), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), emit: raw
    path "${sample}.freebayes.*.log", optional: true, emit: logs

    script:
    def skip = truthyParam(params.skip_freebayes) ? 'true' : 'false'
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"
    source "${scope_env}"
    out="${sample}.freebayes.vcf.gz"
    status_file="${sample}.freebayes_status.txt"
    echo "[FREEBAYES_CALL] sample=${sample} target_mode=\$TARGET_MODE"
    if [[ "${skip}" == "true" ]]; then
      write_empty_vcf "${sample}" "\$out" "FREEBAYES_SKIPPED"
      echo "skipped" > "\$status_file"
      exit 0
    fi
    FREEBAYES_ARGS=(
      freebayes
      -f "${reference_fasta}"
      -m "${params.min_mapq}"
      -q "${params.min_bq}"
      --min-alternate-fraction 0.05
      --min-alternate-count 2
    )
    if [[ "\$TARGET_MODE" != "genome_wide" && -s "${target_merged}" ]]; then
      FREEBAYES_ARGS+=(-t "${target_merged}")
    fi
    FREEBAYES_ARGS+=("${bam}")
    set +e
    ( set -o pipefail; "\${FREEBAYES_ARGS[@]}" | bgzip -c > "\$out" ) \
      > "${sample}.freebayes.stdout.log" \
      2> "${sample}.freebayes.stderr.log"
    status=\$?
    set -e
    if [[ "\$status" -ne 0 || ! -s "\$out" ]]; then
      echo "ERROR: FreeBayes failed for ${sample}; see ${sample}.freebayes.stderr.log" >&2
      exit 1
    fi
    index_vcf "\$out"
    echo "completed" > "\$status_file"
    """
}

process NORMALIZE_FREEBAYES {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.normalize_cpus
    memory params.normalize_memory
    time params.quick_time
    maxForks params.max_normalize_parallel
    publishDir "${params.outdir}/normalized_vcf", mode: 'copy', pattern: '*freebayes.norm.vcf.gz*', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.freebayes.norm.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(raw_vcf), path(raw_tbi), path(status_file), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.freebayes.norm.vcf.gz"), path("${sample}.freebayes.norm.vcf.gz.tbi"), path("${sample}.freebayes_status.txt"), emit: norm
    path "${sample}.freebayes.norm.log", optional: true, emit: logs

    script:
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"
    echo "[NORMALIZE_FREEBAYES] sample=${sample}"
    normalize_vcf "${raw_vcf}" "${sample}.freebayes.norm.vcf.gz" "${reference_fasta}" "${sample}.freebayes.norm.log" "${sample}" "FREEBAYES_NORMALIZE_EMPTY"
    """
}

process BCFTOOLS_CALL {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.bcftools_cpus
    memory params.bcftools_memory
    time params.caller_time
    maxForks params.max_bcftools_parallel
    publishDir "${params.outdir}/bcftools", mode: 'copy', pattern: '*bcftools.vcf.gz*', overwrite: params.overwrite
    publishDir "${params.outdir}/bcftools", mode: 'copy', pattern: '*.bcftools.*.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(bam), path(bai), path(sm_file), val(scope_key), val(scope_unused), path(scope_env), path(target_annot), path(target_merged), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.bcftools.vcf.gz"), path("${sample}.bcftools.vcf.gz.tbi"), path("${sample}.bcftools_status.txt"), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready), emit: raw
    path "${sample}.bcftools.*.log", optional: true, emit: logs

    script:
    def skip = truthyParam(params.skip_bcftools) ? 'true' : 'false'
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"
    source "${scope_env}"
    out="${sample}.bcftools.vcf.gz"
    status_file="${sample}.bcftools_status.txt"
    echo "[BCFTOOLS_CALL] sample=${sample} cpus=${task.cpus} target_mode=\$TARGET_MODE"
    if [[ "${skip}" == "true" ]]; then
      write_empty_vcf "${sample}" "\$out" "BCFTOOLS_SKIPPED"
      echo "skipped" > "\$status_file"
      exit 0
    fi
    MPILEUP_ARGS=(
      bcftools mpileup
      -f "${reference_fasta}"
      -q "${params.min_mapq}"
      -Q "${params.min_bq}"
      -a FORMAT/DP,FORMAT/AD
      -Ou
    )
    if [[ "\$TARGET_MODE" != "genome_wide" && -s "${target_merged}" ]]; then
      MPILEUP_ARGS+=(-R "${target_merged}")
    fi
    MPILEUP_ARGS+=("${bam}")
    set +e
    ( set -o pipefail; "\${MPILEUP_ARGS[@]}" | bcftools call -mv --threads ${task.cpus} -Oz -o "\$out" ) \
      > "${sample}.bcftools.stdout.log" \
      2> "${sample}.bcftools.stderr.log"
    status=\$?
    set -e
    if [[ "\$status" -ne 0 || ! -s "\$out" ]]; then
      echo "ERROR: BCFtools failed for ${sample}; see ${sample}.bcftools.stderr.log" >&2
      exit 1
    fi
    index_vcf "\$out"
    echo "completed" > "\$status_file"
    """
}

process NORMALIZE_BCFTOOLS {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.normalize_cpus
    memory params.normalize_memory
    time params.quick_time
    maxForks params.max_normalize_parallel
    publishDir "${params.outdir}/normalized_vcf", mode: 'copy', pattern: '*bcftools.norm.vcf.gz*', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.bcftools.norm.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(raw_vcf), path(raw_tbi), path(status_file), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.bcftools.norm.vcf.gz"), path("${sample}.bcftools.norm.vcf.gz.tbi"), path("${sample}.bcftools_status.txt"), emit: norm
    path "${sample}.bcftools.norm.log", optional: true, emit: logs

    script:
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"
    echo "[NORMALIZE_BCFTOOLS] sample=${sample}"
    normalize_vcf "${raw_vcf}" "${sample}.bcftools.norm.vcf.gz" "${reference_fasta}" "${sample}.bcftools.norm.log" "${sample}" "BCFTOOLS_NORMALIZE_EMPTY"
    """
}


process CALLER_READY {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus 1
    memory '2 GB'
    time '1h'
    maxForks params.max_finalize_parallel
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.callers.ready.tsv', overwrite: params.overwrite

    input:
    tuple val(sample), path(mutect_vcf), path(mutect_tbi), path(mutect_status_file), path(freebayes_vcf), path(freebayes_tbi), path(freebayes_status_file), path(bcftools_vcf), path(bcftools_tbi), path(bcftools_status_file)

    output:
    tuple val(sample), path(mutect_vcf), path(mutect_tbi), path(mutect_status_file), path(freebayes_vcf), path(freebayes_tbi), path(freebayes_status_file), path(bcftools_vcf), path(bcftools_tbi), path(bcftools_status_file), path("${sample}.callers.ready.tsv"), emit: ready

    script:
    """
    set -Eeuo pipefail
    echo "[CALLER_READY] sample=${sample}: verifying Mutect2 + FreeBayes + BCFtools normalized VCFs are present before merge/nested classifier"

    for f in "${mutect_vcf}" "${mutect_tbi}" "${mutect_status_file}" "${freebayes_vcf}" "${freebayes_tbi}" "${freebayes_status_file}" "${bcftools_vcf}" "${bcftools_tbi}" "${bcftools_status_file}"; do
      if [[ ! -s "\$f" ]]; then
        echo "ERROR: caller output missing or empty for sample=${sample}: \$f" >&2
        exit 1
      fi
    done

    {
      printf 'sample\tmutect2_vcf\tmutect2_status\tfreebayes_vcf\tfreebayes_status\tbcftools_vcf\tbcftools_status\tstatus\n'
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${sample}" \
        "${mutect_vcf}" \
        "\$(cat "${mutect_status_file}" 2>/dev/null || echo unknown)" \
        "${freebayes_vcf}" \
        "\$(cat "${freebayes_status_file}" 2>/dev/null || echo unknown)" \
        "${bcftools_vcf}" \
        "\$(cat "${bcftools_status_file}" 2>/dev/null || echo unknown)" \
        "caller_stages_finished"
    } > "${sample}.callers.ready.tsv"
    cat "${sample}.callers.ready.tsv"
    """
}

process MERGE_CALLERS_PER_SAMPLE {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.finalize_cpus
    memory params.finalize_memory
    time params.quick_time
    maxForks params.max_finalize_parallel
    publishDir "${params.outdir}/reports/per_sample", mode: 'copy', pattern: '*.final_variants.*.tsv', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.merge_callers.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(mutect_vcf), path(mutect_tbi), path(mutect_status_file), path(freebayes_vcf), path(freebayes_tbi), path(freebayes_status_file), path(bcftools_vcf), path(bcftools_tbi), path(bcftools_status_file), path(callers_ready), path(thresholds_tsv), val(scope_key), val(scope_unused), path(scope_env), path(target_annot), path(target_merged), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.final_variants.pre_norm.vcf"), emit: pre_vcf
    path "${sample}.final_variants.*.tsv", emit: reports
    path "${sample}.merge_callers.log", optional: true, emit: logs

    script:
    def includeArtifacts = truthyParam(params.include_obvious_artifacts) ? '--include-obvious-artifacts' : ''
    def sampleMode = truthyParam(params.ffpe) ? 'ffpe' : 'fresh'
    """
    set -Eeuo pipefail
    source "${scope_env}"
    mutect_status="\$(cat "${mutect_status_file}" 2>/dev/null || echo not_run)"
    freebayes_status="\$(cat "${freebayes_status_file}" 2>/dev/null || echo not_run)"
    bcftools_status="\$(cat "${bcftools_status_file}" 2>/dev/null || echo not_run)"
    TARGETS_ARG=""
    if [[ "\$TARGET_MODE" != "genome_wide" && -s "${target_annot}" ]]; then
      TARGETS_ARG="${target_annot}"
    fi

    echo "[MERGE_CALLERS_PER_SAMPLE] sample=${sample} mutect_status=\$mutect_status freebayes_status=\$freebayes_status bcftools_status=\$bcftools_status target_mode=\$TARGET_MODE" | tee "${sample}.merge_callers.log"
    echo "[MERGE_CALLERS_PER_SAMPLE] caller_gate=${callers_ready}" | tee -a "${sample}.merge_callers.log"
    cat "${callers_ready}" >> "${sample}.merge_callers.log" || true

    if [[ -s "${thresholds_tsv}" ]]; then
      read -r _header < "${thresholds_tsv}"
      tab_char="\$(printf '\t')"
      if [[ "${sampleMode}" == "ffpe" ]]; then
        IFS="\$tab_char" read -r th_sample th_target_mode th_cov_source th_mean_depth th_mode th_min_dp th_min_alt th_min_af th_vote th_note th_ffpe_vaf th_ffpe_alt < <(tail -n +2 "${thresholds_tsv}" | head -1)
      else
        IFS="\$tab_char" read -r th_sample th_target_mode th_cov_source th_mean_depth th_mode th_min_dp th_min_alt th_min_af th_vote th_note < <(tail -n +2 "${thresholds_tsv}" | head -1)
      fi
    else
      th_mean_depth="."
      th_mode="manual_fallback_missing_threshold_file"
      th_min_dp="${params.min_dp}"
      th_min_alt="${params.min_alt_reads}"
      th_min_af="${params.min_af}"
      th_vote="${params.vote_threshold}"
      if [[ "${sampleMode}" == "ffpe" ]]; then
        th_ffpe_vaf="${params.ffpe_ct_vaf_keep}"
        th_ffpe_alt="${params.ffpe_ct_min_alt}"
      else
        th_ffpe_vaf="not_applicable"
        th_ffpe_alt="not_applicable"
      fi
      th_note="missing_threshold_file"
    fi

    echo "[MERGE_CALLERS_PER_SAMPLE] effective_thresholds sample=${sample} sample_mode=${sampleMode} mean_depth=\$th_mean_depth mode=\$th_mode min_dp=\$th_min_dp min_alt_reads=\$th_min_alt min_af=\$th_min_af vote_threshold=\$th_vote note=\$th_note" | tee -a "${sample}.merge_callers.log"
    if [[ "${sampleMode}" == "ffpe" ]]; then
      echo "[MERGE_CALLERS_PER_SAMPLE] ffpe_thresholds sample=${sample} ffpe_ct_vaf_keep=\$th_ffpe_vaf ffpe_ct_min_alt=\$th_ffpe_alt" | tee -a "${sample}.merge_callers.log"
    fi

    FINALIZER_ARGS=(
      python "${projectDir}/bin/build_final_per_sample_vcf.py"
      --sample "${sample}"
      --sample-mode "${sampleMode}"
      --mutect-vcf "${mutect_vcf}"
      --mutect-status "\$mutect_status"
      --freebayes-vcf "${freebayes_vcf}"
      --freebayes-status "\$freebayes_status"
      --bcftools-vcf "${bcftools_vcf}"
      --bcftools-status "\$bcftools_status"
      --targets "\$TARGETS_ARG"
      --ref-fai "${reference_fai}"
      --outdir .
      --out-vcf "${sample}.final_variants.pre_norm.vcf"
      --min-dp "\$th_min_dp"
      --min-alt-reads "\$th_min_alt"
      --min-af "\$th_min_af"
      --vote-threshold "\$th_vote"
    )
    if [[ "${sampleMode}" == "ffpe" ]]; then
      FINALIZER_ARGS+=(
        --ffpe-ct-vaf-keep "\$th_ffpe_vaf"
        --ffpe-ct-min-alt "\$th_ffpe_alt"
      )
    fi
    if [[ -n "${includeArtifacts}" ]]; then
      FINALIZER_ARGS+=("${includeArtifacts}")
    fi
    "\${FINALIZER_ARGS[@]}" >> "${sample}.merge_callers.log" 2>&1
    """
}

process NORMALIZE_FINAL_VCF {
    tag "${sample}"
    conda "${projectDir}/envs/lowpass_variants.yml"
    cpus params.normalize_cpus
    memory params.normalize_memory
    time params.quick_time
    maxForks params.max_finalize_parallel
    publishDir "${params.outdir}/final_vcf", mode: 'copy', pattern: '*.final_variants.vcf*', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.final_vcf.norm.log', overwrite: params.overwrite

    input:
    tuple val(sample), path(pre_vcf), path(reference_fasta), path(reference_fai), path(reference_dict), path(reference_ready)

    output:
    tuple val(sample), path("${sample}.final_variants.vcf"), path("${sample}.final_variants.vcf.gz"), path("${sample}.final_variants.vcf.gz.tbi"), emit: final_vcf
    path "${sample}.final_vcf.norm.log", optional: true, emit: logs

    script:
    """
    set -Eeuo pipefail
    source "${projectDir}/bin/nf_common.sh"
    echo "[NORMALIZE_FINAL_VCF] sample=${sample}"
    normalize_vcf \
      "${pre_vcf}" \
      "${sample}.final_variants.vcf.gz" \
      "${reference_fasta}" \
      "${sample}.final_vcf.norm.log" \
      "${sample}" \
      "FINAL_NORMALIZE_EMPTY"
    bcftools view "${sample}.final_variants.vcf.gz" > "${sample}.final_variants.vcf"
    bcftools index --stats "${sample}.final_variants.vcf.gz" >/dev/null
    """
}

process FFPERASE {
    tag 'ffperase_cohort'
    conda "${projectDir}/envs/ffperase_runner.yml"
    cpus params.ffperase_threads
    memory params.ffperase_memory
    time params.ffperase_time
    maxForks 1
    publishDir "${params.outdir}", mode: 'copy', pattern: 'ffperase_classification', overwrite: params.overwrite
    publishDir "${params.outdir}", mode: 'copy', pattern: 'ffperase_status', overwrite: params.overwrite
    publishDir "${params.outdir}", mode: 'copy', pattern: 'ffperase_picard_metrics_summary', overwrite: params.overwrite
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'ffperase.done.txt', overwrite: params.overwrite

    input:
    path final_files
    path bam_files
    path threshold_files
    tuple path(reference_fasta), path(reference_fai), path(reference_dict)

    output:
    path 'ffperase_classification', emit: annotated
    path 'ffperase_status', emit: status
    path 'ffperase_picard_metrics_summary', emit: picard_metrics
    path 'ffperase.done.txt', emit: done

    script:
    def filterArtifacts = truthyParam(params.ffperase_filter_artifacts) ? '--filter-artifacts' : ''
    def failOnError = truthyParam(params.ffperase_fail_on_error) ? 'true' : 'false'
    """
    set -Eeuo pipefail
    shopt -s nullglob

    INPUT_ROOT="\$PWD/ffperase_input"
    mkdir -p "\$INPUT_ROOT/final_vcf" "\$INPUT_ROOT/preprocessed_bam" "\$INPUT_ROOT/coverage"

    stage_regular_file() {
      local src="\$1"
      local dest="\$2"
      # A hard link is a regular file and avoids duplicating WGS BAM data. If
      # work/output paths span filesystems, fall back to a dereferenced copy.
      ln -L "\$src" "\$dest" 2>/dev/null || cp -L "\$src" "\$dest"
    }

    final_count=0
    for f in *.final_variants.vcf.gz *.final_variants.vcf.gz.tbi; do
      # Use regular files here: the FFPErase adapter discovers inputs with
      # `find -type f`, while Nextflow normally stages process inputs as links.
      stage_regular_file "\$PWD/\$f" "\$INPUT_ROOT/final_vcf/\$f"
      if [[ "\$f" == *.final_variants.vcf.gz ]]; then
        final_count=\$((final_count + 1))
      fi
    done

    bam_count=0
    for f in *.preprocessed.bam *.preprocessed.bam.bai; do
      stage_regular_file "\$PWD/\$f" "\$INPUT_ROOT/preprocessed_bam/\$f"
      if [[ "\$f" == *.preprocessed.bam ]]; then
        bam_count=\$((bam_count + 1))
      fi
    done

    [[ "\$final_count" -gt 0 ]] || {
      echo "ERROR: no final VCFs were staged for FFPErase" >&2
      exit 1
    }
    [[ "\$bam_count" -eq "\$final_count" ]] || {
      echo "ERROR: FFPErase input mismatch: final_vcfs=\$final_count bams=\$bam_count" >&2
      exit 1
    }

    threshold_count=0
    for f in *.adaptive_thresholds.tsv; do
      stage_regular_file "\$PWD/\$f" "\$INPUT_ROOT/coverage/\$f"
      threshold_count=\$((threshold_count + 1))
    done
    [[ "\$threshold_count" -eq "\$final_count" ]] || {
      echo "ERROR: FFPErase input mismatch: final_vcfs=\$final_count adaptive_thresholds=\$threshold_count" >&2
      exit 1
    }

    bash "${projectDir}/bin/run_ffperase_single_picard_pileup_nf.sh" \
      --root "\$INPUT_ROOT" \
      --ref "${reference_fasta}" \
      --outroot "\$PWD/ffperase_post" \
      --picard-metrics-root "\$PWD/ffperase_picard_metrics" \
      --features-root "\$PWD/ffperase_post/features_cache" \
      --models-dir "\$PWD/ffperase_post/models" \
      --threads "${params.ffperase_threads}" \
      --classify-threads "${params.ffperase_classify_threads}" \
      --sample-jobs "${params.ffperase_sample_jobs}" \
      --picard-jobs "${params.ffperase_picard_jobs}" \
      --split-pileup "${params.ffperase_split_pileup}" \
      --split-reads "${params.ffperase_split_reads}" \
      --container "${params.ffperase_container}" \
      --engine "${params.ffperase_engine}" \
      --repository "${params.ffperase_repository}" \
      --revision "${params.ffperase_revision}" \
      --fail-on-error "${failOnError}" \
      --no-conda \
      --no-create-conda-env \
      --gatk-bin "\$(command -v gatk)"

    python "${projectDir}/bin/annotate_vcfs_with_ffperase.py" \
      --root "\$INPUT_ROOT" \
      --final-vcf-dir "\$INPUT_ROOT/final_vcf" \
      --ffperase-dir "\$PWD/ffperase_post" \
      --outdir "\$PWD/ffperase_classification" \
      ${filterArtifacts}

    # Publish only durable results. Nested work, feature caches, models, input
    # BAM links, and reference copies remain in the outer task for -resume.
    mkdir -p ffperase_status/logs ffperase_picard_metrics_summary
    if [[ -d ffperase_post/status ]]; then
      cp -a ffperase_post/status/. ffperase_status/
    fi
    for f in ffperase_post/missing_*.tsv ffperase_post/*.log; do
      [[ -f "\$f" ]] && cp -a "\$f" ffperase_status/
    done
    for f in ffperase_post/*/snvs/*.nf_ffperase_*.console.log ffperase_post/*/indels/*.nf_ffperase_*.console.log; do
      [[ -f "\$f" ]] && cp -a "\$f" ffperase_status/logs/
    done
    for f in ffperase_picard_metrics/*/pre_adapter_metrics.tsv ffperase_picard_metrics/*/bait_bias_metrics.tsv; do
      [[ -f "\$f" ]] || continue
      sample_dir="ffperase_picard_metrics_summary/\$(basename "\$(dirname "\$f")")"
      mkdir -p "\$sample_dir"
      cp -a "\$f" "\$sample_dir/"
    done
    for f in ffperase_picard_metrics/picard_metrics_manifest.tsv ffperase_picard_metrics/*.log; do
      [[ -f "\$f" ]] && cp -a "\$f" ffperase_picard_metrics_summary/
    done

    failed_count=0
    status_file_count=0
    for f in ffperase_status/*.status.tsv; do
      status_file_count=\$((status_file_count + 1))
      if awk -F '\t' 'NR == 2 && \$3 ~ /^FAILED/ { failed=1 } END { exit(failed ? 0 : 1) }' "\$f"; then
        failed_count=\$((failed_count + 1))
      fi
    done
    completion_status="completed"
    if [[ "\$failed_count" -gt 0 ]]; then
      completion_status="completed_with_failures"
    fi

    {
      printf 'mode\tffpe\n'
      printf 'samples\t%s\n' "\$final_count"
      printf 'sample_type_status_files\t%s\n' "\$status_file_count"
      printf 'failed_sample_types\t%s\n' "\$failed_count"
      printf 'status\t%s\n' "\$completion_status"
    } > ffperase.done.txt
    """
}


workflow {
    if( truthyParam(params.help) || truthyParam(params.show_help) || truthyParam(params.h) || truthyParam(params.help_full) ) {
        println workflowHelp(truthyParam(params.help_full), truthyParam(params.show_hidden) || truthyParam(params.showHidden), "${projectDir}/assets")
        System.exit(0)
    }

    ffpe_mode = truthyParam(params.ffpe)
    fresh_mode = truthyParam(params.fresh)
    if( ffpe_mode == fresh_mode ) {
        error "Select exactly one sample mode: --ffpe or --fresh"
    }
    if( !params.input ) {
        error "Missing required parameter: --input"
    }
    if( !params.outdir ) {
        error "Missing required parameter: --outdir"
    }
    if( !params.ref ) {
        error "Missing required parameter: --ref"
    }
    if( !(params.duplicate_mode in ['skip', 'mark', 'remove']) ) {
        error "--duplicate_mode must be one of: skip, mark, remove"
    }
    if( !(params.cohort in ['all', 'lymphoma', 'brain']) ) {
        error "--cohort must be one of: all, lymphoma, brain"
    }
    if( !truthyParam(params.dry_run) && !truthyParam(params.skip_bqsr) && !(params.known_sites?.toString()?.trim()) ) {
        error "BQSR is enabled by default: provide --known_sites or use --skip_bqsr"
    }
    if( threadsParam() < 1 || mutect2Cpus() < 1 || freebayesCpus() < 1 || bcftoolsCpus() < 1 ) {
        error "CPU parameters must be >= 1"
    }
    if( baseParallel() < 1 || mutect2Parallel() < 1 || freebayesParallel() < 1 || bcftoolsParallel() < 1 ) {
        error "Parallelism parameters must be >= 1"
    }

    RESOLVE_BAMS()

    samples_ch = RESOLVE_BAMS.out.samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample as String, file(row.bam as String, checkIfExists: true)) }

    if ( truthyParam(params.dry_run) ) {
        samples_ch.view { row -> "DRY_RUN_SAMPLE\t${row[0]}\t${row[1]}" }
    } else {
        reference_source_ch = channel.fromPath(params.ref.toString(), checkIfExists: true)
        PREPARE_REFERENCE(reference_source_ch)

        placeholder_file = file("${projectDir}/assets/no_resource.placeholder", checkIfExists: true)
        scope_bed_name = ''
        scope_gtf_name = ''
        scope_resource_paths = []
        scope_bed_raw = (params.bed ?: params.target_bed ?: '').toString().trim()
        scope_gtf_raw = (params.gtf ?: '').toString().trim()
        if( scope_bed_raw ) {
            scope_bed_file = file(scope_bed_raw, checkIfExists: true)
            scope_bed_name = scope_bed_file.name
            scope_resource_paths << scope_bed_file
        }
        if( scope_gtf_raw ) {
            scope_gtf_file = file(scope_gtf_raw, checkIfExists: true)
            scope_gtf_name = scope_gtf_file.name
            scope_resource_paths << scope_gtf_file
        }
        if( scope_resource_paths.isEmpty() ) {
            scope_resource_paths << placeholder_file
        }
        scope_stage_names = scope_resource_paths.collect { scope_resource -> scope_resource.name }
        if( scope_stage_names.size() != scope_stage_names.toSet().size() ) {
            error "BED/GTF resources have colliding basenames; use uniquely named files"
        }
        scope_resources_ch = channel.value(tuple(scope_bed_name, scope_gtf_name, scope_resource_paths))
        PREPARE_SCOPE(scope_resources_ch)

        known_site_names = []
        known_site_paths = []
        if ( truthyParam(params.skip_bqsr) ) {
            known_site_paths << placeholder_file
        } else {
            params.known_sites.toString().split(',')
                .collect { raw_site -> raw_site.trim() }
                .findAll { raw_site -> raw_site }
                .each { raw_path ->
                def known_vcf = file(raw_path, checkIfExists: true)
                def known_path = known_vcf.toString()
                def index_path
                index_path = ["${known_path}.tbi", "${known_path}.csi", "${known_path}.idx"]
                    .find { candidate -> new java.io.File(candidate).isFile() && new java.io.File(candidate).length() > 0 }
                if( !index_path ) {
                    error "Indexed BQSR known-sites VCF required (.tbi/.csi/.idx): ${raw_path}"
                }
                known_site_names << known_vcf.name
                known_site_paths << known_vcf
                known_site_paths << file(index_path, checkIfExists: true)
            }
        }
        known_stage_names = known_site_paths.collect { known_resource -> known_resource.name }
        if( known_stage_names.size() != known_stage_names.toSet().size() ) {
            error "BQSR resources have colliding basenames; use uniquely named VCF/index files"
        }
        known_sites_bundle_ch = channel.value(tuple(known_site_names.join(','), known_site_paths))

        germline_name = ''
        pon_name = ''
        mutect_resource_paths = []
        if ( !truthyParam(params.skip_mutect2) ) {
            pon_raw = params.panel_of_normals ?: params.pon ?: ''
            mutect_specs = [
                [kind: 'germline', raw: params.germline_resource ?: ''],
                [kind: 'pon', raw: pon_raw]
            ]
            mutect_specs.findAll { spec -> spec.raw.toString().trim() }.each { spec ->
                def resource_vcf = file(spec.raw.toString(), checkIfExists: true)
                def resource_path = resource_vcf.toString()
                def resource_index
                resource_index = ["${resource_path}.tbi", "${resource_path}.csi", "${resource_path}.idx"]
                    .find { candidate -> new java.io.File(candidate).isFile() && new java.io.File(candidate).length() > 0 }
                if( !resource_index ) {
                    error "Indexed Mutect2 ${spec.kind} VCF required (.tbi/.csi/.idx): ${spec.raw}"
                }
                if( spec.kind == 'germline' ) germline_name = resource_vcf.name
                if( spec.kind == 'pon' ) pon_name = resource_vcf.name
                mutect_resource_paths << resource_vcf
                mutect_resource_paths << file(resource_index, checkIfExists: true)
            }
        }
        if( mutect_resource_paths.isEmpty() ) {
            mutect_resource_paths << placeholder_file
        }
        mutect_resource_paths = mutect_resource_paths.unique()
        mutect_stage_names = mutect_resource_paths.collect { mutect_resource -> mutect_resource.name }
        if( mutect_stage_names.size() != mutect_stage_names.toSet().size() ) {
            error "Mutect2 resources have colliding basenames; use uniquely named VCF/index files"
        }
        mutect_resources_ch = channel.value(tuple(germline_name, pon_name, mutect_resource_paths))

        sample_scope_ref_ch = samples_ch
            .combine(PREPARE_SCOPE.out.scope_tuple)
            .combine(PREPARE_REFERENCE.out.bundle)
            .combine(known_sites_bundle_ch)
            .map { sample, bam, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready, known_sites_names, known_site_files ->
                tuple(sample, bam, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready, known_sites_names, known_site_files)
            }

        PREPARE_BAM(sample_scope_ref_ch)

        prepared_scope_ch = PREPARE_BAM.out.prepared
            .combine(PREPARE_SCOPE.out.scope_tuple)
            .map { sample, bam, bai, sm_file, scope_key, scope_unused, scope_env, target_annot, target_merged ->
                tuple(sample, bam, bai, sm_file, scope_key, scope_unused, scope_env, target_annot, target_merged)
            }

        COVERAGE(prepared_scope_ch)
        CALIBRATE_THRESHOLDS(COVERAGE.out.coverage_bundle)

        prepared_scope_ref_ch = prepared_scope_ch
            .combine(PREPARE_REFERENCE.out.bundle)
            .map { sample, bam, bai, sm_file, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready ->
                tuple(sample, bam, bai, sm_file, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready)
            }

        mutect_input_ch = prepared_scope_ref_ch
            .combine(mutect_resources_ch)
            .map { sample, bam, bai, sm_file, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready, germline_resource_name, pon_resource_name, mutect_resource_files ->
                tuple(sample, bam, bai, sm_file, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready, germline_resource_name, pon_resource_name, mutect_resource_files)
            }

        MUTECT2_CALL(mutect_input_ch)
        MUTECT2_ORIENTATION_MODEL(MUTECT2_CALL.out.raw)
        MUTECT2_FILTER(MUTECT2_ORIENTATION_MODEL.out.artifact)
        NORMALIZE_MUTECT2(MUTECT2_FILTER.out.candidate)

        FREEBAYES_CALL(prepared_scope_ref_ch)
        NORMALIZE_FREEBAYES(FREEBAYES_CALL.out.raw)

        BCFTOOLS_CALL(prepared_scope_ref_ch)
        NORMALIZE_BCFTOOLS(BCFTOOLS_CALL.out.raw)

        caller_norm_ch = NORMALIZE_MUTECT2.out.norm
            .join(NORMALIZE_FREEBAYES.out.norm)
            .join(NORMALIZE_BCFTOOLS.out.norm)

        // Explicit per-sample barrier: no merge or post-classification can start until
        // Mutect2, FreeBayes and BCFtools normalized outputs are all present for that sample.
        CALLER_READY(caller_norm_ch)

        caller_threshold_ch = CALLER_READY.out.ready
            .join(CALIBRATE_THRESHOLDS.out.thresholds)

        caller_scope_ch = caller_threshold_ch
            .combine(PREPARE_SCOPE.out.scope_tuple)
            .combine(PREPARE_REFERENCE.out.bundle)
            .map { sample, mut_vcf, mut_tbi, mut_status, fb_vcf, fb_tbi, fb_status, bcf_vcf, bcf_tbi, bcf_status, callers_ready, thresholds_tsv, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready ->
                tuple(sample, mut_vcf, mut_tbi, mut_status, fb_vcf, fb_tbi, fb_status, bcf_vcf, bcf_tbi, bcf_status, callers_ready, thresholds_tsv, scope_key, scope_unused, scope_env, target_annot, target_merged, reference_fasta, reference_fai, reference_dict, reference_ready)
            }

        MERGE_CALLERS_PER_SAMPLE(caller_scope_ch)
        final_pre_ref_ch = MERGE_CALLERS_PER_SAMPLE.out.pre_vcf
            .combine(PREPARE_REFERENCE.out.bundle)
            .map { sample, pre_vcf, reference_fasta, reference_fai, reference_dict, reference_ready ->
                tuple(sample, pre_vcf, reference_fasta, reference_fai, reference_dict, reference_ready)
            }
        NORMALIZE_FINAL_VCF(final_pre_ref_ch)

        if ( ffpe_mode ) {
            ffperase_final_files_ch = NORMALIZE_FINAL_VCF.out.final_vcf
                .map { _sample, _plain_vcf, final_vcf, final_tbi -> [final_vcf, final_tbi] }
                .flatten()
                .collect()

            ffperase_bam_files_ch = PREPARE_BAM.out.prepared
                .map { _sample, bam, bai, _sm_file -> [bam, bai] }
                .flatten()
                .collect()

            ffperase_threshold_files_ch = CALIBRATE_THRESHOLDS.out.thresholds
                .map { _sample, thresholds_tsv -> thresholds_tsv }
                .collect()

            ffperase_reference_ch = PREPARE_REFERENCE.out.bundle
                .map { reference_fasta, reference_fai, reference_dict, _reference_ready ->
                    tuple(reference_fasta, reference_fai, reference_dict)
                }

            FFPERASE(ffperase_final_files_ch, ffperase_bam_files_ch, ffperase_threshold_files_ch, ffperase_reference_ch)
        }
    }
}
