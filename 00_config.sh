#!/usr/bin/env bash
# =============================================================================
# 00_config.sh
# Shared configuration sourced by every step script.
# Edit this file before running anything else.
# =============================================================================

# ---- Paths ------------------------------------------------------------------
export PROJECT_DIR="/home/projects-phoenix/ApoE_Carnitine/MetaG"
export RAW_DIR="${PROJECT_DIR}/raw_reads"

# ---- Scripts directory ------------------------------------------------------
# This must be the directory where all the pipeline scripts AND sample_list.txt live.
export SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# ---- Sample list ------------------------------------------------------------
# Reads sample names from sample_list.txt in the scripts directory.
# One sample name per line. Generate it with the command in the README.
# Individual scripts can still be overridden by passing sample names as arguments.
SAMPLE_LIST="${SCRIPTS_DIR}/sample_list.txt"
if [[ -f "${SAMPLE_LIST}" ]]; then
    export SAMPLES="${SAMPLES:-$(cat "${SAMPLE_LIST}" | tr '\n' ' ')}"
else
    export SAMPLES="${SAMPLES:-}"
fi

# Just ApoE sample list, generated manually
export APOE_SAMPLES="${APOE_SAMPLES:-$(cat "${SCRIPTS_DIR}/apoe_sample_list.txt" | tr '\n' ' ')}"

# ---- Catalog samples (control + carnitine only) -------
CONTROL_LIST="${SCRIPTS_DIR}/control_list.txt"
CARNITINE_LIST="${SCRIPTS_DIR}/carnitine_list.txt"

export CATALOG_SAMPLES=$(awk 1 "${CONTROL_LIST}" "${CARNITINE_LIST}" | tr '\n' ' ')

# ---- Global resources -------------------------------------------------------
export THREADS="50"
export MEMORY="400G"

# ---- Utility scripts (Borton Lab ORG-Data) ----------------------------------
export CONTIG_STATS="/ORG-Data/scripts/quicklooks/contig_stats.pl"
export RENAME_BINS_SCRIPT="rename_bins_like_dram.py"

# ---- Assembly ---------------------------------------------------------------
export KMER_MIN="31"
export KMER_MAX="121"
export KMER_STEP="10"
export MIN_CONTIG_LEN="2500"
export POLYG_TRIM="50"

# ---- Mapping identity filters -----------------------------------------------
export BBMAP_MIN_ID="0.99"
export BOWTIE2_MIN_ID="0.97"
export COVERM_MIN_ID="0.97"
export COVERM_MIN_BREADTH="0.75"

# ---- MAG quality thresholds -------------------------------------------------
export CHECKM2_MIN_COMP="50"
export CHECKM2_MAX_CONT="10"
export DREP_MIN_COMP="50"
export DREP_MAX_CONT="10"
# ANI threshold for dRep clustering. 99% = strain-level dereplication (recommended
# for building a high-resolution MAG database). Default dRep is 95%.
export DREP_ANI="0.99"

# ---- Per-tool thread overrides ----------------------------------------------
export CHECKM2_THREADS="10"
export GTDBTK_CPUS="20"
export BOWTIE2_THREADS="15"

# ---- Metadata file (project-level, written to by multiple steps) ------------
export METADATA="${PROJECT_DIR}/pipeline_metadata.tsv"

# ---- Helper: resolve a raw read path that may be .fastq or .fastq.gz -------
find_raw_read() {
    local base="$1"
    if   [[ -f "${base}.fastq" ]];    then echo "${base}.fastq"
    elif [[ -f "${base}.fastq.gz" ]]; then echo "${base}.fastq.gz"
    else echo ""
    fi
}

# ---- Derived per-sample path helpers ----------------------------------------
set_sample_dirs() {
    local s="$1"
    export SAMPLE_DIR="${PROJECT_DIR}/${s}"
    export TRIMMED_DIR="${SAMPLE_DIR}/trimmed_reads"
    export MEGAHIT_DIR="${SAMPLE_DIR}/megahit_out"
    export MAGS_DIR="${SAMPLE_DIR}/MAGs"
    export GTDB_DIR="${SAMPLE_DIR}/gtdb_v2.7.0_r232"
    export MAP_DIR="${SAMPLE_DIR}/MAG_db_mapping"

    export R1_RAW=$(find_raw_read "${RAW_DIR}/${s}_R1_001")
    export R2_RAW=$(find_raw_read "${RAW_DIR}/${s}_R2_001")

    export SICKLE_R1="${TRIMMED_DIR}/${s}_R1_sickle_trimmed.fastq"
    export SICKLE_R2="${TRIMMED_DIR}/${s}_R2_sickle_trimmed.fastq"
    export TRIMMED_R1="${TRIMMED_DIR}/${s}_R1_trimmed_noplyG.fastq"
    export TRIMMED_R2="${TRIMMED_DIR}/${s}_R2_trimmed_noplyG.fastq"
    export DISCARDED="${TRIMMED_DIR}/${s}_discarded.fastq"

    export FILTERED_SCAFFOLDS="${MEGAHIT_DIR}/${s}_${MIN_CONTIG_LEN}.fa"
    export BINS_DIR="${FILTERED_SCAFFOLDS}.metabat-bins"
    export CHECKM2_OUT="${MEGAHIT_DIR}/checkm2_v1.1.0_fa"
    export CHECKM2_RESULTS="${CHECKM2_OUT}/quality_report.tsv"
}

# ---- Cross-sample paths -----------------------------------------------------
export CROSS_DREP_DIR="${PROJECT_DIR}/cross_sample_dRep"
export CROSS_DEREP_GENOMES="${CROSS_DREP_DIR}/dereplicated_genomes"
export CROSS_RENAMED_DIR="${CROSS_DEREP_GENOMES}/genome_renamed"
export CROSS_BOWTIE_DIR="${PROJECT_DIR}/cross_sample_bowtie_DB"
export CROSS_MAG_DB_FA="${CROSS_BOWTIE_DIR}/all_samples_derep_MAGs.fa"
export COASSEMBLY_DIR="${PROJECT_DIR}/coassemblies"
export POOL_DIR="${PROJECT_DIR}/ALL_MQHQ_MAGs"
export CROSS_BOWTIE_PREFIX="${CROSS_BOWTIE_DIR}/all_samples_derep_MAG_DB"
export COVERM_OUT_DIR="${PROJECT_DIR}/coverm_output"
export MAPPING_SUMMARY="${PROJECT_DIR}/sample_MAG_database_mapping_summary.tsv"
export CON_CARN_MAPPING_SUMMARY="${PROJECT_DIR}/con-carn_MAG_database_mapping_summary.tsv"
