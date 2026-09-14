#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
#SBATCH --job-name=coverm
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_11_coverm_%j.out
#SBATCH --error=slurm_11_coverm_%j.err
# =============================================================================
# 11_coverm.sh  —  Step 17: Between-sample CoverM abundance tables.
# Skips individual metrics if their output file already exists.
# Usage:  sbatch 11_coverm.sh

#This version uses only the carnitine + control possort bams and the carnitine + control dRep MAG database
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

BAM_LIST="${PROJECT_DIR}/possort_bam_list.txt"

if [[ ! -f "${BAM_LIST}" ]]; then
    echo "ERROR: BAM list not found at ${BAM_LIST}. Run 10_bowtie2_map.sh first."
    exit 1
fi
if [[ ! -d "${CROSS_RENAMED_DIR}" ]]; then
    echo "ERROR: cross-sample renamed genomes directory not found. Run 09_cross_drep.sh first."
    exit 1
fi

mapfile -t ALL_BAMS < "${BAM_LIST}"
if [[ "${#ALL_BAMS[@]}" -eq 0 ]]; then
    echo "ERROR: BAM list is empty."
    exit 1
fi

mkdir -p "${COVERM_OUT_DIR}"
echo "Running CoverM with ${#ALL_BAMS[@]} BAM file(s)..."

# ---- reads_per_base ---------------------------------------------------------
if [[ -f "${COVERM_OUT_DIR}/coverm_reads_per_base.txt" ]]; then
    echo "  coverm_reads_per_base.txt already exists — skipping."
else
    echo "  [CoverM] reads_per_base..."
    coverm genome \
        --proper-pairs-only \
        --genome-fasta-extension fa \
        --genome-fasta-directory "${CROSS_RENAMED_DIR}" \
        --bam-files "${ALL_BAMS[@]}" \
        --threads "${BOWTIE2_THREADS}" \
        --min-read-percent-identity-pair "${COVERM_MIN_ID}" \
        --min-covered-fraction 0 \
        -m reads_per_base \
        --output-file "${COVERM_OUT_DIR}/coverm_reads_per_base.txt" \
        2> "${COVERM_OUT_DIR}/reads_per_base_stats.txt"
fi

# ---- min covered fraction ---------------------------------------------------
if [[ -f "${COVERM_OUT_DIR}/coverm_min75.txt" ]]; then
    echo "  coverm_min75.txt already exists — skipping."
else
    echo "  [CoverM] min covered fraction >= ${COVERM_MIN_BREADTH}..."
    coverm genome \
        --proper-pairs-only \
        --genome-fasta-extension fa \
        --genome-fasta-directory "${CROSS_RENAMED_DIR}" \
        --bam-files "${ALL_BAMS[@]}" \
        --threads "${BOWTIE2_THREADS}" \
        --min-read-percent-identity-pair "${COVERM_MIN_ID}" \
        --min-covered-fraction "${COVERM_MIN_BREADTH}" \
        --output-file "${COVERM_OUT_DIR}/coverm_min75.txt" \
        2> "${COVERM_OUT_DIR}/min75_stats.txt"
fi

# ---- trimmed mean -----------------------------------------------------------
if [[ -f "${COVERM_OUT_DIR}/coverm_trimmed_mean.txt" ]]; then
    echo "  coverm_trimmed_mean.txt already exists — skipping."
else
    echo "  [CoverM] trimmed_mean..."
    coverm genome \
        --proper-pairs-only \
        --genome-fasta-extension fa \
        --genome-fasta-directory "${CROSS_RENAMED_DIR}" \
        --bam-files "${ALL_BAMS[@]}" \
        --threads "${BOWTIE2_THREADS}" \
        --min-read-percent-identity-pair "${COVERM_MIN_ID}" \
        -m trimmed_mean \
        --output-file "${COVERM_OUT_DIR}/coverm_trimmed_mean.txt" \
        2> "${COVERM_OUT_DIR}/trimmed_mean_stats.txt"
fi

echo ""
echo "CoverM tables in: ${COVERM_OUT_DIR}"
echo "=== [11_coverm] done ==="
