#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=14-00:00:00
#SBATCH --mem=50gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=C831467393@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

# =============================================================================
# 00c_raw_qc.sh  —  Pre-flight: count raw reads and sequencing yield (Gbp)
#                   for each sample, run FastQC on raw reads, and write stats
#                   into the metadata file. Run BEFORE 01_trim.sh.
#
# Handles both plain .fastq and gzip-compressed .fastq.gz inputs.
# FastQC output goes to: {PROJECT_DIR}/fastqc_raw/{sample}/
#
# Skips read counting if the sample already has raw_reads filled in metadata.
# Skips FastQC if the output HTML files already exist for that sample.
#
# Usage:  bash 00c_raw_qc.sh SAMPLE [SAMPLE2 ...]
#   or    SAMPLES="S1 S2 S3" bash 00c_raw_qc.sh
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

SAMPLES="${*:-${SAMPLES}}"
[[ -z "${SAMPLES}" ]] && { echo "ERROR: no samples specified."; exit 1; }

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: metadata file not found at ${METADATA}."
    echo "Run 00b_init_metadata.sh first."
    exit 1
fi

# ---- Helper: count reads from plain or gzipped fastq (4 lines per record) --
count_reads_fastq() {
    local f="$1"
    if [[ "${f}" == *.gz ]]; then
        zcat "${f}" | awk 'END{print NR/4}'
    else
        awk 'END{print NR/4}' "${f}"
    fi
}

# ---- Helper: sum base pairs across sequence lines of a fastq ----------------
count_bp_fastq() {
    local f="$1"
    if [[ "${f}" == *.gz ]]; then
        zcat "${f}" | awk 'NR%4==2 {bp += length($0)} END {print bp}'
    else
        awk 'NR%4==2 {bp += length($0)} END {print bp}' "${f}"
    fi
}

# ---- Helper: update a single metadata column for a sample row ---------------
update_meta() {
    local sample="$1" col="$2" val="$3"
    python3 - "${METADATA}" "${sample}" "${col}" "${val}" <<'EOF'
import sys, csv, os
tsv, sample, col, val = sys.argv[1:]
rows = []
with open(tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = reader.fieldnames
    for row in reader:
        if row['sample'] == sample:
            row[col] = val
        rows.append(row)
tmp = tsv + '.tmp'
with open(tmp, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
os.replace(tmp, tsv)
EOF
}

# ---- Helper: check if raw_reads is already filled in metadata for a sample --
reads_already_counted() {
    local sample="$1"
    local val
    val=$(awk -F'\t' -v s="${sample}" 'NR>1 && $1==s {print $2; exit}' "${METADATA}")
    # Returns 0 (true) if the value exists and is not NA or empty
    [[ -n "${val}" && "${val}" != "NA" ]]
}

# ---- Helper: check if FastQC output already exists for a sample -------------
fastqc_already_done() {
    local sample="$1"
    local fastqc_dir="${PROJECT_DIR}/fastqc_raw/${sample}"
    # FastQC produces one _fastqc.html per input file — check both R1 and R2
    local r1_stem r2_stem
    r1_stem=$(basename "${R1_RAW}" | sed 's/\.fastq\.gz$//; s/\.fastq$//')
    r2_stem=$(basename "${R2_RAW}" | sed 's/\.fastq\.gz$//; s/\.fastq$//')
    [[ -f "${fastqc_dir}/${r1_stem}_fastqc.html" && \
       -f "${fastqc_dir}/${r2_stem}_fastqc.html" ]]
}

# =============================================================================
# Per-sample loop
# =============================================================================
for SAMPLE in ${SAMPLES}; do
    echo "=== [00c_raw_qc] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"

    # ---- Guard: raw reads must exist ----------------------------------------
    if [[ -z "${R1_RAW}" || -z "${R2_RAW}" ]]; then
        echo "  ERROR: raw reads not found for ${SAMPLE}. Looked for:"
        echo "    ${RAW_DIR}/${SAMPLE}_R1_001.fastq[.gz]"
        echo "    ${RAW_DIR}/${SAMPLE}_R2_001.fastq[.gz]"
        continue
    fi
    echo "  R1: ${R1_RAW}"
    echo "  R2: ${R2_RAW}"

    # ========== Count reads and base pairs ===================================
    if reads_already_counted "${SAMPLE}"; then
        echo "  Read counts already in metadata — skipping."
    else
        echo "  Counting raw reads and base pairs (this may take a few minutes)..."

        R1_READS=$(count_reads_fastq "${R1_RAW}")
        R1_BP=$(count_bp_fastq "${R1_RAW}")
        R2_BP=$(count_bp_fastq "${R2_RAW}")

        TOTAL_PAIRS="${R1_READS}"
        RAW_GBP=$(awk -v bp1="${R1_BP}" -v bp2="${R2_BP}" \
                  'BEGIN{printf "%.4f", (bp1 + bp2) / 1e9}')

        echo "  Raw read pairs:  ${TOTAL_PAIRS}"
        echo "  Raw yield (Gbp): ${RAW_GBP}"

        # Write stub row or update existing
        SAMPLE_IN_META=$(awk -F'\t' -v s="${SAMPLE}" 'NR>1 && $1==s {print 1; exit}' "${METADATA}")
        if [[ -z "${SAMPLE_IN_META}" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${SAMPLE}" \
                "${TOTAL_PAIRS}" \
                "${RAW_GBP}" \
                "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" \
                >> "${METADATA}"
            echo "  New row written to metadata."
        else
            update_meta "${SAMPLE}" "raw_reads"   "${TOTAL_PAIRS}"
            update_meta "${SAMPLE}" "RAW_GBP_SEQ" "${RAW_GBP}"
            echo "  Existing metadata row updated."
        fi
    fi

    # ========== FastQC ========================================================
    FASTQC_DIR="${PROJECT_DIR}/fastqc_raw/${SAMPLE}"
    mkdir -p "${FASTQC_DIR}"

    if fastqc_already_done "${SAMPLE}"; then
        echo "  FastQC output already exists — skipping."
    else
        echo "  [FastQC] running on R1 and R2..."
        fastqc \
            --outdir "${FASTQC_DIR}" \
            --threads 2 \
            "${R1_RAW}" \
            "${R2_RAW}"
        echo "  FastQC reports: ${FASTQC_DIR}/"
    fi

    echo "=== [00c_raw_qc] ${SAMPLE} done ==="
    echo ""
done

# =============================================================================
# MultiQC summary across all samples (runs once after per-sample FastQC)
# =============================================================================
MULTIQC_DIR="${PROJECT_DIR}/fastqc_raw/multiqc_raw"
echo "=== [00c_raw_qc] Running MultiQC across all samples ==="
mkdir -p "${MULTIQC_DIR}"

multiqc \
    "${PROJECT_DIR}/fastqc_raw" \
    --outdir "${MULTIQC_DIR}" \
    --filename "multiqc_raw_reads" \
    --force

echo "  MultiQC report: ${MULTIQC_DIR}/multiqc_raw_reads.html"
echo "=== [00c_raw_qc] all done ==="
