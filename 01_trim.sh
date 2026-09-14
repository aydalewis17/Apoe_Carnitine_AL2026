#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --job-name=trim
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_01_trim_%j.out
#SBATCH --error=slurm_01_trim_%j.err
# =============================================================================
# 01_trim.sh  —  Steps 1 + 1.2: sickle + bbduk trimming.
# Runs FastQC on post-sickle reads, then bbduk for adapter + poly-G removal.
# Skips if final trimmed reads already exist.
# Trimmed reads are left uncompressed.
# Usage:  sbatch 01_trim.sh SAMPLE [SAMPLE2 ...]
# =============================================================================
 
set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"
 
SAMPLES="${*:-${SAMPLES}}"
[[ -z "${SAMPLES}" ]] && { echo "ERROR: no samples specified."; exit 1; }
 
count_reads() {
    awk 'END{print NR/4}' "$1"
}
 
count_bp() {
    awk 'NR%4==2{bp+=length($0)} END{print bp}' "$1"
}
 
# Returns 0 if both R1 and R2 post-sickle FastQC HTML files exist
post_sickle_fastqc_done() {
    local fqc_dir="$1" prefix="$2"
    local r1_html r2_html
    r1_html=$(find "${fqc_dir}" -name "${prefix}_R1_sickle_trimmed_fastqc.html" 2>/dev/null | head -1)
    r2_html=$(find "${fqc_dir}" -name "${prefix}_R2_sickle_trimmed_fastqc.html" 2>/dev/null | head -1)
    [[ -n "${r1_html}" && -n "${r2_html}" ]]
}
 
for SAMPLE in ${SAMPLES}; do
    echo "=== [01_trim] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"
    mkdir -p "${TRIMMED_DIR}"
 
    if [[ -z "${R1_RAW}" || -z "${R2_RAW}" ]]; then
        echo "  ERROR: raw reads not found for ${SAMPLE}."
        continue
    fi
 
    # ---- Skip if trimmed reads already exist --------------------------------
    if [[ -f "${TRIMMED_R1}" && -f "${TRIMMED_R2}" ]]; then
        echo "  Trimmed reads already exist — skipping."
        echo "=== [01_trim] ${SAMPLE} done ==="
        continue
    fi
 
    # ---- Decompress if gzipped (sickle requires plain fastq) ----------------
    if [[ "${R1_RAW}" == *.gz ]]; then
        echo "  Decompressing raw reads..."
        R1_PLAIN="${RAW_DIR}/${SAMPLE}_R1_001.fastq"
        R2_PLAIN="${RAW_DIR}/${SAMPLE}_R2_001.fastq"
        gunzip -k "${R1_RAW}"
        gunzip -k "${R2_RAW}"
        DECOMP=true
    else
        R1_PLAIN="${R1_RAW}"
        R2_PLAIN="${R2_RAW}"
        DECOMP=false
    fi
 
    # ========== Step 1: sickle ===============================================
    # Only run if sickle outputs don't already exist (allows resuming after a
    # crash between sickle and bbduk)
    if [[ -f "${SICKLE_R1}" && -f "${SICKLE_R2}" ]]; then
        echo "  Sickle outputs already exist — skipping sickle."
    else
        echo "  [sickle] quality trimming..."
        sickle pe \
            -f "${R1_PLAIN}" -r "${R2_PLAIN}" \
            -t sanger \
            -o "${SICKLE_R1}" -p "${SICKLE_R2}" \
            -s "${DISCARDED}"
    fi
 
    if [[ "${DECOMP}" == true ]]; then
        rm -f "${R1_PLAIN}" "${R2_PLAIN}"
        echo "  Decompressed copies removed."
    fi
 
    # ========== FastQC on post-sickle reads ===================================
    FASTQC_DIR="${TRIMMED_DIR}/fastqc"
    mkdir -p "${FASTQC_DIR}"
 
    if post_sickle_fastqc_done "${FASTQC_DIR}" "${SAMPLE}"; then
        echo "  Post-sickle FastQC already done — skipping."
    else
        echo "  [FastQC] running on post-sickle reads..."
        fastqc --threads 2 --outdir "${FASTQC_DIR}" "${SICKLE_R1}" "${SICKLE_R2}"
        echo "  Post-sickle FastQC reports: ${FASTQC_DIR}/"
    fi
 
    # ========== Step 1.2: bbduk ==============================================
    echo "  [bbduk] adapter and poly-G removal..."
    bbduk.sh \
        threads=10 overwrite=t \
        in1="${SICKLE_R1}" in2="${SICKLE_R2}" \
        ref=/opt/bbtools/bbmap/resources/adapters.fa \
        tpe tbo \
        trimpolygright="${POLYG_TRIM}" \
        out1="${TRIMMED_R1}" out2="${TRIMMED_R2}"
 
    # ========== Read count sanity check =======================================
    TRIM_R1_READS=$(count_reads "${TRIMMED_R1}")
    TRIM_R2_READS=$(count_reads "${TRIMMED_R2}")
    if [[ "${TRIM_R1_READS}" -ne "${TRIM_R2_READS}" ]]; then
        echo "  WARNING: R1/R2 read count mismatch — R1=${TRIM_R1_READS}, R2=${TRIM_R2_READS}"
    fi
 
    # ========== Update metadata ==============================================
    echo "  Counting trimmed reads..."
    TRIM_READS="${TRIM_R1_READS}"
    R1_BP=$(count_bp "${TRIMMED_R1}")
    R2_BP=$(count_bp "${TRIMMED_R2}")
    TRIM_GBP=$(awk -v bp1="${R1_BP}" -v bp2="${R2_BP}" 'BEGIN{printf "%.4f", (bp1+bp2)/1e9}')
 
    SAMPLE_IN_META=$(awk -F'\t' -v s="${SAMPLE}" 'NR>1 && $1==s {print 1; exit}' "${METADATA}")
    if [[ -n "${SAMPLE_IN_META}" ]]; then
        python3 - "${METADATA}" "${SAMPLE}" "${TRIM_READS}" "${TRIM_GBP}" <<'PYEOF'
import sys, csv, os
tsv, sample, tr, tg = sys.argv[1:]
rows = []
with open(tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    fieldnames = reader.fieldnames
    for row in reader:
        if row['sample'] == sample:
            row['trimmed_reads']   = tr
            row['TRIMMED_GBP_SEQ'] = tg
        rows.append(row)
tmp = tsv + '.tmp'
with open(tmp, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
os.replace(tmp, tsv)
PYEOF
        echo "  Metadata updated."
    else
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${SAMPLE}" "NA" "NA" "${TRIM_READS}" "${TRIM_GBP}" \
            "NA" "NA" "NA" "NA" "NA" "NA" "NA" >> "${METADATA}"
        echo "  New metadata stub row written."
    fi
 
    echo "  trimmed reads: ${TRIM_READS}  |  trimmed Gbp: ${TRIM_GBP}"
    echo "=== [01_trim] ${SAMPLE} done ==="
done