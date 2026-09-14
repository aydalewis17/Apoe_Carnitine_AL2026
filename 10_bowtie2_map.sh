#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=450gb
#SBATCH --time=336:00:00
#SBATCH --job-name=bowtie2_map
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_10_bowtie2_map_%j.out
#SBATCH --error=slurm_10_bowtie2_map_%j.err
# =============================================================================
# 10_bowtie2_map.sh  —  Step 16: Build shared Bowtie2 index, map all samples.
# Skips index build if index files already exist.
# Skips per-sample mapping if POSSORT BAM already exists.
# Usage:  sbatch 10_bowtie2_map.sh   (reads carnitine_list.txt and control_list.txt via 00_config.sh)
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

[[ -z "${CATALOG_SAMPLES}" ]] && { echo "ERROR: CATALOG_SAMPLES is empty. Check that carnitine_list.txt and control_list.txt exists in ${SCRIPTS_DIR}."; exit 1; }

if [[ ! -f "${CROSS_MAG_DB_FA}" ]]; then
    echo "ERROR: cross-sample database FASTA not found. Run 09_cross_drep.sh first."
    exit 1
fi

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

# ========== Build shared Bowtie2 index =======================================
if ls "${CROSS_BOWTIE_PREFIX}".*.bt2 &>/dev/null 2>&1; then
    echo "[bowtie2-build] index already exists — skipping."
else
    echo "[bowtie2-build] building shared index..."
    bowtie2-build "${CROSS_MAG_DB_FA}" "${CROSS_BOWTIE_PREFIX}" --threads "${BOWTIE2_THREADS}"
    echo "  Index built: ${CROSS_BOWTIE_PREFIX}"
fi

# ---- Init mapping summary ---------------------------------------------------
printf '%s\t%s\t%s\t%s\n' \
    "sample" "total_trimmed_pairs" "mapped_pairs" "percent_pairs_mapped" \
    > "${CON_CARN_MAPPING_SUMMARY}"

# ========== Map each sample ==================================================
ALL_POSSORT_BAMS=()

for SAMPLE in ${CATALOG_SAMPLES}; do
    echo "=== [10_bowtie2_map] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"
    mkdir -p "${MAP_DIR}"

    if [[ ! -f "${TRIMMED_R1}" || ! -f "${TRIMMED_R2}" ]]; then
        echo "  WARNING: trimmed reads not found for ${SAMPLE}, skipping."
        continue
    fi

    POSSORT_BAM="${MAP_DIR}/${SAMPLE}_mapped_crossDB_97id_POSSORT.bam"
    SAM_FILE="${MAP_DIR}/${SAMPLE}_mapped_crossDB.sam"
    BAM_FILE="${MAP_DIR}/${SAMPLE}_mapped_crossDB.bam"
    FILT_BAM="${MAP_DIR}/${SAMPLE}_mapped_crossDB_97id.bam"

    # ---- Skip if POSSORT BAM already exists ---------------------------------
    if [[ -f "${POSSORT_BAM}" ]]; then
        echo "  POSSORT BAM already exists — skipping mapping."
        # Still collect for CoverM and update summary
        TOTAL_PAIRS=$(awk 'END{print NR/4}' "${TRIMMED_R1}")
        MAPPED_READS=$(samtools view -c -F 4 "${POSSORT_BAM}")
        MAPPED_PAIRS=$(( MAPPED_READS / 2 ))
        PCT=$(awk -v mp="${MAPPED_PAIRS}" -v tp="${TOTAL_PAIRS}" \
              'BEGIN{printf "%.2f", (tp>0 ? mp/tp*100 : 0)}')
        printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${TOTAL_PAIRS}" "${MAPPED_PAIRS}" "${PCT}" >> "${CON_CARN_MAPPING_SUMMARY}"
        update_meta "${SAMPLE}" "percent_reads_mapped" "${PCT}"
        ALL_POSSORT_BAMS+=("${POSSORT_BAM}")
        echo "  ${SAMPLE}: ${PCT}% pairs mapped (from existing BAM)"
        echo "=== [10_bowtie2_map] ${SAMPLE} done ==="
        continue
    fi

    echo "  [bowtie2] mapping reads..."
    bowtie2 \
        -D 10 -R 2 -N 0 -L 22 -i S,0,2.50 \
        -p "${BOWTIE2_THREADS}" \
        -x "${CROSS_BOWTIE_PREFIX}" \
        -S "${SAM_FILE}" \
        -1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}"

    echo "  [samtools view] SAM -> BAM..."
    samtools view -bS "${SAM_FILE}" > "${BAM_FILE}"

    echo "  [reformat.sh] identity filter >= ${BOWTIE2_MIN_ID}..."
    reformat.sh \
        -Xmx"${MEMORY}" \
        minidfilter="${BOWTIE2_MIN_ID}" \
        in="${BAM_FILE}" out="${FILT_BAM}" \
        pairedonly=t primaryonly=t

    echo "  [samtools sort] position sorting..."
    samtools sort -@ "${THREADS}" -o "${POSSORT_BAM}" "${FILT_BAM}"

    rm -f "${SAM_FILE}" "${BAM_FILE}" "${FILT_BAM}"

    TOTAL_PAIRS=$(awk 'END{print NR/4}' "${TRIMMED_R1}")
    MAPPED_READS=$(samtools view -c -F 4 "${POSSORT_BAM}")
    MAPPED_PAIRS=$(( MAPPED_READS / 2 ))
    PCT=$(awk -v mp="${MAPPED_PAIRS}" -v tp="${TOTAL_PAIRS}" \
          'BEGIN{printf "%.2f", (tp>0 ? mp/tp*100 : 0)}')

    printf '%s\t%s\t%s\t%s\n' "${SAMPLE}" "${TOTAL_PAIRS}" "${MAPPED_PAIRS}" "${PCT}" >> "${CON_CARN_MAPPING_SUMMARY}"
    update_meta "${SAMPLE}" "percent_reads_mapped" "${PCT}"

    ALL_POSSORT_BAMS+=("${POSSORT_BAM}")
    echo "  ${SAMPLE}: ${PCT}% pairs mapped"
    echo "=== [10_bowtie2_map] ${SAMPLE} done ==="
done

BAM_LIST="${PROJECT_DIR}/possort_bam_list.txt"
printf '%s\n' "${ALL_POSSORT_BAMS[@]}" > "${BAM_LIST}"
echo ""
echo "POSSORT BAM list saved to: ${BAM_LIST}"
echo "Mapping summary saved to:  ${CON_CARN_MAPPING_SUMMARY}"
echo "=== [10_bowtie2_map] all samples done ==="
