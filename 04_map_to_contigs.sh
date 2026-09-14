#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=450gb
#SBATCH --time=336:00:00
#SBATCH --job-name=map_contigs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_04_map_to_contigs_%j.out
#SBATCH --error=slurm_04_map_to_contigs_%j.err
# =============================================================================
# 04_map_to_contigs.sh  —  Steps 4-6.1: BBMap -> SAM -> BAM -> sort -> filter.
# Skips if the filtered 99% BAM already exists.
# Usage:  sbatch 04_map_to_contigs.sh SAMPLE [SAMPLE2 ...]
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

SAMPLES="${*:-${SAMPLES}}"
[[ -z "${SAMPLES}" ]] && { echo "ERROR: no samples specified."; exit 1; }

for SAMPLE in ${SAMPLES}; do
    echo "=== [04_map_to_contigs] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"

    FILT_BAM="${MEGAHIT_DIR}/${SAMPLE}_B_mapped99per.sorted.bam"

    if [[ ! -f "${FILTERED_SCAFFOLDS}" ]]; then
        echo "  ERROR: filtered scaffolds not found. Run 03_filter_contigs.sh first."
        continue
    fi

    # ---- Skip if filtered BAM already exists --------------------------------
    if [[ -f "${FILT_BAM}" ]]; then
        echo "  Filtered BAM already exists — skipping."
        echo "=== [04_map_to_contigs] ${SAMPLE} done ==="
        continue
    fi

    cd "${MEGAHIT_DIR}"

    echo "  [BBMap] mapping reads to contigs..."
    bbmap.sh \
        -Xmx"${MEMORY}" threads="${THREADS}" \
        minid=90 overwrite=t \
        ref="${FILTERED_SCAFFOLDS}" \
        in1="${TRIMMED_R1}" in2="${TRIMMED_R2}" \
        out="${SAMPLE}_B_mapped.sam"

    echo "  [samtools view] SAM -> BAM..."
    samtools view -@ "${THREADS}" -bS "${SAMPLE}_B_mapped.sam" > "${SAMPLE}_B_mapped.bam"

    echo "  [samtools sort] sorting..."
    samtools sort -T "${SAMPLE}.sorted" -o "${SAMPLE}_B_mapped.sorted.bam" \
        "${SAMPLE}_B_mapped.bam" -@ "${THREADS}"

    echo "  [reformat.sh] identity filter >= ${BBMAP_MIN_ID}..."
    reformat.sh \
        -Xmx"${MEMORY}" \
        minidfilter="${BBMAP_MIN_ID}" \
        in="${SAMPLE}_B_mapped.sorted.bam" \
        out="${FILT_BAM}" \
        pairedonly=t primaryonly=t

    rm -f "${SAMPLE}_B_mapped.sam" "${SAMPLE}_B_mapped.bam"

    echo "  Filtered BAM: ${FILT_BAM}"
    echo "=== [04_map_to_contigs] ${SAMPLE} done ==="
done
