#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=450gb
#SBATCH --time=336:00:00
#SBATCH --job-name=assemble
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_02_assemble_%j.out
#SBATCH --error=slurm_02_assemble_%j.err
# =============================================================================
# 02_assemble.sh  —  Step 2: MEGAHIT assembly.
# Skips if final.contigs.fa already exists.
# Usage:  sbatch 02_assemble.sh SAMPLE [SAMPLE2 ...]
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

SAMPLES="${*:-${SAMPLES}}"
[[ -z "${SAMPLES}" ]] && { echo "ERROR: no samples specified."; exit 1; }

for SAMPLE in ${SAMPLES}; do
    echo "=== [02_assemble] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"

    if [[ ! -f "${TRIMMED_R1}" || ! -f "${TRIMMED_R2}" ]]; then
        echo "  ERROR: trimmed reads not found. Run 01_trim.sh first."
        continue
    fi

    # ---- Skip if assembly already complete ----------------------------------
    if [[ -f "${MEGAHIT_DIR}/final.contigs.fa" ]]; then
        echo "  Assembly already exists — skipping."
        echo "=== [02_assemble] ${SAMPLE} done ==="
        continue
    fi

    # MEGAHIT will fail if output dir already exists (even partially)
    if [[ -d "${MEGAHIT_DIR}" ]]; then
        echo "  Removing incomplete MEGAHIT directory from previous run..."
        rm -rf "${MEGAHIT_DIR}"
    fi

    echo "  [megahit] assembling..."
    megahit \
        -1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}" \
        --k-min "${KMER_MIN}" --k-max "${KMER_MAX}" --k-step "${KMER_STEP}" \
        -m 0.4 -t "${THREADS}" \
        -o "${MEGAHIT_DIR}"

    echo "  Assembly: ${MEGAHIT_DIR}/final.contigs.fa"
    echo "=== [02_assemble] ${SAMPLE} done ==="
done
