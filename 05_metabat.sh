#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=200gb
#SBATCH --time=48:00:00
#SBATCH --job-name=metabat
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_05_metabat_%j.out
#SBATCH --error=slurm_05_metabat_%j.err
# =============================================================================
# 05_metabat.sh  —  Step 7: MetaBAT binning.
# Skips if bins directory already exists and contains .fa files.
# Usage:  sbatch 05_metabat.sh SAMPLE [SAMPLE2 ...]
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

SAMPLES="${*:-${SAMPLES}}"
[[ -z "${SAMPLES}" ]] && { echo "ERROR: no samples specified."; exit 1; }

for SAMPLE in ${SAMPLES}; do
    echo "=== [05_metabat] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"

    FILT_BAM="${MEGAHIT_DIR}/${SAMPLE}_B_mapped.sorted.bam" #AL Note: using unfiltered BAM file here, filtering unnecessary

    if [[ ! -f "${FILT_BAM}" ]]; then
        echo "  ERROR: filtered BAM not found. Run 04_map_to_contigs.sh first."
        continue
    fi

    # ---- Skip if bins already exist -----------------------------------------
    if [[ -d "${BINS_DIR}" ]] && ls "${BINS_DIR}"/*.fa &>/dev/null; then
        N_BINS=$(ls "${BINS_DIR}"/*.fa | wc -l)
        echo "  Bins already exist (${N_BINS} bins) — skipping."
        echo "=== [05_metabat] ${SAMPLE} done ==="
        continue
    fi

    cd "${MEGAHIT_DIR}"
    echo "  [MetaBAT] binning contigs..."
    runMetaBat.sh "${FILTERED_SCAFFOLDS}" "${FILT_BAM}"

    N_BINS=$(ls "${BINS_DIR}"/*.fa 2>/dev/null | wc -l)
    echo "  Bins produced: ${N_BINS}"
    echo "  Bin directory: ${BINS_DIR}"
    echo "=== [05_metabat] ${SAMPLE} done ==="
done
