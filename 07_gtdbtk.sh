#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=200gb
#SBATCH --time=336:00:00
#SBATCH --job-name=gtdbtk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_07_gtdbtk_%j.out
#SBATCH --error=slurm_07_gtdbtk_%j.err
# =============================================================================
# 07_gtdbtk.sh  —  Step 11: GTDB-Tk taxonomy assignment.
# Skips if GTDB-Tk output directory already exists and contains summary files.
# Usage:  sbatch 07_gtdbtk.sh SAMPLE [SAMPLE2 ...]
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

SAMPLES="${*:-${SAMPLES}}"
[[ -z "${SAMPLES}" ]] && { echo "ERROR: no samples specified."; exit 1; }

for SAMPLE in ${SAMPLES}; do
    echo "=== [07_gtdbtk] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"

    N_MAGS=$(ls "${MAGS_DIR}"/*.fa 2>/dev/null | wc -l)
    if [[ "${N_MAGS}" -eq 0 ]]; then
        echo "  No MAGs found in ${MAGS_DIR} — skipping."
        continue
    fi

    # ---- Skip if GTDB-Tk already completed ----------------------------------
    # gtdbtk always writes at least one summary TSV on completion
    if ls "${GTDB_DIR}"/gtdbtk.*.summary.tsv &>/dev/null 2>&1; then
        echo "  GTDB-Tk output already exists — skipping."
        echo "=== [07_gtdbtk] ${SAMPLE} done ==="
        continue
    fi

    mkdir -p "${GTDB_DIR}"
    echo "  [GTDB-Tk] classifying ${N_MAGS} MAGs..."

    source /home/opt/Miniconda3/miniconda3/bin/activate gtdbtk_v2.7.0
    gtdbtk classify_wf \
        -x fa \
        --genome_dir "${MAGS_DIR}" \
        --out_dir "${GTDB_DIR}" \
        --cpus "${GTDBTK_CPUS}"
    conda deactivate

    echo "  GTDB-Tk output: ${GTDB_DIR}"
    echo "=== [07_gtdbtk] ${SAMPLE} done ==="
done
