#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=200gb
#SBATCH --time=336:00:00
#SBATCH --job-name=pre_drep_gtdbtk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_08b_pre_drep_gtdbtk_%j.out
#SBATCH --error=slurm_08b_pre_drep_gtdbtk_%j.err
# =============================================================================
# 08b_pre_drep_gtdbtk.sh  —  GTDB-Tk taxonomy on all pre-dRep MQ/HQ MAGs.
#
# Runs AFTER 02e_pool_mqhq_mags.sh and BEFORE 09_cross_drep.sh.
# Classifies all 703 carnitine and control MAGs in ALL_MQHQ_MAGs/ to provide taxonomy for the
# full pre-dRep pool. Useful for exploration, provenance tracking, and
# confirming that organisms of interest (e.g. g__Eubacterium_E) are present
# before dereplication potentially collapses representatives.
#
# NOTE: This is NOT the authoritative taxonomy for final reporting.
#       Use 09.1_cross_drep_gtdbtk.sh output (post-dRep representatives)
#       for publication figures.
#
# Skips if output directory already exists and contains summary files.
# Usage:  sbatch 08b_pre_drep_gtdbtk.sh
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

PREDREP_GTDB_DIR="${PROJECT_DIR}/ALL_MQHQ_MAGs_gtdb_v2.7.0_r232"

# ---- Preflight checks -------------------------------------------------------
if [[ ! -d "${POOL_DIR}" ]]; then
    echo "ERROR: ALL_MQHQ_MAGs/ not found at ${POOL_DIR}."
    echo "  Run 02e_pool_mqhq_mags.sh first."
    exit 1
fi

N_MAGS=$(find "${POOL_DIR}" -maxdepth 1 -name '*.fa' | wc -l)
if [[ "${N_MAGS}" -eq 0 ]]; then
    echo "ERROR: No .fa files found in ${POOL_DIR}."
    echo "  Run 02e_pool_mqhq_mags.sh first."
    exit 1
fi

# ---- Skip if already completed ----------------------------------------------
if ls "${PREDREP_GTDB_DIR}"/gtdbtk.*.summary.tsv &>/dev/null 2>&1; then
    echo "GTDB-Tk output already exists — skipping."
    echo "  Output: ${PREDREP_GTDB_DIR}"
    exit 0
fi

mkdir -p "${PREDREP_GTDB_DIR}"

echo "=== [08b_pre_drep_gtdbtk] ==="
echo "  Input:   ${POOL_DIR}"
echo "  Output:  ${PREDREP_GTDB_DIR}"
echo "  MAGs:    ${N_MAGS}"
echo ""

set +u
source /home/opt/Miniconda3/miniconda3/bin/activate gtdbtk_v2.7.0
set -u

gtdbtk classify_wf \
    -x fa \
    --genome_dir "${POOL_DIR}" \
    --out_dir "${PREDREP_GTDB_DIR}" \
    --cpus "${GTDBTK_CPUS}"

set +u
conda deactivate
set -u

echo ""
echo "  GTDB-Tk complete. Summary files:"
ls "${PREDREP_GTDB_DIR}"/gtdbtk.*.summary.tsv 2>/dev/null || echo "  (no summary files found — check for errors above)"
echo "=== [08b_pre_drep_gtdbtk] done ==="
