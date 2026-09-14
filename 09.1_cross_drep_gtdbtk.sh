#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=200gb
#SBATCH --time=336:00:00
#SBATCH --job-name=cross_gtdbtk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_09.1_cross_drep_gtdbtk_%j.out
#SBATCH --error=slurm_09.1_cross_drep_gtdbtk_%j.err
# =============================================================================
# 09.1_cross_drep_gtdbtk.sh  —  GTDB-Tk taxonomy on cross-sample dereplicated MAGs.
#
# Runs AFTER 09_cross_drep.sh. Classifies the final non-redundant MAG catalog
# (the dereplicated_genomes/ directory from cross-sample dRep). This is the
# taxonomy that should be reported in publications, as it reflects the actual
# representative genomes in the MAG database rather than the pre-dRep set.
#
# NOTE: 07_gtdbtk.sh runs GTDB-Tk per-sample on MQ/HQ MAGs before dRep —
# that run is useful for exploration but should NOT be used for final reporting.
# This script is the authoritative taxonomy for the project.
#
# Skips if GTDB-Tk output directory already exists and contains summary files.
# Usage:  sbatch 09.1_cross_drep_gtdbtk.sh
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"
CROSS_GTDB_DIR="${CROSS_DREP_DIR}/gtdb_v2.7.0_r232"


# ---- Preflight checks -------------------------------------------------------
if [[ ! -d "${CROSS_DEREP_GENOMES}" ]]; then
    echo "ERROR: cross-sample dereplicated_genomes directory not found."
    echo "  Expected: ${CROSS_DEREP_GENOMES}"
    echo "  Run 09_cross_drep.sh first."
    exit 1
fi

N_MAGS=$(ls "${CROSS_DEREP_GENOMES}"/*.fa 2>/dev/null | wc -l)
if [[ "${N_MAGS}" -eq 0 ]]; then
    echo "ERROR: no .fa files found in ${CROSS_DEREP_GENOMES}."
    echo "  Run 09_cross_drep.sh first."
    exit 1
fi

# ---- Skip if already completed ----------------------------------------------
if ls "${CROSS_GTDB_DIR}"/gtdbtk.*.summary.tsv &>/dev/null 2>&1; then
    echo "GTDB-Tk output already exists — skipping."
    echo "  Output: ${CROSS_GTDB_DIR}"
    exit 0
fi

mkdir -p "${CROSS_GTDB_DIR}"

echo "=== [09.1_cross_drep_gtdbtk] ==="
echo "  Input:   ${CROSS_DEREP_GENOMES}"
echo "  Output:  ${CROSS_GTDB_DIR}"
echo "  MAGs:    ${N_MAGS}"
echo ""

source /home/opt/Miniconda3/miniconda3/bin/activate gtdbtk_v2.7.0
gtdbtk classify_wf \
    -x fa \
    --genome_dir "${CROSS_DEREP_GENOMES}" \
    --out_dir "${CROSS_GTDB_DIR}" \
    --cpus "${GTDBTK_CPUS}"
conda deactivate

echo ""
echo "  GTDB-Tk complete. Summary files:"
ls "${CROSS_GTDB_DIR}"/gtdbtk.*.summary.tsv 2>/dev/null || echo "  (no summary files found — check for errors above)"
echo "=== [09.1_cross_drep_gtdbtk] done ==="
