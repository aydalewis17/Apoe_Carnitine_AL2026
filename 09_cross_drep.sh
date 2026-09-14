#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=200gb
#SBATCH --time=48:00:00
#SBATCH --job-name=cross_drep
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aydalewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_09_cross_drep_%j.out
#SBATCH --error=slurm_09_cross_drep_%j.err
# =============================================================================
# 09_cross_drep.sh  —  Phase 3: Cross-sample dRep + rename + concatenate.
#
# PURPOSE
# -------
# Dereplicates all MQ/HQ MAGs from ALL assembly strategies into a single
# non-redundant MAG catalog at 99% ANI (strain-level dereplication).
#
# INPUT
# -----
# ALL_MQHQ_MAGs/ — the pooled directory produced by 02e_pool_mqhq_mags.sh.
# Contains 938 pre-dRep MQ/HQ MAGs from five sources:
#   (1) Per-sample assemblies        (460 MAGs, 10 samples)
#   (2) Control co-assembly          (83 MAGs)
#   (3) Carnitine co-assembly        (102 MAGs)
#   (4) Iterative round 1            (32 MAGs, control + carnitine)
#   (5) Iterative round 2            (26 MAGs, control + carnitine)
#
# This script must NOT be run on per-sample MAGs only. Always use the full
# ALL_MQHQ_MAGs/ pool so that co-assembly and iterative MAGs are included
# in dereplication. (Note: Phase 3 was run once on per-sample MAGs only
# as an interim step before Phase 2 was complete — those results are
# superseded by this run on the full pool.)
#
# GENOME_INFO.CSV
# ---------------
# dRep requires CheckM2 quality scores via --genomeInfo so it uses
# pre-computed completeness/contamination rather than re-running its own
# internal CheckM assessment. This script builds genome_info.csv from all
# five CheckM2 quality_report.tsv files before running dRep.
#
# The genome column must match the BARE FILENAME of each .fa file passed
# to dRep via -g (e.g. "con1_c1r_bin.3.fa", not the full path).
# This matches the naming convention used by 02e_pool_mqhq_mags.sh.
#
# DREP ANI
# --------
# dRep runs at 99% ANI (DREP_ANI="0.99" in 00_config.sh) for strain-level
# dereplication. This is more stringent than the dRep default (95%) and
# ensures that closely related strains are kept as separate representatives
# rather than collapsed. This follows the workshop recommendation.
#
# STEPS
# -----
#   1. Verify ALL_MQHQ_MAGs/ exists and is non-empty
#   2. Build genome_info.csv from all five CheckM2 sources
#   3. dRep dereplicate at 99% ANI with --genomeInfo
#   4. Rename contig headers with rename_bins_like_dram.py
#   5. Concatenate renamed MAGs into shared database FASTA
#
# WHAT RUNS AFTER THIS
# --------------------
#   09.1_cross_drep_gtdbtk.sh  — GTDB-Tk taxonomy on dereplicated representatives
#   10_bowtie2_map.sh           — map all 13 samples to shared MAG catalog
#   11_coverm.sh                — between-sample abundance tables
#   singleM_pipe.py             — re-run to complete Steps 2+3 (appraise)
#
# SKIP LOGIC
# ----------
# genome_info.csv, dRep, rename, and concatenation each checked independently.
#
# Usage: sbatch 09_cross_drep.sh
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

# Variables that must be in 00_config.sh (add if missing):
#   CATALOG_SAMPLES   — from acontrol_list.txt and carnitine_list.txt
#   DREP_ANI       — "0.99"
#   POOL_DIR       — "${PROJECT_DIR}/ALL_MQHQ_MAGs"
#   COASSEMBLY_DIR — "${PROJECT_DIR}/coassemblies"

CONTROL_LIST="${SCRIPTS_DIR}/control_list.txt"
CARNITINE_LIST="${SCRIPTS_DIR}/carnitine_list.txt"
for f in "${CONTROL_LIST}" "${CARNITINE_LIST}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required sample list not found: ${f}"
        exit 1
    fi
done
CATALOG_SAMPLES=$(awk 1 "${CONTROL_LIST}" "${CARNITINE_LIST}" | tr '\n' ' ')

[[ -z "${DREP_ANI:-}" ]] && {
    echo "ERROR: DREP_ANI not set in 00_config.sh. Add: export DREP_ANI=\"0.99\""
    exit 1
}
[[ -z "${POOL_DIR:-}" ]] && {
    echo "ERROR: POOL_DIR not set in 00_config.sh. Add: export POOL_DIR=\"\${PROJECT_DIR}/ALL_MQHQ_MAGs\""
    exit 1
}

mkdir -p "${CROSS_DREP_DIR}" "${CROSS_BOWTIE_DIR}"

# =============================================================================
# Step 1: Verify ALL_MQHQ_MAGs/ exists and is non-empty
# =============================================================================
echo "=== [09_cross_drep] Verifying MAG pool ==="
if [[ ! -d "${POOL_DIR}" ]]; then
    echo "ERROR: ALL_MQHQ_MAGs/ not found at ${POOL_DIR}."
    echo "       Run 02e_pool_mqhq_mags.sh first."
    exit 1
fi

mapfile -t ALL_MAG_FILES < <(find "${POOL_DIR}" -maxdepth 1 -name '*.fa' | sort)
N_TOTAL="${#ALL_MAG_FILES[@]}"

if [[ "${N_TOTAL}" -eq 0 ]]; then
    echo "ERROR: No .fa files found in ${POOL_DIR}."
    echo "       Run 02e_pool_mqhq_mags.sh first."
    exit 1
fi
echo "  Found ${N_TOTAL} MAGs in ${POOL_DIR}"

# =============================================================================
# Step 2: Build genome_info.csv from all five CheckM2 sources
#
# genome column = bare filename matching what's in POOL_DIR
# (e.g. "con1_c1r_bin.3.fa" for per-sample, "control_bin.5.fa" for co-assembly,
#  "control_iter1_bin.2.fa" for iterative)
#
# Naming convention set by 02e_pool_mqhq_mags.sh:
#   Per-sample:   {SAMPLE}_{binname}.fa       (e.g. con1_c1r_bin.3.fa)
#   Co-assembly:  {GROUP}_{binname}.fa        (e.g. control_bin.5.fa)
#   Iterative:    {GROUP}_iter{N}_{binname}.fa (e.g. control_iter1_bin.2.fa)
# =============================================================================
GENOME_INFO="${CROSS_DREP_DIR}/genome_info.csv"

if [[ -f "${GENOME_INFO}" ]]; then
    N_ROWS=$(tail -n +2 "${GENOME_INFO}" | wc -l)
    echo "=== [09_cross_drep] genome_info.csv already exists (${N_ROWS} MAGs) — skipping build. ==="
else
    echo "=== [09_cross_drep] Building genome_info.csv from CheckM2 results ==="
    echo "genome,completeness,contamination" > "${GENOME_INFO}"
    N_WRITTEN=0
    N_WARN=0

    # Source 1: Per-sample MAGs
    echo "  Processing per-sample CheckM2 results..."
    for SAMPLE in ${CATALOG_SAMPLES}; do
        TSV="${PROJECT_DIR}/${SAMPLE}/megahit_out/checkm2_v1.1.0_fa/quality_report.tsv"
        MAG_DIR="${PROJECT_DIR}/${SAMPLE}/MAGs"
        if [[ ! -f "${TSV}" ]]; then
            echo "  WARNING: no CheckM2 results for ${SAMPLE}"
            N_WARN=$((N_WARN + 1))
            continue
        fi
        while IFS=$'\t' read -r name comp cont _rest; do
            [[ "${name}" == "Name" ]] && continue
            BARE="${SAMPLE}_${name}.fa"
            if [[ -f "${POOL_DIR}/${BARE}" ]]; then
                echo "${BARE},${comp},${cont}" >> "${GENOME_INFO}"
                N_WRITTEN=$((N_WRITTEN + 1))
            fi
        done < "${TSV}"
    done

    # Source 2 & 3: Co-assembly MAGs (control and carnitine)
    echo "  Processing co-assembly CheckM2 results..."
    for GROUP in control carnitine; do
        TSV="${COASSEMBLY_DIR}/${GROUP}/checkm2_v1.1.0_fa/quality_report.tsv"
        if [[ ! -f "${TSV}" ]]; then
            echo "  WARNING: no co-assembly CheckM2 results for ${GROUP}"
            N_WARN=$((N_WARN + 1))
            continue
        fi
        while IFS=$'\t' read -r name comp cont _rest; do
            [[ "${name}" == "Name" ]] && continue
            BARE="${GROUP}_${name}.fa"
            if [[ -f "${POOL_DIR}/${BARE}" ]]; then
                echo "${BARE},${comp},${cont}" >> "${GENOME_INFO}"
                N_WRITTEN=$((N_WRITTEN + 1))
            fi
        done < "${TSV}"
    done

    # Source 4 & 5: Iterative MAGs (rounds 1 and 2, control and carnitine)
    echo "  Processing iterative assembly CheckM2 results..."
    for GROUP in control carnitine; do
        for ITER in 1 2; do
            TSV="${COASSEMBLY_DIR}/${GROUP}/iter${ITER}/checkm2_v1.1.0_fa/quality_report.tsv"
            if [[ ! -f "${TSV}" ]]; then
                echo "  WARNING: no iter${ITER} CheckM2 results for ${GROUP}"
                N_WARN=$((N_WARN + 1))
                continue
            fi
            while IFS=$'\t' read -r name comp cont _rest; do
                [[ "${name}" == "Name" ]] && continue
                BARE="${GROUP}_iter${ITER}_${name}.fa"
                if [[ -f "${POOL_DIR}/${BARE}" ]]; then
                    echo "${BARE},${comp},${cont}" >> "${GENOME_INFO}"
                    N_WRITTEN=$((N_WRITTEN + 1))
                fi
            done < "${TSV}"
        done
    done

    echo "  genome_info.csv complete: ${N_WRITTEN} MAGs written"
    [[ "${N_WARN}" -gt 0 ]] && echo "  WARNING: ${N_WARN} CheckM2 source(s) had missing results."

    # Sanity check — every MAG in the pool should have a genome_info entry
    N_CSV=$(tail -n +2 "${GENOME_INFO}" | wc -l)
    if [[ "${N_CSV}" -ne "${N_TOTAL}" ]]; then
        echo "  WARNING: genome_info.csv has ${N_CSV} entries but pool has ${N_TOTAL} MAGs."
        echo "           MAGs without CheckM2 entries will be assessed by dRep's internal CheckM."
    fi
fi

# =============================================================================
# Step 3: dRep at 99% ANI
# =============================================================================
if [[ -d "${CROSS_DEREP_GENOMES}" ]] && find "${CROSS_DEREP_GENOMES}" -name '*.fa' | grep -q .; then
    echo "=== [09_cross_drep] dereplicated_genomes already exists — skipping dRep. ==="
else
    echo "=== [09_cross_drep] Running dRep at ${DREP_ANI} ANI on ${N_TOTAL} MAGs ==="
    source /home/opt/Miniconda3/miniconda3/bin/activate drep_3.4.2
    dRep dereplicate "${CROSS_DREP_DIR}" \
        -p "${THREADS}" \
        -comp "${DREP_MIN_COMP}" \
        -con "${DREP_MAX_CONT}" \
        -sa "${DREP_ANI}" \
        --genomeInfo "${GENOME_INFO}" \
        -g "${ALL_MAG_FILES[@]}"
    conda deactivate

    N_DEREP=$(find "${CROSS_DEREP_GENOMES}" -name '*.fa' | wc -l)
    echo "  dRep complete: ${N_DEREP} representative MAGs in ${CROSS_DEREP_GENOMES}"
fi

# =============================================================================
# Step 4: Rename contig headers with rename_bins_like_dram.py
#
# Required for DRAM2 compatibility. Rewrites contig headers in each MAG
# to include the bin name, preventing header collisions when MAGs are
# concatenated into the shared database FASTA.
# =============================================================================
if [[ -d "${CROSS_RENAMED_DIR}" ]] && find "${CROSS_RENAMED_DIR}" -name '*.fa' | grep -q .; then
    echo "=== [09_cross_drep] Renamed genomes already exist — skipping rename. ==="
else
    echo "=== [09_cross_drep] Renaming contig headers... ==="
    [[ -d "${CROSS_RENAMED_DIR}" ]] && rm -rf "${CROSS_RENAMED_DIR}"
    set +u
    source /opt/Miniconda2/miniconda2/bin/activate scripts
    set -u
    python /ORG-Data/scripts/rename_bins_like_dram.py \
        -i "${CROSS_DEREP_GENOMES}/*.fa" \
        -o "${CROSS_RENAMED_DIR}"
    set +u
    conda deactivate
    set -u
    N_RENAMED=$(find "${CROSS_RENAMED_DIR}" -name '*.fa' | wc -l)
    echo "  Renamed ${N_RENAMED} MAGs -> ${CROSS_RENAMED_DIR}"
fi

# =============================================================================
# Step 5: Concatenate renamed MAGs into shared database FASTA
#
# This is the reference used by 10_bowtie2_map.sh for all-sample mapping
# and by singleM_pipe.py Step 2 for MAG-level community profiling.
# =============================================================================
if [[ -f "${CROSS_MAG_DB_FA}" ]]; then
    echo "=== [09_cross_drep] Shared MAG database already exists — skipping concatenation. ==="
else
    echo "=== [09_cross_drep] Concatenating renamed MAGs into shared database... ==="
    > "${CROSS_MAG_DB_FA}"
    while IFS= read -r -d '' f; do
        cat "${f}" >> "${CROSS_MAG_DB_FA}"
    done < <(find "${CROSS_RENAMED_DIR}" -maxdepth 1 -name '*.fa' -print0 | sort -z)
    echo "  Shared MAG database: ${CROSS_MAG_DB_FA}"
fi

echo ""
echo "========================================================="
echo " 09_cross_drep complete."
N_FINAL=$(find "${CROSS_DEREP_GENOMES}" -name '*.fa' | wc -l)
echo " Input MAGs:              ${N_TOTAL}"
echo " Dereplicated MAGs:       ${N_FINAL}"
echo " Renamed MAGs:            ${CROSS_RENAMED_DIR}"
echo " Shared database FASTA:   ${CROSS_MAG_DB_FA}"
echo ""
echo " NEXT STEPS:"
echo "   09.1_cross_drep_gtdbtk.sh  — taxonomy on dereplicated representatives"
echo "   10_bowtie2_map.sh           — map all samples to shared catalog"
echo "   11_coverm.sh                — between-sample abundance tables"
echo "   singleM_pipe.py             — re-run to complete appraise (Steps 2+3)"
echo "========================================================="
echo "=== [09_cross_drep] done ==="
