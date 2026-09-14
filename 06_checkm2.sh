#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --job-name=checkm2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_06_checkm2_%j.out
#SBATCH --error=slurm_06_checkm2_%j.err
# =============================================================================
# 06_checkm2.sh  —  Steps 9-10: CheckM2 QC + filter/copy MQ/HQ MAGs.
# Skips CheckM2 if quality_report.tsv already exists.
# Skips MAG copy if MAGs/ already contains .fa files.
# Usage:  sbatch 06_checkm2.sh SAMPLE [SAMPLE2 ...]
# =============================================================================

set -euo pipefail
source "/home/projects-phoenix/ApoE_Carnitine/MetaG/scripts/00_config.sh"

SAMPLES="${*:-${SAMPLES}}"
[[ -z "${SAMPLES}" ]] && { echo "ERROR: no samples specified."; exit 1; }

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

for SAMPLE in ${SAMPLES}; do
    echo "=== [06_checkm2] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"
    mkdir -p "${MAGS_DIR}"

    if [[ ! -d "${BINS_DIR}" ]]; then
        echo "  ERROR: bins directory not found. Run 05_metabat.sh first."
        continue
    fi

    N_BINS=$(ls "${BINS_DIR}"/*.fa 2>/dev/null | wc -l)
    echo "  Total bins from MetaBAT: ${N_BINS}"
    update_meta "${SAMPLE}" "bins" "${N_BINS}"

    # ========== Step 9: CheckM2 ==============================================
    if [[ -f "${CHECKM2_RESULTS}" ]]; then
        echo "  CheckM2 results already exist — skipping CheckM2."
    else
        echo "  [CheckM2] assessing bin quality..."
        source /home/opt/Miniconda3/miniconda3/bin/activate checkm2_v1.1.0
        checkm2 predict \
            -x fa \
            --input "${BINS_DIR}" \
            --output-directory "${CHECKM2_OUT}" \
            --threads "${CHECKM2_THREADS}"
        conda deactivate
    fi

    # ========== Step 10: filter and copy MQ/HQ MAGs ==========================
    if ls "${MAGS_DIR}"/*.fa &>/dev/null; then
        N_MQHQ=$(ls "${MAGS_DIR}"/*.fa | wc -l)
        echo "  MQ/HQ MAGs already copied (${N_MQHQ} MAGs) — skipping copy."
        update_meta "${SAMPLE}" "MQHQ_bins" "${N_MQHQ}"
    else
        echo "  Filtering for completeness >= ${CHECKM2_MIN_COMP}% and contamination <= ${CHECKM2_MAX_CONT}%..."
        N_MQHQ=0
        while IFS=$'\t' read -r name comp cont _rest; do
            [[ "${name}" == "Name" ]] && continue
            passes=$(awk -v c="${comp}" -v co="${cont}" -v mc="${CHECKM2_MIN_COMP}" -v mx="${CHECKM2_MAX_CONT}" \
                     'BEGIN{print (c+0 >= mc+0 && co+0 <= mx+0) ? "yes" : "no"}')
            if [[ "${passes}" == "yes" ]]; then
                src="${BINS_DIR}/${name}.fa"
                dst="${MAGS_DIR}/${SAMPLE}_${name}.fa"
                [[ ! -f "${src}" ]] && src="${BINS_DIR}/${name}"
                if [[ -f "${src}" ]]; then
                    cp -p "${src}" "${dst}"
                    N_MQHQ=$((N_MQHQ + 1))
                fi
            fi
        done < "${CHECKM2_RESULTS}"

        [[ "${N_MQHQ}" -eq 0 ]] && echo "  WARNING: no MQ/HQ MAGs passed for ${SAMPLE}."
        update_meta "${SAMPLE}" "MQHQ_bins" "${N_MQHQ}"
        echo "  MQ/HQ MAGs: ${N_MQHQ} -> ${MAGS_DIR}"
    fi

    echo "  Metadata updated."
    echo "=== [06_checkm2] ${SAMPLE} done ==="
done
