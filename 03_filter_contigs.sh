#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=50gb
#SBATCH --time=12:00:00
#SBATCH --job-name=filter_contigs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ayda.lewis@colostate.edu
#SBATCH --partition=borton-hi,borton-low
#SBATCH --output=slurm_03_filter_contigs_%j.out
#SBATCH --error=slurm_03_filter_contigs_%j.err
# =============================================================================
# 03_filter_contigs.sh  —  Steps 3 + 8: pullseq length filter + contig stats.
# Skips if filtered scaffolds already exist.
# Usage:  sbatch 03_filter_contigs.sh SAMPLE [SAMPLE2 ...]
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
    echo "=== [03_filter_contigs] ${SAMPLE} ==="
    set_sample_dirs "${SAMPLE}"

    RAW_CONTIGS="${MEGAHIT_DIR}/final.contigs.fa"

    if [[ ! -f "${RAW_CONTIGS}" ]]; then
        echo "  ERROR: assembly not found. Run 02_assemble.sh first."
        continue
    fi

    # ---- Skip if filtered scaffolds already exist ---------------------------
    if [[ -f "${FILTERED_SCAFFOLDS}" ]]; then
        echo "  Filtered scaffolds already exist — skipping."
        echo "=== [03_filter_contigs] ${SAMPLE} done ==="
        continue
    fi

    # ========== pullseq ======================================================
    echo "  [pullseq] filtering contigs >= ${MIN_CONTIG_LEN} bp..."
    pullseq -i "${RAW_CONTIGS}" -m "${MIN_CONTIG_LEN}" > "${FILTERED_SCAFFOLDS}"

    # ========== contig_stats.pl ==============================================
    echo "  [contig_stats.pl] computing assembly stats..."
    STATS_OUT="${MEGAHIT_DIR}/${SAMPLE}_final.contigs_STATS"
    "${CONTIG_STATS}" -i "${RAW_CONTIGS}" -o "${STATS_OUT}"

    STATS_FILT="${MEGAHIT_DIR}/${SAMPLE}_filtered_contigs_STATS"
    "${CONTIG_STATS}" -i "${FILTERED_SCAFFOLDS}" -o "${STATS_FILT}"

    # ========== Parse and update metadata ====================================
    N_CONTIGS_ALL=$(grep -c '^>' "${RAW_CONTIGS}")
    N_CONTIGS_FILT=$(grep -c '^>' "${FILTERED_SCAFFOLDS}")

    N50=$(grep 'N50:' "${STATS_FILT}.summary.txt" | awk '{print $2}')
    LONGEST=$(grep '^1\.' "${STATS_FILT}.summary.txt" | awk '{print $2}')

    if [[ -z "${N50}" ]]; then
        N50=$(awk '/^>/{if(seq)print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq)print length(seq)}' \
              "${FILTERED_SCAFFOLDS}" | sort -rn | \
              awk '{sum+=$1; lens[NR]=$1} END{half=sum/2; t=0; for(i=1;i<=NR;i++){t+=lens[i]; if(t>=half){print lens[i]; exit}}}')
    fi
    if [[ -z "${LONGEST}" ]]; then
        LONGEST=$(awk '/^>/{if(seq)print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq)print length(seq)}' \
                  "${FILTERED_SCAFFOLDS}" | sort -rn | head -1)
    fi

    update_meta "${SAMPLE}" "CONTIGS_ASSEMBLY" "${N_CONTIGS_ALL}"
    update_meta "${SAMPLE}" "CONTIGS_GT2500"   "${N_CONTIGS_FILT}"
    update_meta "${SAMPLE}" "N50"              "${N50}"
    update_meta "${SAMPLE}" "longest_contig"   "${LONGEST}"

    echo "  contigs (all): ${N_CONTIGS_ALL}  |  >2500bp: ${N_CONTIGS_FILT}  |  N50: ${N50}  |  longest: ${LONGEST}"
    echo "  Metadata updated."
    echo "=== [03_filter_contigs] ${SAMPLE} done ==="
done
