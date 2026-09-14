#!/usr/bin/env bash
# =============================================================================
# 00b_init_metadata.sh
# Initialize the pipeline metadata TSV with a header row.
# Run ONCE at the start of a new project before any per-sample steps.
# =============================================================================

source "$(dirname "$0")/00_config.sh"

mkdir -p "${PROJECT_DIR}"

if [[ -f "${METADATA}" ]]; then
    echo "Metadata file already exists at: ${METADATA}"
    echo "Delete it manually if you want to start fresh."
    exit 0
fi

printf '%s\t' \
    "sample" \
    "raw_reads" \
    "RAW_GBP_SEQ" \
    "trimmed_reads" \
    "TRIMMED_GBP_SEQ" \
    "CONTIGS_ASSEMBLY" \
    "CONTIGS_GT2500" \
    "N50" \
    "longest_contig" \
    "bins" \
    "MQHQ_bins" \
    "percent_reads_mapped" \
    > "${METADATA}"
printf '\n' >> "${METADATA}"

echo "Metadata file initialized: ${METADATA}"
