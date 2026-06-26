#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-trait_program_gene_model}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/trait_program_gene_model}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"

TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"
MANIFEST_OUT="${MANIFEST_OUT:-${META_DIR}/manifest_all.tsv}"

mkdir -p "${META_DIR}"

mapfile -t ALL_IDS < <(awk -F '\t' 'NR > 1 { sub(/\r$/, "", $2); print $2 }' "${FILE_ID_MAP}")

{
    printf 'figure_id\tfigure_kind\tsource_id\ttable_long_path\ttable_concordant_long_path\ttable_program_path\tgene_predictions_path\tprogram_rank_path\tregulator_coefficients_path\tpermutation_path\tplot_pdf\tplot_png\tstatus\n'
    for source_id in "${ALL_IDS[@]}"; do
        status="missing"
        [ -f "${STATUS_DIR}/${source_id}.ok" ] && status="ok"
        [ -f "${STATUS_DIR}/${source_id}.failed" ] && status="failed"
        [ -f "${STATUS_DIR}/${source_id}.running" ] && [ "${status}" = "missing" ] && status="running"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${source_id}" \
            "trait_program_gene_model" \
            "${source_id}" \
            "${TABLES_DIR}/${source_id}_long.tsv" \
            "${TABLES_DIR}/${source_id}_concordant_long.tsv" \
            "${TABLES_DIR}/${source_id}_programs.tsv" \
            "${TABLES_DIR}/${source_id}_gene_predictions.tsv" \
            "${TABLES_DIR}/${source_id}_program_rank.tsv" \
            "${TABLES_DIR}/${source_id}_regulator_coefficients.tsv" \
            "${TABLES_DIR}/${source_id}_permutation.tsv" \
            "${PLOTS_DIR}/${source_id}.pdf" \
            "${PLOTS_DIR}/${source_id}.png" \
            "${status}"
    done
} > "${MANIFEST_OUT}"

echo "Manifest written: ${MANIFEST_OUT}"
