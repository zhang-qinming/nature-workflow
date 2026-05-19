#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-gwas_manhattan}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/gwas_manhattan}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"

TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"
MANIFEST_OUT="${MANIFEST_OUT:-${META_DIR}/manifest_all.tsv}"

mkdir -p "${META_DIR}"

mapfile -t ALL_IDS < <(awk -F '\t' 'NR > 1 { sub(/\r$/, "", $1); print $1 }' "${FILE_ID_MAP}")

{
    printf 'figure_id\tfigure_kind\tsource_id\ttable_path\tvariants_path\thits_path\tplot_pdf\tplot_png\tstatus\n'
    for source_id in "${ALL_IDS[@]}"; do
        variants_path="${TABLES_DIR}/${source_id}_variants.tsv"
        hits_path="${TABLES_DIR}/${source_id}_hits.tsv"
        plot_pdf="${PLOTS_DIR}/${source_id}.pdf"
        plot_png="${PLOTS_DIR}/${source_id}.png"
        status="missing"
        [ -f "${STATUS_DIR}/${source_id}.ok" ] && status="ok"
        [ -f "${STATUS_DIR}/${source_id}.failed" ] && status="failed"
        [ -f "${STATUS_DIR}/${source_id}.running" ] && [ "${status}" = "missing" ] && status="running"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${source_id}" \
            "gwas_manhattan" \
            "${source_id}" \
            "${variants_path}" \
            "${variants_path}" \
            "${hits_path}" \
            "${plot_pdf}" \
            "${plot_png}" \
            "${status}"
    done
} > "${MANIFEST_OUT}"

echo "Manifest written: ${MANIFEST_OUT}"
