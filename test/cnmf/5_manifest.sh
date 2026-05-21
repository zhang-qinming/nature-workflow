#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cnmf}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cnmf}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"

TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"
MANIFEST_OUT="${MANIFEST_OUT:-${META_DIR}/manifest_all.tsv}"

mkdir -p "${META_DIR}"

mapfile -t ALL_IDS < <(awk -F '\t' 'NR > 1 { sub(/\r$/, "", $2); print $2 }' "${FILE_ID_MAP}")

{
    printf 'figure_id\tfigure_kind\tsource_id\ttable_path\tplot_pdf\tplot_png\tstatus\n'
    for source_id in "${ALL_IDS[@]}"; do
        table_path="${TABLES_DIR}/program_regulator/${source_id}.tsv"
        plot_pdf="${PLOTS_DIR}/program_regulator/${source_id}.pdf"
        plot_png=""
        status="missing"
        [ -f "${STATUS_DIR}/${source_id}.ok" ] && status="ok"
        [ -f "${STATUS_DIR}/${source_id}.failed" ] && status="failed"
        [ -f "${STATUS_DIR}/${source_id}.running" ] && [ "${status}" = "missing" ] && status="running"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${source_id}" \
            "cnmf_program_regulator_scatter" \
            "${source_id}" \
            "${table_path}" \
            "${plot_pdf}" \
            "${plot_png}" \
            "${status}"
    done
} > "${MANIFEST_OUT}"

echo "Manifest written: ${MANIFEST_OUT}"
