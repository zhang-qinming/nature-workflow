#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TASK_NAME="${TASK_NAME:-cross_trait_scatter}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cross_trait}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"

TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"
MANIFEST_OUT="${MANIFEST_OUT:-${META_DIR}/manifest_all.tsv}"

mkdir -p "${META_DIR}"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "Manifest not found: ${MANIFEST_PATH}" >&2
    echo "Run: bash test/cross_trait_scatter/1_generate.sh" >&2
    exit 1
fi

{
    printf 'figure_id\tfigure_kind\tx_id\ty_id\ttable_path\tplot_pdf\tplot_png\tstatus\n'
    while IFS=$'\t' read -r output_id id_x id_y script_path; do
        [ "${output_id}" = "output_id" ] && continue
        output_id="${output_id%$'\r'}"
        id_x="${id_x%$'\r'}"
        id_y="${id_y%$'\r'}"
        status="missing"
        [ -f "${STATUS_DIR}/${output_id}.ok" ] && status="ok"
        [ -f "${STATUS_DIR}/${output_id}.failed" ] && status="failed"
        [ -f "${STATUS_DIR}/${output_id}.running" ] && [ "${status}" = "missing" ] && status="running"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${output_id}" \
            "cross_trait" \
            "${id_x}" \
            "${id_y}" \
            "${TABLES_DIR}/${output_id}.tsv" \
            "${PLOTS_DIR}/${output_id}.pdf" \
            "${PLOTS_DIR}/${output_id}.png" \
            "${status}"
    done < "${MANIFEST_PATH}"
} > "${MANIFEST_OUT}"

echo "Manifest written: ${MANIFEST_OUT}"
