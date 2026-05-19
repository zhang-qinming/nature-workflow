#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-gwas_manhattan123}"
JOB_NAME="${JOB_NAME:-gman123}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/gwas_manhattan}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"

mapfile -t ALL_IDS < <(tail -n +2 "${FILE_ID_MAP}" | cut -f1)
EXPECTED_TOTAL="${#ALL_IDS[@]}"

ok=0
failed=0
running=0
pdf=0
table_main=0
table_hits=0
failed_files=()

for source_id in "${ALL_IDS[@]}"; do
    if [ -f "${STATUS_DIR}/${source_id}.ok" ]; then
        ((ok++)) || true
    elif [ -f "${STATUS_DIR}/${source_id}.failed" ]; then
        ((failed++)) || true
        failed_files+=("${STATUS_DIR}/${source_id}.failed")
    elif [ -f "${STATUS_DIR}/${source_id}.running" ]; then
        ((running++)) || true
    fi

    [ -f "${OUTPUT_DIR}/plots/${source_id}.pdf" ] && ((pdf++)) || true
    [ -f "${OUTPUT_DIR}/tables/${source_id}.tsv" ] && ((table_main++)) || true
    [ -f "${OUTPUT_DIR}/tables/${source_id}_hits.tsv" ] && ((table_hits++)) || true
done

processed=$((ok + failed + running))
pending=$((EXPECTED_TOTAL - processed))
pd_count=0
r_count=0
if command -v squeue >/dev/null 2>&1; then
    pd_count="$(squeue -u "${USER}" -h -t PD -n "${JOB_NAME}" 2>/dev/null | wc -l)"
    r_count="$(squeue -u "${USER}" -h -t R -n "${JOB_NAME}" 2>/dev/null | wc -l)"
fi

echo "============================================"
echo "GWAS Manhattan 123 status"
echo "============================================"
echo "Task name: ${TASK_NAME}"
echo "Job name:  ${JOB_NAME}"
echo "Expected:  ${EXPECTED_TOTAL}"
echo "OK:        ${ok}"
echo "Failed:    ${failed}"
echo "Running:   ${running}"
echo "Pending:   ${pending}"
echo "PD queue:  ${pd_count}"
echo "R queue:   ${r_count}"
echo ""
echo "PDF:       ${pdf}"
echo "Main TSV:  ${table_main}"
echo "Hits TSV:  ${table_hits}"

if [ "${failed}" -gt 0 ]; then
    echo ""
    echo "Recent failures:"
    ls -1t "${failed_files[@]}" 2>/dev/null | head -10 | while read -r path; do
        echo "--- $(basename "${path}" .failed) ---"
        cat "${path}"
    done
fi
