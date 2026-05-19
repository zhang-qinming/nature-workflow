#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-gwas_manhattan}"
JOB_NAME="${JOB_NAME:-gman}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/gwas_manhattan}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"
RECONCILE_STALE="${RECONCILE_STALE:-1}"

mkdir -p "${STATUS_DIR}" "${FAILURE_DIR}"

mapfile -t ALL_IDS < <(awk -F '\t' 'NR > 1 { sub(/\r$/, "", $1); print $1 }' "${FILE_ID_MAP}")
EXPECTED_TOTAL="${#ALL_IDS[@]}"

ok=0
failed=0
running=0
stale=0
pdf=0
table_main=0
table_hits=0
failed_files=()
stale_ids=()
reconciled_stale=0

job_still_active() {
    local job_id="$1"
    if ! command -v squeue >/dev/null 2>&1; then
        return 0
    fi
    if [ -z "${job_id}" ]; then
        return 1
    fi
    squeue -j "${job_id}" -h 2>/dev/null | grep -q .
}

for source_id in "${ALL_IDS[@]}"; do
    if [ -f "${STATUS_DIR}/${source_id}.ok" ]; then
        ((ok++)) || true
    elif [ -f "${STATUS_DIR}/${source_id}.failed" ]; then
        ((failed++)) || true
        failed_files+=("${STATUS_DIR}/${source_id}.failed")
    elif [ -f "${STATUS_DIR}/${source_id}.running" ]; then
        running_job_id="$(tr -d '[:space:]' < "${STATUS_DIR}/${source_id}.running" 2>/dev/null || true)"
        if job_still_active "${running_job_id}"; then
            ((running++)) || true
        else
            stale_ids+=("${source_id}")
            if [ "${RECONCILE_STALE}" = "1" ]; then
                printf 'time=%s\nreason=stale_or_cancelled\njob_id=%s\nmessage=running marker exists but job is no longer visible in squeue\n' \
                    "$(date -Iseconds)" "${running_job_id}" > "${STATUS_DIR}/${source_id}.failed"
                cp "${STATUS_DIR}/${source_id}.failed" "${FAILURE_DIR}/${source_id}.failed"
                rm -f "${STATUS_DIR}/${source_id}.running"
                failed_files+=("${STATUS_DIR}/${source_id}.failed")
                ((failed++)) || true
                ((reconciled_stale++)) || true
            else
                ((stale++)) || true
            fi
        fi
    fi

    [ -f "${OUTPUT_DIR}/plots/${source_id}.pdf" ] && ((pdf++)) || true
    [ -f "${OUTPUT_DIR}/tables/${source_id}_variants.tsv" ] && ((table_main++)) || true
    [ -f "${OUTPUT_DIR}/tables/${source_id}_hits.tsv" ] && ((table_hits++)) || true
done

processed=$((ok + failed + running + stale))
pending=$((EXPECTED_TOTAL - processed))
pd_count=0
r_count=0
if command -v squeue >/dev/null 2>&1; then
    pd_count="$(squeue -u "${USER}" -h -t PD -n "${JOB_NAME}" 2>/dev/null | wc -l)"
    r_count="$(squeue -u "${USER}" -h -t R -n "${JOB_NAME}" 2>/dev/null | wc -l)"
fi

echo "============================================"
echo "GWAS Manhattan status"
echo "============================================"
echo "Task name: ${TASK_NAME}"
echo "Job name:  ${JOB_NAME}"
echo "Expected:  ${EXPECTED_TOTAL}"
echo "OK:        ${ok}"
echo "Failed:    ${failed}"
echo "Running:   ${running}"
echo "Stale:     ${stale}"
echo "Reconciled:${reconciled_stale}"
echo "Pending:   ${pending}"
echo "PD queue:  ${pd_count}"
echo "R queue:   ${r_count}"
echo ""
echo "PDF:       ${pdf}"
echo "Variants TSV: ${table_main}"
echo "Hits TSV:     ${table_hits}"

if [ "${failed}" -gt 0 ]; then
    echo ""
    echo "Recent failures:"
    ls -1t "${failed_files[@]}" 2>/dev/null | head -10 | while read -r path; do
        echo "--- $(basename "${path}" .failed) ---"
        cat "${path}"
    done
fi

if [ "${stale}" -gt 0 ]; then
    echo ""
    echo "Stale IDs:"
    printf '%s\n' "${stale_ids[@]}" | head -20
fi
