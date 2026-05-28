#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TASK_NAME="${TASK_NAME:-cross_trait_heatmap}"
JOB_NAME="${JOB_NAME:-ctheat}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cross_trait_heatmap}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"
RECONCILE_STALE="${RECONCILE_STALE:-1}"

mkdir -p "${STATUS_DIR}" "${FAILURE_DIR}"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "Manifest not found: ${MANIFEST_PATH}" >&2
    echo "Run: bash test/cross_trait_heatmap/1_generate.sh" >&2
    exit 1
fi

ok=0
failed=0
running=0
stale=0
pairs=0
matrix=0
effects=0
pdf=0
png=0
failed_files=()
stale_ids=()
reconciled_stale=0
expected_total=0

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

while IFS=$'\t' read -r output_id lof_ids method script_path; do
    [ "${output_id}" = "output_id" ] && continue
    output_id="${output_id%$'\r'}"
    ((expected_total++)) || true

    if [ -f "${STATUS_DIR}/${output_id}.ok" ]; then
        ((ok++)) || true
    elif [ -f "${STATUS_DIR}/${output_id}.failed" ]; then
        ((failed++)) || true
        failed_files+=("${STATUS_DIR}/${output_id}.failed")
    elif [ -f "${STATUS_DIR}/${output_id}.running" ]; then
        running_job_id="$(tr -d '[:space:]' < "${STATUS_DIR}/${output_id}.running" 2>/dev/null || true)"
        if job_still_active "${running_job_id}"; then
            ((running++)) || true
        else
            stale_ids+=("${output_id}")
            if [ "${RECONCILE_STALE}" = "1" ]; then
                printf 'time=%s\nreason=stale_or_cancelled\njob_id=%s\nmessage=running marker exists but job is no longer visible in squeue\n' \
                    "$(date -Iseconds)" "${running_job_id}" > "${STATUS_DIR}/${output_id}.failed"
                cp "${STATUS_DIR}/${output_id}.failed" "${FAILURE_DIR}/${output_id}.failed"
                rm -f "${STATUS_DIR}/${output_id}.running"
                failed_files+=("${STATUS_DIR}/${output_id}.failed")
                ((failed++)) || true
                ((reconciled_stale++)) || true
            else
                ((stale++)) || true
            fi
        fi
    fi

    [ -f "${OUTPUT_DIR}/tables/${output_id}_pairs.tsv" ] && ((pairs++)) || true
    [ -f "${OUTPUT_DIR}/tables/${output_id}_matrix.tsv" ] && ((matrix++)) || true
    [ -f "${OUTPUT_DIR}/tables/${output_id}_effects.tsv" ] && ((effects++)) || true
    [ -f "${OUTPUT_DIR}/plots/${output_id}.pdf" ] && ((pdf++)) || true
    [ -f "${OUTPUT_DIR}/plots/${output_id}.png" ] && ((png++)) || true
done < "${MANIFEST_PATH}"

processed=$((ok + failed + running + stale))
pending=$((expected_total - processed))
pd_count=0
r_count=0
if command -v squeue >/dev/null 2>&1; then
    pd_count="$(squeue -u "${USER}" -h -t PD -n "${JOB_NAME}" 2>/dev/null | wc -l)"
    r_count="$(squeue -u "${USER}" -h -t R -n "${JOB_NAME}" 2>/dev/null | wc -l)"
fi

echo "============================================"
echo "Cross-trait heatmap status"
echo "============================================"
echo "Task name: ${TASK_NAME}"
echo "Job name:  ${JOB_NAME}"
echo "Expected:  ${expected_total}"
echo "OK:        ${ok}"
echo "Failed:    ${failed}"
echo "Running:   ${running}"
echo "Stale:     ${stale}"
echo "Reconciled:${reconciled_stale}"
echo "Pending:   ${pending}"
echo "PD queue:  ${pd_count}"
echo "R queue:   ${r_count}"
echo ""
echo "Pairs TSV:   ${pairs}"
echo "Matrix TSV:  ${matrix}"
echo "Effects TSV: ${effects}"
echo "PDF:         ${pdf}"
echo "PNG:         ${png}"

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
