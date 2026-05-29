#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TASK_NAME="${TASK_NAME:-cross_trait_heatmap}"
JOB_NAME="${JOB_NAME:-ctheat}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"

STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
LOGS_DIR="${LOGS_DIR:-${RUN_ROOT}/logs}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"

DRY_RUN="${DRY_RUN:-0}"
MAX_PENDING="${MAX_PENDING:-100}"
SQUEUE_USER="${SQUEUE_USER:-${USER}}"

mkdir -p "${STATUS_DIR}" "${LOGS_DIR}" "${FAILURE_DIR}"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "Manifest not found: ${MANIFEST_PATH}" >&2
    echo "Run: bash test/cross_trait_heatmap/1_generate.sh" >&2
    exit 1
fi

pending_count() {
    squeue -u "${SQUEUE_USER}" -h -t PD -n "${JOB_NAME}" 2>/dev/null | wc -l
}

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

START_ID="${START_ID:-}"
SKIP=true
[ -z "${START_ID}" ] && SKIP=false

SUBMITTED=0
SKIPPED=0
FAILED_SUBMIT=0

while IFS=$'\t' read -r lof_id script_path; do
    [ "${lof_id}" = "lof_id" ] && continue
    lof_id="${lof_id%$'\r'}"
    script_path="${script_path%$'\r'}"

    if [ "${SKIP}" = true ]; then
        if [ "${lof_id}" = "${START_ID}" ]; then
            SKIP=false
        else
            continue
        fi
    fi

    if [ -f "${STATUS_DIR}/${lof_id}.ok" ]; then
        ((SKIPPED++)) || true
        continue
    fi

    if [ -f "${STATUS_DIR}/${lof_id}.running" ]; then
        running_job_id="$(tr -d '[:space:]' < "${STATUS_DIR}/${lof_id}.running" 2>/dev/null || true)"
        if job_still_active "${running_job_id}"; then
            ((SKIPPED++)) || true
            continue
        fi
        rm -f "${STATUS_DIR}/${lof_id}.running"
    fi

    if [ ! -f "${script_path}" ]; then
        printf 'time=%s\nreason=missing_generated_script\nscript_path=%s\n' \
            "$(date -Iseconds)" "${script_path}" > "${STATUS_DIR}/${lof_id}.failed"
        cp "${STATUS_DIR}/${lof_id}.failed" "${FAILURE_DIR}/${lof_id}.failed"
        ((FAILED_SUBMIT++)) || true
        continue
    fi

    while [ "$(pending_count)" -ge "${MAX_PENDING}" ]; do
        sleep 10
    done

    if [ "${DRY_RUN}" = "1" ]; then
        echo "[DRY RUN] sbatch ${script_path}"
        ((SUBMITTED++)) || true
        continue
    fi

    submit_output=""
    if submit_output="$(sbatch --parsable "${script_path}" 2>&1)"; then
        echo "${submit_output}" > "${STATUS_DIR}/${lof_id}.running"
        rm -f "${FAILURE_DIR}/${lof_id}.failed"
        echo "[OK] ${lof_id} -> Job ${submit_output}"
        ((SUBMITTED++)) || true
    else
        echo "[FAIL] ${lof_id} submit failed"
        printf 'time=%s\nreason=submit_failed\nmessage=%s\n' \
            "$(date -Iseconds)" "${submit_output}" > "${STATUS_DIR}/${lof_id}.failed"
        cp "${STATUS_DIR}/${lof_id}.failed" "${FAILURE_DIR}/${lof_id}.failed"
        ((FAILED_SUBMIT++)) || true
    fi
done < "${MANIFEST_PATH}"

echo ""
echo "Task name: ${TASK_NAME}"
echo "Job name: ${JOB_NAME}"
echo "Submitted: ${SUBMITTED}"
echo "Skipped: ${SKIPPED}"
echo "Submit failed: ${FAILED_SUBMIT}"
echo "Check: bash test/cross_trait_heatmap/3_check.sh"
