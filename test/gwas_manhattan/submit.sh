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
BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"

STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
LOGS_DIR="${LOGS_DIR:-${RUN_ROOT}/logs}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"

START_ID="${START_ID:-}"
DRY_RUN="${DRY_RUN:-0}"
MAX_PENDING="${MAX_PENDING:-100}"
SQUEUE_USER="${SQUEUE_USER:-${USER}}"

mkdir -p "${STATUS_DIR}" "${LOGS_DIR}" "${FAILURE_DIR}"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "Manifest not found: ${MANIFEST_PATH}" >&2
    echo "Run: bash test/gwas_manhattan/generate.sh" >&2
    exit 1
fi

pending_count() {
    squeue -u "${SQUEUE_USER}" -h -t PD -n "${JOB_NAME}" 2>/dev/null | wc -l
}

SKIP=true
[ -z "${START_ID}" ] && SKIP=false

SUBMITTED=0
SKIPPED=0
FAILED_SUBMIT=0

while IFS=$'\t' read -r gwas_id script_path; do
    [ "${gwas_id}" = "gwas_id" ] && continue

    if [ "${SKIP}" = true ]; then
        if [ "${gwas_id}" = "${START_ID}" ]; then
            SKIP=false
        else
            continue
        fi
    fi

    if [ -f "${STATUS_DIR}/${gwas_id}.ok" ] || [ -f "${STATUS_DIR}/${gwas_id}.running" ]; then
        ((SKIPPED++)) || true
        continue
    fi

    if [ ! -f "${script_path}" ]; then
        printf 'time=%s\nreason=missing_generated_script\nscript_path=%s\n' \
            "$(date -Iseconds)" "${script_path}" > "${STATUS_DIR}/${gwas_id}.failed"
        cp "${STATUS_DIR}/${gwas_id}.failed" "${FAILURE_DIR}/${gwas_id}.failed"
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
        echo "${submit_output}" > "${STATUS_DIR}/${gwas_id}.running"
        rm -f "${FAILURE_DIR}/${gwas_id}.failed"
        echo "[OK] ${gwas_id} -> Job ${submit_output}"
        ((SUBMITTED++)) || true
    else
        echo "[FAIL] ${gwas_id} submit failed"
        printf 'time=%s\nreason=submit_failed\nmessage=%s\n' \
            "$(date -Iseconds)" "${submit_output}" > "${STATUS_DIR}/${gwas_id}.failed"
        cp "${STATUS_DIR}/${gwas_id}.failed" "${FAILURE_DIR}/${gwas_id}.failed"
        ((FAILED_SUBMIT++)) || true
    fi
done < "${MANIFEST_PATH}"

echo ""
echo "Task name: ${TASK_NAME}"
echo "Job name: ${JOB_NAME}"
echo "Submitted: ${SUBMITTED}"
echo "Skipped: ${SKIPPED}"
echo "Submit failed: ${FAILED_SUBMIT}"
echo "Check: bash test/gwas_manhattan/check.sh"
