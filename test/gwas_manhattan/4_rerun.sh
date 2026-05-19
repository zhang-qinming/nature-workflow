#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"

TASK_NAME="${TASK_NAME:-gwas_manhattan}"
JOB_NAME="${JOB_NAME:-gman}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"

STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"

DRY_RUN="${DRY_RUN:-0}"
MAX_PENDING="${MAX_PENDING:-100}"
RERUN_MEM="${RERUN_MEM:-32G}"
RERUN_CPUS="${RERUN_CPUS:-2}"

if [ ! -f "${MANIFEST_PATH}" ]; then
    echo "Manifest not found: ${MANIFEST_PATH}" >&2
    echo "Run: bash test/gwas_manhattan/1_generate.sh" >&2
    exit 1
fi

declare -A SCRIPT_BY_ID
while IFS=$'\t' read -r gwas_id script_path; do
    if [ "${gwas_id}" = "gwas_id" ]; then
        continue
    fi
    gwas_id="${gwas_id%$'\r'}"
    script_path="${script_path%$'\r'}"
    SCRIPT_BY_ID["${gwas_id}"]="${script_path}"
done < "${MANIFEST_PATH}"

FAILED_IDS=()
if [ "$#" -gt 0 ]; then
    FAILED_IDS=("$@")
else
    for path in "${STATUS_DIR}"/*.failed; do
        [ -f "${path}" ] || continue
        FAILED_IDS+=("$(basename "${path}" .failed)")
    done
fi

if [ "${#FAILED_IDS[@]}" -eq 0 ]; then
    echo "No failed IDs to rerun."
    exit 0
fi

pending_count() {
    squeue -u "${USER}" -h -t PD -n "${JOB_NAME}" 2>/dev/null | wc -l
}

SUBMITTED=0

for gwas_id in "${FAILED_IDS[@]}"; do
    script_path="${SCRIPT_BY_ID[${gwas_id}]:-}"
    if [ -z "${script_path}" ] || [ ! -f "${script_path}" ]; then
        echo "[SKIP] missing generated script for ${gwas_id}"
        continue
    fi

    rm -f "${STATUS_DIR}/${gwas_id}.ok" "${STATUS_DIR}/${gwas_id}.failed" "${STATUS_DIR}/${gwas_id}.running" "${FAILURE_DIR}/${gwas_id}.failed"

    while [ "$(pending_count)" -ge "${MAX_PENDING}" ]; do
        sleep 10
    done

    if [ "${DRY_RUN}" = "1" ]; then
        echo "[DRY RUN] sbatch --mem=${RERUN_MEM} --cpus-per-task=${RERUN_CPUS} ${script_path}"
        ((SUBMITTED++)) || true
        continue
    fi

    submit_output=""
    if submit_output="$(sbatch --parsable --mem="${RERUN_MEM}" --cpus-per-task="${RERUN_CPUS}" "${script_path}" 2>&1)"; then
        echo "${submit_output}" > "${STATUS_DIR}/${gwas_id}.running"
        echo "[OK] ${gwas_id} -> Job ${submit_output}"
        ((SUBMITTED++)) || true
    else
        printf 'time=%s\nreason=submit_failed\nmessage=%s\n' \
            "$(date -Iseconds)" "${submit_output}" > "${STATUS_DIR}/${gwas_id}.failed"
        cp "${STATUS_DIR}/${gwas_id}.failed" "${FAILURE_DIR}/${gwas_id}.failed"
        echo "[FAIL] ${gwas_id} submit failed"
    fi
done

echo "Rerun submitted: ${SUBMITTED}"
