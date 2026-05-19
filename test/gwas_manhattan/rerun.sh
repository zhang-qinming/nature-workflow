#!/usr/bin/env bash
# ============================================================================
# 重跑失败的 GWAS ID
#
# 用法:
#   bash test/gwas_manhattan/rerun.sh               # 重跑全部失败
#   bash test/gwas_manhattan/rerun.sh GCST90079736  # 重跑指定 ID
#   DRY_RUN=1 bash test/gwas_manhattan/rerun.sh     # 只列
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/gwas_manhattan}"
RUN_SCRIPT="${SCRIPT_DIR}/run_one.sh"
STATUS_DIR="${STATUS_DIR:-${OUTPUT_DIR}/status}"
LOGS_DIR="${LOGS_DIR:-${OUTPUT_DIR}/logs}"
FAILURE_DIR="${FAILURE_DIR:-${OUTPUT_DIR}/failed}"
DRY_RUN="${DRY_RUN:-0}"
MAX_CONCURRENT="${MAX_CONCURRENT:-200}"

# 收集要重跑的 ID
FAILED_IDS=()
if [ $# -gt 0 ]; then
    FAILED_IDS=("$@")
else
    mapfile -t ALL_IDS < <(tail -n +2 "${FILE_ID_MAP}" | cut -f1)
    for source_id in "${ALL_IDS[@]}"; do
        if [ -f "${STATUS_DIR}/${source_id}.failed" ]; then
            FAILED_IDS+=("${source_id}")
        fi
    done
fi

if [ ${#FAILED_IDS[@]} -eq 0 ]; then
    echo "没有需要重跑的任务。"
    exit 0
fi

mkdir -p "${STATUS_DIR}" "${LOGS_DIR}" "${FAILURE_DIR}"

echo "重跑 ${#FAILED_IDS[@]} 个失败 ID"
echo ""

SUBMITTED=0
BATCH_COUNT=0

for gwas_id in "${FAILED_IDS[@]}"; do
    # 清除旧状态
    rm -f "${STATUS_DIR}/${gwas_id}.ok" "${STATUS_DIR}/${gwas_id}.failed"

    if [ "${BATCH_COUNT}" -ge "${MAX_CONCURRENT}" ]; then
        echo "  ... 已提交 ${SUBMITTED}, 等待队列下降 ..."
        while true; do
            RUNNING=$(squeue -n gman --noheader 2>/dev/null | wc -l)
            [ "${RUNNING}" -lt $((MAX_CONCURRENT / 2)) ] && break
            sleep 10
        done
        BATCH_COUNT=0
    fi

    if [ "${DRY_RUN}" = "1" ]; then
        echo "[DRY RUN] sbatch --export=GWAS_ID=${gwas_id} ${RUN_SCRIPT}"
        ((SUBMITTED++)) || true
        ((BATCH_COUNT++)) || true
    else
        if job_id=$(sbatch --parsable \
            --chdir="${PROJECT_ROOT}" \
            --output="${LOGS_DIR}/gman_%A.out" \
            --error="${LOGS_DIR}/gman_%A.err" \
            --export="ALL,GWAS_ID=${gwas_id},PROJECT_ROOT=${PROJECT_ROOT},OUTPUT_ROOT=${OUTPUT_ROOT},OUTPUT_DIR=${OUTPUT_DIR},STATUS_DIR=${STATUS_DIR}" \
            "${RUN_SCRIPT}" 2>&1); then
            echo "[OK] ${gwas_id} -> Job ${job_id}"
            echo "${job_id}" > "${STATUS_DIR}/${gwas_id}.running"
            rm -f "${FAILURE_DIR}/${gwas_id}.failed"
            ((SUBMITTED++)) || true
            ((BATCH_COUNT++)) || true
        else
            echo "[FAIL] ${gwas_id} 提交失败"
            printf '%s | SUBMIT_FAILED | %s\n' "$(date -Iseconds)" "${job_id}" > "${STATUS_DIR}/${gwas_id}.failed"
            cp "${STATUS_DIR}/${gwas_id}.failed" "${FAILURE_DIR}/${gwas_id}.failed"
        fi
    fi
done

echo ""
echo "已提交 ${SUBMITTED} 个重跑任务"
