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
RUN_SCRIPT="${SCRIPT_DIR}/run_one.sh"
STATUS_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/status"
DRY_RUN="${DRY_RUN:-0}"
MAX_CONCURRENT="${MAX_CONCURRENT:-200}"

# 收集要重跑的 ID
FAILED_IDS=()
if [ $# -gt 0 ]; then
    FAILED_IDS=("$@")
else
    for f in "${STATUS_DIR}"/*.failed; do
        [ -f "$f" ] || continue
        FAILED_IDS+=("$(basename "$f" .failed)")
    done
fi

if [ ${#FAILED_IDS[@]} -eq 0 ]; then
    echo "没有需要重跑的任务。"
    exit 0
fi

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
        ((BATCH_COUNT++)) || true
    else
        if job_id=$(sbatch --parsable --export="ALL,GWAS_ID=${gwas_id}" "${RUN_SCRIPT}" 2>/dev/null); then
            echo "[OK] ${gwas_id} -> Job ${job_id}"
            echo "${job_id}" > "${STATUS_DIR}/${gwas_id}.running"
            ((SUBMITTED++)) || true
            ((BATCH_COUNT++)) || true
        else
            echo "[FAIL] ${gwas_id} 提交失败"
        fi
    fi
done

echo ""
echo "已提交 ${SUBMITTED} 个重跑任务"
