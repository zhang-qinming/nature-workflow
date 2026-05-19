#!/usr/bin/env bash
# ============================================================================
# 提交全部 GWAS Manhattan 任务 — 1 ID = 1 Slurm Job
#
# 用法:
#   bash test/gwas_manhattan/submit.sh
#   DRY_RUN=1 bash test/gwas_manhattan/submit.sh  # 只列不提交
#
# 可调参数:
#   MAX_CONCURRENT  同时运行的最大 Job 数 (默认 200，超过后等待)
#   START_ID        从指定 GWAS ID 开始 (默认: 从头开始)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${PROJECT_ROOT}/configs/path.file_id_map.tsv"
RUN_SCRIPT="${SCRIPT_DIR}/run_one.sh"
STATUS_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/status"
LOGS_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/logs"

MAX_CONCURRENT="${MAX_CONCURRENT:-200}"
DRY_RUN="${DRY_RUN:-0}"
START_ID="${START_ID:-}"

mkdir -p "${STATUS_DIR}" "${LOGS_DIR}"

# ---- 读取全部 GWAS ID ----
echo "==> 读取 GWAS IDs..."
mapfile -t ALL_IDS < <(tail -n +2 "${FILE_ID_MAP}" | cut -f1)
echo "  共 ${#ALL_IDS[@]} 个"

# 如果指定了起始 ID，跳过之前的
SKIP=true
if [ -z "${START_ID}" ]; then
    SKIP=false
fi

SUBMITTED=0
SKIPPED=0
FAILED_SUBMIT=0

echo ""
echo "==> 开始提交 (最多 ${MAX_CONCURRENT} 并发)..."
echo ""

BATCH_COUNT=0

for gwas_id in "${ALL_IDS[@]}"; do
    # 起始 ID 过滤
    if [ "${SKIP}" = true ]; then
        if [ "${gwas_id}" = "${START_ID}" ]; then
            SKIP=false
        else
            ((SKIPPED++)) || true
            continue
        fi
    fi

    # 跳过已成功 / 运行中
    if [ -f "${STATUS_DIR}/${gwas_id}.ok" ] || [ -f "${STATUS_DIR}/${gwas_id}.running" ]; then
        continue
    fi

    # 每攒够 MAX_CONCURRENT 个，等队列水位降到一半再继续
    if [ "${BATCH_COUNT}" -ge "${MAX_CONCURRENT}" ]; then
        echo "  ... 已提交 ${SUBMITTED}, 等待队列降到 $((MAX_CONCURRENT / 2)) ..."
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
        if job_id=$(sbatch --parsable --export="ALL,GWAS_ID=${gwas_id}" "${RUN_SCRIPT}" 2>/dev/null); then
            echo "[${SUBMITTED}] ${gwas_id} -> Job ${job_id}"
            echo "${job_id}" > "${STATUS_DIR}/${gwas_id}.running"
            ((SUBMITTED++)) || true
            ((BATCH_COUNT++)) || true
        else
            echo "[FAIL] ${gwas_id} 提交失败"
            ((FAILED_SUBMIT++)) || true
        fi
    fi
done

echo ""
echo "============================================"
echo "  提交完毕"
echo "============================================"
echo "  提交成功:    ${SUBMITTED}"
echo "  提交失败:    ${FAILED_SUBMIT}"
echo "  跳过:        ${SKIPPED}"
echo ""
echo "  查看状态:    bash test/gwas_manhattan/check.sh"
echo "  重跑失败:    bash test/gwas_manhattan/rerun.sh"
echo "  队列中:      squeue -n gman"
