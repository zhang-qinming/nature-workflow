#!/usr/bin/env bash
# ============================================================================
# 重新提交失败的 GWAS Manhattan chunk
#
# 用法:
#   bash test/batch/rerun_failures.sh
#   bash test/batch/rerun_failures.sh --dry-run
#   bash test/batch/rerun_failures.sh 0001 0005   # 只重跑指定 chunk
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BATCH_DIR="${SCRIPT_DIR}"
CONFIGS_DIR="${BATCH_DIR}/configs"

OUTPUT_BASE="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all"
LOGS_DIR="${OUTPUT_BASE}/logs/gwas_manhattan"
STATUS_DIR="${OUTPUT_BASE}/status/gwas_manhattan"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"

PARTITION="${PARTITION:-cu,fat,batch01,privority}"
CPUS="${CPUS:-2}"
MEM="${MEM:-9G}"
TIME="${TIME:-12:00:00}"

DRY_RUN=0
RERUN_IDS=()

for arg in "${@}"; do
    case "${arg}" in
        --dry-run) DRY_RUN=1 ;;
        --*) echo "未知参数: ${arg}"; exit 1 ;;
        *) RERUN_IDS+=("${arg}") ;;
    esac
done

if [ ${#RERUN_IDS[@]} -gt 0 ]; then
    FAILED_IDS=("${RERUN_IDS[@]}")
    echo "手动指定重跑 ${#FAILED_IDS[@]} 个 chunk: ${FAILED_IDS[*]}"
else
    FAILED_IDS=()
    if [ -d "${STATUS_DIR}" ]; then
        for fail_file in "${STATUS_DIR}"/chunk_*.failed; do
            [ -f "${fail_file}" ] || continue
            chunk_idx="$(basename "${fail_file}" .failed | sed 's/chunk_//')"
            FAILED_IDS+=("${chunk_idx}")
        done
    fi
    if [ ${#FAILED_IDS[@]} -eq 0 ]; then
        echo "没有失败的 chunk。"
        exit 0
    fi
    echo "检测到 ${#FAILED_IDS[@]} 个失败的 chunk: ${FAILED_IDS[*]}"
fi

echo ""

RESUBMIT_DIR="${BATCH_DIR}/resubmit_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${RESUBMIT_DIR}"

JOB_IDS=()

for chunk_idx in "${FAILED_IDS[@]}"; do
    CHUNK_CONFIG="${CONFIGS_DIR}/config.chunk_${chunk_idx}.yaml"
    if [ ! -f "${CHUNK_CONFIG}" ]; then
        echo "[WARN] 配置文件不存在，跳过: ${CHUNK_CONFIG}"
        continue
    fi

    # 清除旧状态
    rm -f "${STATUS_DIR}/chunk_${chunk_idx}.ok" "${STATUS_DIR}/chunk_${chunk_idx}.failed"

    sbatch_script="${RESUBMIT_DIR}/rerun_chunk_${chunk_idx}.sh"
    job_name="gwas_mht_r${chunk_idx}"

    cat > "${sbatch_script}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${job_name}
#SBATCH --error=${LOGS_DIR}/rerun_chunk_${chunk_idx}.err
#SBATCH --output=${LOGS_DIR}/rerun_chunk_${chunk_idx}.out
#SBATCH --mem=${MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

echo "重跑 chunk ${chunk_idx} | \$(date -Iseconds)"
source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"
cd "${PROJECT_ROOT}"

if paper-pipeline run --config "${CHUNK_CONFIG}" figures-gwas-manhattan; then
    echo "\$(date -Iseconds) | SUCCESS (rerun)" > "${STATUS_DIR}/chunk_${chunk_idx}.ok"
    rm -f "${STATUS_DIR}/chunk_${chunk_idx}.failed"
    echo "chunk ${chunk_idx} 重跑成功"
else
    EXIT_CODE=\$?
    echo "\$(date -Iseconds) | EXIT=\${EXIT_CODE} (rerun)" > "${STATUS_DIR}/chunk_${chunk_idx}.failed"
    echo "chunk ${chunk_idx} 再次失败 (exit=\${EXIT_CODE})"
    exit \${EXIT_CODE}
fi
EOF

    chmod +x "${sbatch_script}"

    if [ "${DRY_RUN}" = "1" ]; then
        echo "[DRY RUN] 将提交: sbatch ${sbatch_script}"
    else
        job_id=$(sbatch --parsable "${sbatch_script}")
        JOB_IDS+=("${job_id}")
        echo "[OK] chunk ${chunk_idx} -> Job ${job_id}"
    fi
done

echo ""
if [ "${DRY_RUN}" = "1" ]; then
    echo "[DRY RUN] 未实际提交。去掉 --dry-run 后重新运行。"
elif [ ${#JOB_IDS[@]} -gt 0 ]; then
    echo "已提交 ${#JOB_IDS[@]} 个重跑任务: ${JOB_IDS[*]}"
    echo "监控: squeue -j $(echo "${JOB_IDS[@]}" | tr ' ' ',')"
fi
