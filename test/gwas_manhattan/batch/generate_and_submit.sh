#!/usr/bin/env bash
# ============================================================================
# GWAS Manhattan 批量生成 + 提交脚本
#
# 从 file_id_map 提取全部 GWAS ID，按每块 50 个分块，提交 Slurm 数组任务。
#
# 用法:
#   cd /path/to/paper-pipeline
#   bash test/gwas_manhattan/batch/generate_and_submit.sh
#
#   可覆盖的环境变量:
#     PROJECT_ROOT      项目根目录 (默认: 脚本自动探测)
#     CHUNK_SIZE        每块 GWAS ID 数 (默认: 50)
#     CONDA_SH          conda.sh 路径 (默认: ~/miniconda3/etc/profile.d/conda.sh)
#     CONTROL_ENV       paper-pipeline 所在 conda 环境 (默认: paper-pipeline-control)
#     DRY_RUN           1=只生成不提交 (默认: 0)
#
# 输出:
#   chunks/              分块 ID 文件
#   configs/              每块的独立 YAML 配置
#   run_chunk.sbatch      生成的 Slurm 数组脚本
#   logs/                 Slurm 输出日志
#   status/              各任务完成状态
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"
BATCH_DIR="${SCRIPT_DIR}"

CHUNK_SIZE="${CHUNK_SIZE:-50}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
DRY_RUN="${DRY_RUN:-0}"

FILE_ID_MAP="${PROJECT_ROOT}/configs/path.file_id_map.tsv"
BASE_CONFIG="${BATCH_DIR}/config.base.yaml"

# 中间工作文件留在 test/gwas_manhattan/batch 下
CHUNKS_DIR="${BATCH_DIR}/chunks"
CONFIGS_DIR="${BATCH_DIR}/configs"

# 输出统一放在 figure_all/outputs 下
OUTPUT_BASE="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs"
LOGS_DIR="${OUTPUT_BASE}/logs/gwas_manhattan"
STATUS_DIR="${OUTPUT_BASE}/status/gwas_manhattan"

# ---- Slurm 资源配置 ----
PARTITION="${PARTITION:-cu,fat,batch01,privority}"
CPUS="${CPUS:-2}"
MEM="${MEM:-9G}"
TIME="${TIME:-12:00:00}"
JOB_NAME="${JOB_NAME:-gwas_mht}"

echo "============================================"
echo "  GWAS Manhattan 批量生成器"
echo "============================================"
echo "PROJECT_ROOT: ${PROJECT_ROOT}"
echo "BATCH_DIR:    ${BATCH_DIR}"
echo "CHUNK_SIZE:   ${CHUNK_SIZE}"
echo ""

# ---- 清理旧数据 ----
rm -rf "${CHUNKS_DIR}" "${CONFIGS_DIR}" "${STATUS_DIR}"
mkdir -p "${CHUNKS_DIR}" "${CONFIGS_DIR}" "${LOGS_DIR}" "${STATUS_DIR}"

# ---- 1. 提取全部 GWAS ID ----
echo "==> [1/4] 提取 GWAS IDs..."
GWAS_IDS_FILE="${BATCH_DIR}/gwas_ids.txt"
tail -n +2 "${FILE_ID_MAP}" | cut -f1 > "${GWAS_IDS_FILE}"
TOTAL_IDS=$(wc -l < "${GWAS_IDS_FILE}")
echo "  共 ${TOTAL_IDS} 个 GWAS IDs"

# ---- 2. 分块 ----
echo "==> [2/4] 分块 (每块 ${CHUNK_SIZE} 个)..."
split -l "${CHUNK_SIZE}" -d -a 4 "${GWAS_IDS_FILE}" "${CHUNKS_DIR}/chunk_"
# 重命名为 chunk_0000.txt ... chunk_0048.txt
for f in "${CHUNKS_DIR}"/chunk_*; do
    idx="${f##*/chunk_}"
    mv "$f" "${CHUNKS_DIR}/chunk_${idx}.txt"
done
NUM_CHUNKS=$(ls -1 "${CHUNKS_DIR}"/chunk_*.txt 2>/dev/null | wc -l)
echo "  共 ${NUM_CHUNKS} 个分块 (Slurm 数组: 0-$((NUM_CHUNKS - 1)))"

# ---- 3. 生成 per-chunk 配置文件 ----
echo "==> [3/4] 生成 per-chunk 配置文件..."

# 加载 conda 确保 PyYAML 可用
if [ -f "${CONDA_SH}" ]; then
    source "${CONDA_SH}"
    conda activate "${CONTROL_ENV}"
fi

python3 - "${PROJECT_ROOT}" "${BASE_CONFIG}" "${CHUNKS_DIR}" "${CONFIGS_DIR}" <<'PYEOF'
import sys, yaml
from pathlib import Path

project_root = Path(sys.argv[1])
base_config_path = Path(sys.argv[2])
chunks_dir = Path(sys.argv[3])
configs_dir = Path(sys.argv[4])

with open(base_config_path) as f:
    base = yaml.safe_load(f)

# 将 shared_inputs 中的相对路径转为绝对路径，确保从任意位置加载 config 都能找到
shared = base.get("workflows", {}).get("figures", {}).get("shared_inputs", {})
for key in ("file_id_map", "gene_map", "gene_annotation", "geneset_dir"):
    val = shared.get(key)
    if val and isinstance(val, str) and not str(val).startswith("/"):
        shared[key] = str((project_root / val).resolve())

for chunk_file in sorted(chunks_dir.glob("chunk_*.txt")):
    idx = chunk_file.stem.replace("chunk_", "")
    ids = [line.strip() for line in chunk_file.read_text().splitlines() if line.strip()]

    config = dict(base)
    config["project_root"] = str(project_root)
    config["workflows"] = dict(base["workflows"])
    config["workflows"]["figures"] = dict(base["workflows"]["figures"])
    config["workflows"]["figures"]["shared_inputs"] = dict(shared)
    config["workflows"]["figures"]["gwas_manhattan"] = {
        "parameters": {
            "gwas_ids": ids,
            "flank_bp": 50000,
            "label_p_threshold": 1e-30,
            "genomewide_threshold": 5e-8,
            "highlight_genesets": [
                "HALLMARK_HEME_METABOLISM",
                "Hematopoiesisgenes",
                "mitotic_cell_cycle",
                "positive_macromolecule_synthesis",
            ],
        }
    }

    out_path = configs_dir / f"config.chunk_{idx}.yaml"
    with open(out_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

print(f"  生成了 {len(list(configs_dir.glob('config.chunk_*.yaml')))} 个配置文件")
PYEOF

# ---- 4. 生成 Slurm 数组脚本 ----
echo "==> [4/4] 生成 Slurm 数组脚本..."
SBATCH_FILE="${BATCH_DIR}/run_chunk.sbatch"

cat > "${SBATCH_FILE}" <<SBEOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --error=${LOGS_DIR}/chunk_%A_%a.err
#SBATCH --output=${LOGS_DIR}/chunk_%A_%a.out
#SBATCH --mem=${MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --array=0-$((NUM_CHUNKS - 1))
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

TASK_ID=\${SLURM_ARRAY_TASK_ID}
CHUNK_CONFIG="${CONFIGS_DIR}/config.chunk_\${TASK_ID}.yaml"
STATUS_OK="${STATUS_DIR}/chunk_\${TASK_ID}.ok"
STATUS_FAIL="${STATUS_DIR}/chunk_\${TASK_ID}.failed"

echo "============================================"
echo "  Array Task \${TASK_ID} 开始"
echo "  Config: \${CHUNK_CONFIG}"
echo "  Time:   \$(date -Iseconds)"
echo "============================================"

if [ ! -f "\${CHUNK_CONFIG}" ]; then
    echo "[FAIL] 配置文件不存在: \${CHUNK_CONFIG}"
    echo "\$(date -Iseconds) | MISSING CONFIG: \${CHUNK_CONFIG}" > "\${STATUS_FAIL}"
    exit 1
fi

# 加载 conda
source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"

# 切到项目根目录
cd "${PROJECT_ROOT}"

# 执行
if paper-pipeline run --config "\${CHUNK_CONFIG}" figures-gwas-manhattan; then
    echo ""
    echo "============================================"
    echo "  Array Task \${TASK_ID} 完成"
    echo "  Time: \$(date -Iseconds)"
    echo "============================================"
    echo "\$(date -Iseconds) | SUCCESS" > "\${STATUS_OK}"
    rm -f "\${STATUS_FAIL}"
else
    EXIT_CODE=\$?
    echo ""
    echo "============================================"
    echo "  Array Task \${TASK_ID} 失败 (exit=\${EXIT_CODE})"
    echo "  Time: \$(date -Iseconds)"
    echo "============================================"
    echo "\$(date -Iseconds) | EXIT=\${EXIT_CODE}" > "\${STATUS_FAIL}"
    exit \${EXIT_CODE}
fi
SBEOF

chmod +x "${SBATCH_FILE}"

echo ""
echo "============================================"
echo "  生成完毕"
echo "============================================"
echo ""
echo "  总 GWAS IDs:  ${TOTAL_IDS}"
echo "  分块数:       ${NUM_CHUNKS} (每块 ${CHUNK_SIZE} 个)"
echo "  Slurm 脚本:   ${SBATCH_FILE}"
echo "  输出根目录:   ${OUTPUT_BASE}/"
echo "    plots/tables: ${OUTPUT_BASE}/gwas_manhattan/"
echo "    logs:         ${LOGS_DIR}/"
echo "    status:       ${STATUS_DIR}/"
echo ""

if [ "${DRY_RUN}" = "1" ]; then
    echo "[DRY RUN] 跳过提交。手动提交:"
    echo "  sbatch ${SBATCH_FILE}"
else
    echo "提交 Slurm 数组任务..."
    sbatch "${SBATCH_FILE}"
    echo ""
    echo "查看状态:"
    echo "  squeue -n ${JOB_NAME}"
    echo "  bash test/gwas_manhattan/batch/check_status.sh"
    echo ""
    echo "重新提交失败的 chunk:"
    echo "  bash test/gwas_manhattan/batch/rerun_failures.sh"
fi
