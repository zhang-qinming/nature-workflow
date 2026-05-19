#!/usr/bin/env bash
# ============================================================================
# Slurm 单任务模板 — 处理 1 个 GWAS ID 的曼哈顿图
# 由 submit.sh 调用: sbatch --export=GWAS_ID=xxx run_one.sh
# ============================================================================
#SBATCH --job-name=gman
#SBATCH --output=/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/logs/gman_%A.out
#SBATCH --error=/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/logs/gman_%A.err
#SBATCH --mem=9G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=cu,fat,batch01,privority
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/chencao/qinminzhang/Nature/mine/code

set -euo pipefail

GWAS_ID="${GWAS_ID:?GWAS_ID not set}"
PROJECT_ROOT="${PROJECT_ROOT:-/gpfs/chencao/qinminzhang/Nature/mine/code}"
BASE_CONFIG="${PROJECT_ROOT}/test/gwas_manhattan/config.base.yaml"
TEMP_CONFIG="/tmp/config_gman_${GWAS_ID}.yaml"
STATUS_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/status"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"

echo "============================================"
echo "  GWAS ID: ${GWAS_ID}"
echo "  Job ID:  ${SLURM_JOB_ID:-local}"
echo "  Host:    $(hostname)"
echo "  Time:    $(date -Iseconds)"
echo "============================================"

# ---- 加载 conda ----
if [ -f "${CONDA_SH}" ]; then
    source "${CONDA_SH}"
    conda activate "${CONTROL_ENV}"
fi

# ---- 生成单 ID 配置 ----
python3 - "${GWAS_ID}" "${PROJECT_ROOT}" "${BASE_CONFIG}" "${TEMP_CONFIG}" <<'PYEOF'
import sys, yaml
from pathlib import Path

gwas_id    = sys.argv[1]
proj_root  = Path(sys.argv[2])
base_path  = Path(sys.argv[3])
temp_path  = Path(sys.argv[4])

with open(base_path) as f:
    config = yaml.safe_load(f)

config["project_root"] = str(proj_root)

# 将共享输入中的相对路径转绝对
shared = config.get("workflows", {}).get("figures", {}).get("shared_inputs", {})
for key in ("file_id_map", "gene_map", "gene_annotation", "geneset_dir"):
    val = shared.get(key)
    if val and isinstance(val, str) and not str(val).startswith("/"):
        shared[key] = str((proj_root / val).resolve())

config["workflows"]["figures"]["shared_inputs"] = shared
config["workflows"]["figures"]["gwas_manhattan"] = {
    "parameters": {
        "gwas_ids": [gwas_id],
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

temp_path.parent.mkdir(parents=True, exist_ok=True)
with open(temp_path, "w") as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
print(f"  Config written: {temp_path}")
PYEOF

# ---- 执行 ----
echo "  Running paper-pipeline..."
mkdir -p "${STATUS_DIR}"

if paper-pipeline run --config "${TEMP_CONFIG}" figures-gwas-manhattan; then
    echo "$(date -Iseconds) | SUCCESS" > "${STATUS_DIR}/${GWAS_ID}.ok"
    rm -f "${STATUS_DIR}/${GWAS_ID}.failed" "${STATUS_DIR}/${GWAS_ID}.running"
    echo "  => SUCCESS"
else
    EXIT_CODE=$?
    echo "$(date -Iseconds) | EXIT=${EXIT_CODE}" > "${STATUS_DIR}/${GWAS_ID}.failed"
    rm -f "${STATUS_DIR}/${GWAS_ID}.running"
    echo "  => FAILED (exit=${EXIT_CODE})"
    exit ${EXIT_CODE}
fi

# ---- 清理 ----
rm -f "${TEMP_CONFIG}"
