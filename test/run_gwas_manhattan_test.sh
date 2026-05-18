#!/usr/bin/env bash
# ============================================================================
# GWAS Manhattan 图测试脚本
# 在 Linux 服务器上运行，测试 5 个 GWAS trait 的曼哈顿图生成。
#
# 用法:
#   PROJECT_ROOT=/path/to/paper-pipeline bash test/run_gwas_manhattan_test.sh
#
#   或直接 cd 到项目根目录:
#   cd /path/to/paper-pipeline && bash test/run_gwas_manhattan_test.sh
#
# 输入:
#   分 trait 跑的文件: 通过 file_id_map (configs/path.file_id_map.tsv) 引用
#     实际 GWAS 文件路径: /gpfs/chencao/qinminzhang/workflow/ldsc/GWAS/all/GCST*.txt.gz
#   只用跑一次的参考文件:
#     /gpfs/chencao/qinminzhang/Nature/mine/code/outputs/
#       ├── genes.protein_coding.v39.gtf
#       ├── gencode_v41_gname_gid_ALL_sorted_onlyID
#       └── geneset/
#
# 输出:
#   /gpfs/chencao/qinminzhang/workflow/catalog_lof/fgtest/outputs/gwas_manhattan/
#       ├── tables/
#       ├── plots/   (*.pdf, *.png)
#       └── meta/    (manifest.tsv)
# ============================================================================

set -euo pipefail

# ---- 项目根目录 ----
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
echo "==> PROJECT_ROOT: ${PROJECT_ROOT}"

CONFIG_FILE="${PROJECT_ROOT}/test/gwas_manhattan_test.yaml"
FILE_ID_MAP="${PROJECT_ROOT}/configs/path.file_id_map.tsv"

# ---- 共享参考文件路径 ----
REF_BASE="/gpfs/chencao/qinminzhang/Nature/mine/code/outputs"
GENE_ANNOTATION="${REF_BASE}/genes.protein_coding.v39.gtf"
GENE_MAP="${REF_BASE}/gencode_v41_gname_gid_ALL_sorted_onlyID"
GENESET_DIR="${REF_BASE}/geneset"

# ---- 需要检查的 GWAS IDs ----
GWAS_IDS=(
    GCST90079736
    GCST90079737
    GCST90080060
    GCST90079734
    GCST90081171
)

# ---- 输出目录 ----
OUTPUT_BASE="/gpfs/chencao/qinminzhang/workflow/catalog_lof/fgtest/outputs"

# ============================================================================
echo ""
echo "============ 环境检查 ============"

# 检查 conda
if ! command -v conda &>/dev/null; then
    echo "[FAIL] conda 未找到，请先加载 conda 环境"
    exit 1
fi
echo "[ OK ] conda: $(conda --version 2>/dev/null || echo 'installed')"

# 检查 paper-pipeline
if ! conda run -n paper-pipeline-control paper-pipeline --help &>/dev/null; then
    echo "[WARN] paper-pipeline 未安装，正在以 editable 模式安装..."
    conda run -n paper-pipeline-control pip install -e "${PROJECT_ROOT}" --no-deps
    echo "[ OK ] paper-pipeline 安装完成"
else
    echo "[ OK ] paper-pipeline 已安装"
fi

# 检查配置文件
echo ""
echo "============ 检查输入文件 ============"
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "[FAIL] 配置文件不存在: ${CONFIG_FILE}"
    exit 1
fi
echo "[ OK ] 配置文件: ${CONFIG_FILE}"

if [ ! -f "${FILE_ID_MAP}" ]; then
    echo "[FAIL] file_id_map 不存在: ${FILE_ID_MAP}"
    exit 1
fi
echo "[ OK ] file_id_map: ${FILE_ID_MAP} (共 $(wc -l < "${FILE_ID_MAP}") 行)"

# 检查参考文件
echo ""
echo "--- 参考文件 (只用跑一次) ---"
check_file() {
    local path="$1"
    local label="$2"
    if [ -f "${path}" ]; then
        echo "[ OK ] ${label}: ${path}"
    elif [ -d "${path}" ]; then
        local count=$(ls -1 "${path}" | wc -l)
        echo "[ OK ] ${label}: ${path} (${count} 个文件)"
    else
        echo "[FAIL] ${label}: ${path} 不存在"
        return 1
    fi
}

REF_OK=true
check_file "${GENE_ANNOTATION}" "gene_annotation" || REF_OK=false
check_file "${GENE_MAP}"        "gene_map"        || REF_OK=false
check_file "${GENESET_DIR}"     "geneset_dir"     || REF_OK=false

if [ "${REF_OK}" = false ]; then
    echo ""
    echo "[WARN] 部分参考文件缺失，请检查 ${REF_BASE} 目录内容:"
    ls -la "${REF_BASE}/" 2>/dev/null || echo "  (无法列出目录)"
    echo ""
    echo "  期望的参考文件:"
    echo "    - genes.protein_coding.v39.gtf     (基因注释 GTF)"
    echo "    - gencode_v41_gname_gid_ALL_sorted_onlyID  (基因 ID 映射)"
    echo "    - geneset/                          (基因集目录，含 HALLMARK_*.txt 等)"
    echo ""
    echo "  如果文件在子目录中 (如 GWAS/genes.protein_coding.v39.gtf)，"
    echo "  请修改本脚本或配置文件中的路径。"
    exit 1
fi

# 检查 GWAS 文件 (通过 file_id_map 查找)
echo ""
echo "--- GWAS 文件 (分 trait 跑) ---"
GWAS_OK=true
for gwas_id in "${GWAS_IDS[@]}"; do
    gwas_path=$(awk -F'\t' -v id="${gwas_id}" '$1==id {print $3; exit}' "${FILE_ID_MAP}")
    if [ -z "${gwas_path}" ]; then
        echo "[FAIL] ${gwas_id}: 在 file_id_map 中未找到"
        GWAS_OK=false
    elif [ -f "${gwas_path}" ]; then
        size=$(du -h "${gwas_path}" | cut -f1)
        echo "[ OK ] ${gwas_id}: ${gwas_path} (${size})"
    else
        echo "[FAIL] ${gwas_id}: ${gwas_path} 文件不存在"
        GWAS_OK=false
    fi
done

if [ "${GWAS_OK}" = false ]; then
    echo ""
    echo "[WARN] 部分 GWAS 文件缺失。可用的 GWAS IDs (前 20 个):"
    head -21 "${FILE_ID_MAP}" | awk -F'\t' '{print "  " $1 "  ->  " $3}'
    exit 1
fi

# 检查输出目录
echo ""
echo "--- 输出目录 ---"
mkdir -p "${OUTPUT_BASE}"
echo "[ OK ] 输出根目录: ${OUTPUT_BASE}"

# ============================================================================
echo ""
echo "============ 开始运行 ============"
echo ""

cd "${PROJECT_ROOT}"

# Step 1: plan (预览)
echo "--- Step 1: plan (预览任务) ---"
conda run -n paper-pipeline-control paper-pipeline plan \
    --config "${CONFIG_FILE}" \
    figures-gwas-manhattan

echo ""
echo "--- Step 2: run (执行) ---"
conda run -n paper-pipeline-control paper-pipeline run \
    --config "${CONFIG_FILE}" \
    figures-gwas-manhattan

# ============================================================================
echo ""
echo "============ 结果检查 ============"
OUTPUT_DIR="${OUTPUT_BASE}/gwas_manhattan"
TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"

for dir in "${TABLES_DIR}" "${PLOTS_DIR}" "${META_DIR}"; do
    if [ -d "${dir}" ]; then
        count=$(find "${dir}" -type f | wc -l)
        echo "[ OK ] ${dir} (${count} 个文件)"
    else
        echo "[WARN] ${dir} 未生成"
    fi
done

echo ""
echo "--- 生成的文件 ---"
if [ -d "${OUTPUT_DIR}" ]; then
    find "${OUTPUT_DIR}" -type f | sort | while read -r f; do
        echo "  ${f}"
    done
else
    echo "  (输出目录未生成，请检查错误日志)"
    exit 1
fi

echo ""
echo "============ 测试完成 ============"
echo "PDF 曼哈顿图: ${PLOTS_DIR}/*.pdf"
echo "PNG 曼哈顿图: ${PLOTS_DIR}/*.png"
echo "TSV 数据表:   ${TABLES_DIR}/*.tsv"
echo "清单文件:     ${META_DIR}/manifest.tsv"
