#!/usr/bin/env bash
# ============================================================================
# 检查 GWAS Manhattan 任务状态
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/gwas_manhattan}"
STATUS_DIR="${STATUS_DIR:-${OUTPUT_DIR}/status}"

mapfile -t ALL_IDS < <(tail -n +2 "${FILE_ID_MAP}" | cut -f1)
EXPECTED_TOTAL="${#ALL_IDS[@]}"

ok=0
failed=0
running=0
failed_files=()

for source_id in "${ALL_IDS[@]}"; do
    if [ -f "${STATUS_DIR}/${source_id}.ok" ]; then
        ((ok++)) || true
    elif [ -f "${STATUS_DIR}/${source_id}.failed" ]; then
        ((failed++)) || true
        failed_files+=("${STATUS_DIR}/${source_id}.failed")
    elif [ -f "${STATUS_DIR}/${source_id}.running" ]; then
        ((running++)) || true
    fi
done

processed=$((ok + failed + running))
pending=$((EXPECTED_TOTAL - processed))

echo "============================================"
echo "  GWAS Manhattan 任务状态"
echo "============================================"
echo "  目标总数:  ${EXPECTED_TOTAL}"
echo "  已完成:  ${ok}"
echo "  失败:    ${failed}"
echo "  运行中:  ${running}"
echo "  未开始:  ${pending}"
echo "  -------------------"
echo "  已处理:  ${processed}"

# 输出文件
pdf=0
tsv=0
hits=0
for source_id in "${ALL_IDS[@]}"; do
    [ -f "${OUTPUT_DIR}/plots/${source_id}.pdf" ] && ((pdf++)) || true
    [ -f "${OUTPUT_DIR}/tables/${source_id}.tsv" ] && ((tsv++)) || true
    [ -f "${OUTPUT_DIR}/tables/${source_id}_hits.tsv" ] && ((hits++)) || true
done
echo ""
echo "  已生成 PDF: ${pdf}"
echo "  已生成主 TSV: ${tsv}"
echo "  已生成 hits TSV: ${hits}"

# 失败详情
if [ "${failed}" -gt 0 ]; then
    echo ""
    echo "----------------------------------------"
    echo "  失败详情 (最后 10 个):"
    echo "----------------------------------------"
    ls -1t "${failed_files[@]}" 2>/dev/null | head -10 | while read -r f; do
        id=$(basename "${f}" .failed)
        info=$(cat "${f}")
        echo "  ${id}  ${info}"
    done
fi

# 集群队列
if command -v squeue &>/dev/null; then
    q=$(squeue -n gman --noheader 2>/dev/null | wc -l)
    echo ""
    echo "  集群队列中: ${q}"
fi

echo ""
echo "============================================"
if [ "${failed}" -gt 0 ]; then
    echo "  重跑失败: bash test/gwas_manhattan/rerun.sh"
elif [ "${EXPECTED_TOTAL}" -gt 0 ] && [ "${ok}" -eq "${EXPECTED_TOTAL}" ] && [ "${failed}" -eq 0 ] && [ "${running}" -eq 0 ]; then
    echo "  全部完成! 构建 manifest:"
    echo "    bash test/gwas_manhattan/manifest.sh"
fi
echo "============================================"
