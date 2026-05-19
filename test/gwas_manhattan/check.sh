#!/usr/bin/env bash
# ============================================================================
# 检查 GWAS Manhattan 任务状态
# ============================================================================
set -euo pipefail

STATUS_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/status"
OUTPUT_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/gwas_manhattan"

ok=$(ls -1 "${STATUS_DIR}"/*.ok      2>/dev/null | wc -l)
failed=$(ls -1 "${STATUS_DIR}"/*.failed 2>/dev/null | wc -l)
running=$(ls -1 "${STATUS_DIR}"/*.running 2>/dev/null | wc -l)
total=$((ok + failed + running))

echo "============================================"
echo "  GWAS Manhattan 任务状态"
echo "============================================"
echo "  已完成:  ${ok}"
echo "  失败:    ${failed}"
echo "  运行中:  ${running}"
echo "  -------------------"
echo "  已处理:  $((ok + failed + running))"

# 输出文件
if [ -d "${OUTPUT_DIR}/plots" ]; then
    pdf=$(ls -1 "${OUTPUT_DIR}/plots"/*.pdf 2>/dev/null | wc -l)
    tsv=$(ls -1 "${OUTPUT_DIR}/tables"/*.tsv 2>/dev/null | wc -l)
    echo ""
    echo "  已生成 PDF: ${pdf}"
    echo "  已生成 TSV: ${tsv}"
fi

# 失败详情
if [ "${failed}" -gt 0 ]; then
    echo ""
    echo "----------------------------------------"
    echo "  失败详情 (最后 10 个):"
    echo "----------------------------------------"
    ls -1t "${STATUS_DIR}"/*.failed 2>/dev/null | head -10 | while read f; do
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
elif [ "${ok}" -gt 0 ] && [ "${failed}" -eq 0 ] && [ "${running}" -eq 0 ]; then
    echo "  全部完成! 构建 manifest:"
    echo "    bash test/gwas_manhattan/manifest.sh"
fi
echo "============================================"
