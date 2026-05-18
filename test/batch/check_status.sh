#!/usr/bin/env bash
# ============================================================================
# 检查 GWAS Manhattan 批量任务状态
#
# 用法:
#   bash test/batch/check_status.sh
#
# 输出:
#   成功/失败/等待中的 chunk 数量
#   失败 chunk 的详细信息
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATUS_DIR="${SCRIPT_DIR}/status"
CHUNKS_DIR="${SCRIPT_DIR}/chunks"

if [ ! -d "${STATUS_DIR}" ]; then
    echo "状态目录不存在: ${STATUS_DIR}"
    echo "请先运行 generate_and_submit.sh 提交任务"
    exit 1
fi

echo "============================================"
echo "  GWAS Manhattan 批量任务状态"
echo "============================================"
echo ""

# 计数
total_chunks=$(ls -1 "${CHUNKS_DIR}"/chunk_*.txt 2>/dev/null | wc -l)
ok_count=$(ls -1 "${STATUS_DIR}"/chunk_*.ok 2>/dev/null | wc -l)
failed_count=$(ls -1 "${STATUS_DIR}"/chunk_*.failed 2>/dev/null | wc -l)
pending_count=$((total_chunks - ok_count - failed_count))

echo "  总 chunk 数:  ${total_chunks}"
echo "  已完成:       ${ok_count}"
echo "  失败:         ${failed_count}"
echo "  等待/运行中:  ${pending_count}"
echo ""

# 进度条
if [ "${total_chunks}" -gt 0 ]; then
    completed=$((ok_count + failed_count))
    pct=$((completed * 100 / total_chunks))
    bar_len=50
    filled=$((pct * bar_len / 100))
    bar=""
    for ((i = 0; i < filled; i++)); do bar="${bar}#"; done
    for ((i = filled; i < bar_len; i++)); do bar="${bar}-"; done
    echo "  [${bar}] ${completed}/${total_chunks} (${pct}%)"
    echo ""
fi

# 检查输出目录
OUTPUT_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/gwas_manhattan"
if [ -d "${OUTPUT_DIR}/plots" ]; then
    plot_count=$(ls -1 "${OUTPUT_DIR}/plots"/*.pdf 2>/dev/null | wc -l)
    table_count=$(ls -1 "${OUTPUT_DIR}/tables"/*.tsv 2>/dev/null | wc -l)
    echo "  已生成 PDF:   ${plot_count}"
    echo "  已生成 TSV:   ${table_count}"
    echo ""
fi

# 失败详情
if [ "${failed_count}" -gt 0 ]; then
    echo "----------------------------------------"
    echo "  失败 chunk 详情:"
    echo "----------------------------------------"
    for fail_file in "${STATUS_DIR}"/chunk_*.failed; do
        [ -f "${fail_file}" ] || continue
        chunk_name="$(basename "${fail_file}" .failed)"
        content="$(cat "${fail_file}")"
        echo "  ${chunk_name}: ${content}"
    done
    echo ""
fi

# 检查还在排队/运行的 chunk
if [ "${pending_count}" -gt 0 ] && command -v squeue &>/dev/null; then
    echo "----------------------------------------"
    echo "  队列中的任务 (squeue):"
    echo "----------------------------------------"
    squeue -n gwas_mht 2>/dev/null | head -20 || echo "  (无法获取 squeue)"
    echo ""
fi

echo "============================================"
if [ "${failed_count}" -gt 0 ]; then
    echo "  有 ${failed_count} 个 chunk 失败。重新提交:"
    echo "    bash test/batch/rerun_failures.sh"
elif [ "${pending_count}" -gt 0 ]; then
    echo "  还有 ${pending_count} 个 chunk 在等待/运行中..."
elif [ "${ok_count}" -eq "${total_chunks}" ]; then
    echo "  全部 chunk 完成! 构建最终 manifest:"
    echo "    bash test/batch/build_final_manifest.sh"
fi
echo "============================================"
