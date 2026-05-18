#!/usr/bin/env bash
# ============================================================================
# 检查 GWAS Manhattan 批量任务状态
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_BASE="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs"
STATUS_DIR="${OUTPUT_BASE}/status/gwas_manhattan"
CHUNKS_DIR="${SCRIPT_DIR}/chunks"
LOGS_DIR="${OUTPUT_BASE}/logs/gwas_manhattan"
OUTPUT_DIR="${OUTPUT_BASE}/gwas_manhattan"

echo "============================================"
echo "  GWAS Manhattan 批量任务状态"
echo "============================================"
echo ""

total_chunks=$(ls -1 "${CHUNKS_DIR}"/chunk_*.txt 2>/dev/null | wc -l)

if [ "${total_chunks}" -eq 0 ]; then
    echo "  未找到分块文件，请先运行 generate_and_submit.sh"
    exit 1
fi

ok_count=$(ls -1 "${STATUS_DIR}"/chunk_*.ok 2>/dev/null | wc -l)
failed_count=$(ls -1 "${STATUS_DIR}"/chunk_*.failed 2>/dev/null | wc -l)
pending_count=$((total_chunks - ok_count - failed_count))

echo "  总 chunk 数:  ${total_chunks}"
echo "  已完成:       ${ok_count}"
echo "  失败:         ${failed_count}"
echo "  等待/运行中:  ${pending_count}"
echo ""

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

# 输出目录统计
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
    # 显示最后几行错误日志
    for fail_file in "${STATUS_DIR}"/chunk_*.failed; do
        [ -f "${fail_file}" ] || continue
        chunk_name="$(basename "${fail_file}" .failed)"
        err_log="${LOGS_DIR}/${chunk_name}.err"
        if [ -f "${err_log}" ]; then
            echo "  --- ${chunk_name} 错误日志 (最后 10 行) ---"
            tail -10 "${err_log}" | sed 's/^/    /'
            echo ""
        fi
    done
fi

# 检查队列
if [ "${pending_count}" -gt 0 ] && command -v squeue &>/dev/null; then
    echo "----------------------------------------"
    echo "  队列中的任务:"
    echo "----------------------------------------"
    squeue -n gwas_mht 2>/dev/null | head -20 || true
    echo ""
fi

echo "============================================"
if [ "${failed_count}" -gt 0 ]; then
    echo "  有 ${failed_count} 个 chunk 失败。重跑:"
    echo "    bash test/gwas_manhattan/batch/rerun_failures.sh"
elif [ "${pending_count}" -gt 0 ]; then
    echo "  还有 ${pending_count} 个 chunk 在等待/运行中..."
elif [ "${ok_count}" -eq "${total_chunks}" ]; then
    echo "  全部完成! 构建最终 manifest:"
    echo "    bash test/gwas_manhattan/batch/build_final_manifest.sh"
fi
echo "============================================"
