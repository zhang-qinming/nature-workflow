#!/usr/bin/env bash
# ============================================================================
# 全部 chunk 完成后构建最终 manifest + 完整性检查
# ============================================================================

set -euo pipefail

OUTPUT_BASE="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs"
OUTPUT_DIR="${OUTPUT_BASE}/gwas_manhattan"
TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"
STATUS_DIR="${OUTPUT_BASE}/status/gwas_manhattan"

echo "============================================"
echo "  构建最终合并 Manifest"
echo "============================================"

# 检查状态
if [ -d "${STATUS_DIR}" ]; then
    total=$(ls -1 "${STATUS_DIR}"/chunk_*.ok "${STATUS_DIR}"/chunk_*.failed 2>/dev/null | wc -l)
    ok=$(ls -1 "${STATUS_DIR}"/chunk_*.ok 2>/dev/null | wc -l)
    failed=$(ls -1 "${STATUS_DIR}"/chunk_*.failed 2>/dev/null | wc -l)
    echo "  状态: ${ok} 成功 / ${failed} 失败 / ${total} 总"
fi

if [ ! -d "${TABLES_DIR}" ]; then
    echo "[FAIL] tables 目录不存在: ${TABLES_DIR}"
    exit 1
fi

TABLE_COUNT=$(ls -1 "${TABLES_DIR}"/*.tsv 2>/dev/null | wc -l)
PLOT_PDF_COUNT=$(ls -1 "${PLOTS_DIR}"/*.pdf 2>/dev/null | wc -l)
PLOT_PNG_COUNT=$(ls -1 "${PLOTS_DIR}"/*.png 2>/dev/null | wc -l)

echo "  Tables (TSV): ${TABLE_COUNT}"
echo "  Plots  (PDF): ${PLOT_PDF_COUNT}"
echo "  Plots  (PNG): ${PLOT_PNG_COUNT}"
echo ""

# 生成 manifest_all.tsv
MANIFEST_PATH="${META_DIR}/manifest_all.tsv"
mkdir -p "${META_DIR}"

{
    printf 'figure_id\tfigure_kind\tsource_id\ttable_path\tplot_pdf\tplot_png\n'
    for table in "${TABLES_DIR}"/*.tsv; do
        [ -f "${table}" ] || continue
        source_id="$(basename "${table}" .tsv)"
        plot_pdf="${PLOTS_DIR}/${source_id}.pdf"
        plot_png="${PLOTS_DIR}/${source_id}.png"

        pdf_exists=""
        png_exists=""
        [ -f "${plot_pdf}" ] && pdf_exists="${plot_pdf}"
        [ -f "${plot_png}" ] && png_exists="${plot_png}"

        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${source_id}" \
            "gwas_manhattan" \
            "${source_id}" \
            "${table}" \
            "${pdf_exists}" \
            "${png_exists}"
    done
} > "${MANIFEST_PATH}"

# 统计缺失
MISSING_PDF=0
MISSING_PNG=0
for table in "${TABLES_DIR}"/*.tsv; do
    [ -f "${table}" ] || continue
    source_id="$(basename "${table}" .tsv)"
    [ -f "${PLOTS_DIR}/${source_id}.pdf" ] || ((MISSING_PDF++)) || true
    [ -f "${PLOTS_DIR}/${source_id}.png" ] || ((MISSING_PNG++)) || true
done

echo "============================================"
echo "  完成"
echo "============================================"
echo "  Manifest:  ${MANIFEST_PATH} (${TABLE_COUNT} 条记录)"
echo "  Tables:    ${TABLE_COUNT}"
echo "  PDF:       ${PLOT_PDF_COUNT}"
echo "  PNG:       ${PLOT_PNG_COUNT}"

if [ "${MISSING_PDF}" -gt 0 ] || [ "${MISSING_PNG}" -gt 0 ]; then
    echo ""
    echo "  [WARN] 不完整: 缺 ${MISSING_PDF} PDF / ${MISSING_PNG} PNG"
    echo "  运行 rerun_failures.sh 重跑失败的 chunk"
fi
echo "============================================"
