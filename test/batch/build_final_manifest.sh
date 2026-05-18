#!/usr/bin/env bash
# ============================================================================
# 全部 chunk 完成后，构建最终的合并 manifest
#
# 用法:
#   bash test/batch/build_final_manifest.sh
#
# 输出:
#   /gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/gwas_manhattan/meta/manifest_all.tsv
# ============================================================================

set -euo pipefail

OUTPUT_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/gwas_manhattan"
TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"

echo "============================================"
echo "  构建最终合并 Manifest"
echo "============================================"

if [ ! -d "${TABLES_DIR}" ]; then
    echo "[FAIL] tables 目录不存在: ${TABLES_DIR}"
    echo "请等待至少一个 chunk 完成后再运行。"
    exit 1
fi

# 统计输出文件
TABLE_COUNT=$(ls -1 "${TABLES_DIR}"/*.tsv 2>/dev/null | wc -l)
PLOT_PDF_COUNT=$(ls -1 "${PLOTS_DIR}"/*.pdf 2>/dev/null | wc -l)
PLOT_PNG_COUNT=$(ls -1 "${PLOTS_DIR}"/*.png 2>/dev/null | wc -l)

echo "  Tables (TSV): ${TABLE_COUNT}"
echo "  Plots  (PDF): ${PLOT_PDF_COUNT}"
echo "  Plots  (PNG): ${PLOT_PNG_COUNT}"

# 生成 manifest_all.tsv
MANIFEST_PATH="${META_DIR}/manifest_all.tsv"
mkdir -p "${META_DIR}"

echo "  生成: ${MANIFEST_PATH}"

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

echo ""
echo "============================================"
echo "  完成!"
echo "============================================"
echo "  Manifest:  ${MANIFEST_PATH}"
echo "  Tables:    ${TABLE_COUNT} 个"
echo "  PDF:       ${PLOT_PDF_COUNT} 个"
echo "  PNG:       ${PLOT_PNG_COUNT} 个"
echo ""

# 统计缺失的 (有 table 但没有 plot 的)
MISSING_PDF=0
MISSING_PNG=0
for table in "${TABLES_DIR}"/*.tsv; do
    [ -f "${table}" ] || continue
    source_id="$(basename "${table}" .tsv)"
    [ -f "${PLOTS_DIR}/${source_id}.pdf" ] || ((MISSING_PDF++)) || true
    [ -f "${PLOTS_DIR}/${source_id}.png" ] || ((MISSING_PNG++)) || true
done

if [ "${MISSING_PDF}" -gt 0 ] || [ "${MISSING_PNG}" -gt 0 ]; then
    echo "[WARN] 存在不完整的输出:"
    echo "  缺少 PDF: ${MISSING_PDF} 个"
    echo "  缺少 PNG: ${MISSING_PNG} 个"
    echo "  请用 rerun_failures.sh 重跑失败的 chunk"
fi
