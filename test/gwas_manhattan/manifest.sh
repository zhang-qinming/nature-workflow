#!/usr/bin/env bash
# ============================================================================
# 构建最终合并 Manifest (全部完成后运行)
# ============================================================================
set -euo pipefail

OUTPUT_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/gwas_manhattan"
TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"
STATUS_DIR="/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/status"

echo "============================================"
echo "  GWAS Manhattan - 最终 Manifest"
echo "============================================"

# 状态总览
ok=$(ls -1 "${STATUS_DIR}"/*.ok 2>/dev/null | wc -l)
failed=$(ls -1 "${STATUS_DIR}"/*.failed 2>/dev/null | wc -l)
echo "  状态: ${ok} 成功 / ${failed} 失败"
echo ""

table_n=$(ls -1 "${TABLES_DIR}"/*.tsv 2>/dev/null | wc -l)
pdf_n=$(ls -1 "${PLOTS_DIR}"/*.pdf 2>/dev/null | wc -l)
png_n=$(ls -1 "${PLOTS_DIR}"/*.png 2>/dev/null | wc -l)
echo "  Tables: ${table_n}  |  PDF: ${pdf_n}  |  PNG: ${png_n}"
echo ""

# 构建 manifest
mkdir -p "${META_DIR}"
MANIFEST="${META_DIR}/manifest_all.tsv"

{
    printf 'figure_id\tfigure_kind\tsource_id\ttable_path\tplot_pdf\tplot_png\n'
    for t in "${TABLES_DIR}"/*.tsv; do
        [ -f "$t" ] || continue
        sid=$(basename "$t" .tsv)
        pdf="${PLOTS_DIR}/${sid}.pdf";  [ -f "$pdf" ] || pdf=""
        png="${PLOTS_DIR}/${sid}.png";  [ -f "$png" ] || png=""
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sid" "gwas_manhattan" "$sid" "$t" "$pdf" "$png"
    done
} > "${MANIFEST}"

# 缺失检查
missing_pdf=0; missing_png=0
for t in "${TABLES_DIR}"/*.tsv; do
    [ -f "$t" ] || continue
    sid=$(basename "$t" .tsv)
    [ -f "${PLOTS_DIR}/${sid}.pdf" ] || ((missing_pdf++)) || true
    [ -f "${PLOTS_DIR}/${sid}.png" ] || ((missing_png++)) || true
done

echo "  Manifest: ${MANIFEST}"
if [ "${missing_pdf}" -gt 0 ] || [ "${missing_png}" -gt 0 ]; then
    echo "  [WARN] 不完整: 缺 ${missing_pdf} PDF / ${missing_png} PNG"
fi

if [ "${failed}" -gt 0 ]; then
    echo ""
    echo "  仍有 ${failed} 个失败，重跑:"
    echo "    bash test/gwas_manhattan/rerun.sh"
fi
echo "============================================"
