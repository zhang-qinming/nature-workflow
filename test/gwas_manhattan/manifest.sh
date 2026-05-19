#!/usr/bin/env bash
# ============================================================================
# 构建最终合并 Manifest (全部完成后运行)
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/gwas_manhattan}"
TABLES_DIR="${OUTPUT_DIR}/tables"
PLOTS_DIR="${OUTPUT_DIR}/plots"
META_DIR="${OUTPUT_DIR}/meta"
STATUS_DIR="${STATUS_DIR:-${OUTPUT_DIR}/status}"

mapfile -t ALL_IDS < <(tail -n +2 "${FILE_ID_MAP}" | cut -f1)
EXPECTED_TOTAL="${#ALL_IDS[@]}"

ok=0
failed=0
for source_id in "${ALL_IDS[@]}"; do
    [ -f "${STATUS_DIR}/${source_id}.ok" ] && ((ok++)) || true
    [ -f "${STATUS_DIR}/${source_id}.failed" ] && ((failed++)) || true
done

echo "============================================"
echo "  GWAS Manhattan - 最终 Manifest"
echo "============================================"

# 状态总览
echo "  状态: ${ok} 成功 / ${failed} 失败 / ${EXPECTED_TOTAL} 总数"
echo ""

table_n=0
hits_n=0
pdf_n=0
png_n=0
for source_id in "${ALL_IDS[@]}"; do
    [ -f "${TABLES_DIR}/${source_id}.tsv" ] && ((table_n++)) || true
    [ -f "${TABLES_DIR}/${source_id}_hits.tsv" ] && ((hits_n++)) || true
    [ -f "${PLOTS_DIR}/${source_id}.pdf" ] && ((pdf_n++)) || true
    [ -f "${PLOTS_DIR}/${source_id}.png" ] && ((png_n++)) || true
done
echo "  Tables: ${table_n}  |  Hits: ${hits_n}  |  PDF: ${pdf_n}  |  PNG: ${png_n}"
echo ""

# 构建 manifest
mkdir -p "${META_DIR}"
MANIFEST="${META_DIR}/manifest_all.tsv"

{
    printf 'figure_id\tfigure_kind\tsource_id\ttable_path\tplot_pdf\tplot_png\n'
    for sid in "${ALL_IDS[@]}"; do
        t="${TABLES_DIR}/${sid}.tsv"
        [ -f "$t" ] || continue
        pdf="${PLOTS_DIR}/${sid}.pdf"; [ -f "$pdf" ] || pdf=""
        png="${PLOTS_DIR}/${sid}.png"; [ -f "$png" ] || png=""
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sid" "gwas_manhattan" "$sid" "$t" "$pdf" "$png"
    done
} > "${MANIFEST}"

# 缺失检查
missing_table=0
missing_pdf=0; missing_png=0
for sid in "${ALL_IDS[@]}"; do
    [ -f "${TABLES_DIR}/${sid}.tsv" ] || ((missing_table++)) || true
    [ -f "${PLOTS_DIR}/${sid}.pdf" ] || ((missing_pdf++)) || true
    [ -f "${PLOTS_DIR}/${sid}.png" ] || ((missing_png++)) || true
done

echo "  Manifest: ${MANIFEST}"
if [ "${missing_table}" -gt 0 ] || [ "${missing_pdf}" -gt 0 ] || [ "${missing_png}" -gt 0 ]; then
    echo "  [WARN] 不完整: 缺 ${missing_table} TSV / ${missing_pdf} PDF / ${missing_png} PNG"
fi

if [ "${failed}" -gt 0 ]; then
    echo ""
    echo "  仍有 ${failed} 个失败，重跑:"
    echo "    bash test/gwas_manhattan/rerun.sh"
fi
echo "============================================"
