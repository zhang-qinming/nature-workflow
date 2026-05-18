#!/usr/bin/env bash
# 一键运行所有图生成器
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Generating all figure scripts ==="
echo ""

for gen in \
    cnmf.sh \
    burden_volcano.sh \
    gene_level_scatter.sh \
    program_rankings.sh \
    gwas_manhattan.sh \
    program_heatmap.sh \
    gene_level_qq.sh \
    cross_trait_heatmap.sh \
    cnmf_program_top_genes.sh \
    cnmf_program_enrichment.sh \
; do
    echo "--- $(date)  Running ${gen} ---"
    bash "${SCRIPT_DIR}/${gen}"
    echo ""
done

echo "=== $(date)  All done ==="
echo "Submit scripts are under: ${SCRIPT_DIR}/../scripts/"
ls "${BATCH_ROOT}/submit_"*.sh
