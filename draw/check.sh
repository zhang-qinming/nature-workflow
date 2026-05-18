#!/usr/bin/env bash
# 检查所有 figure 的输入文件是否存在
# 用法: bash draw/check.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

RED='\033[31m'; GREEN='\033[32m'; YELLOW='\033[33m'; NC='\033[0m'

check() {
    local label="$1" path="$2" cat="$3"
    if [[ -e "$path" ]]; then
        printf "  ${GREEN}✓${NC} %-35s %s\n" "$label" "[${cat}]"
        return 0
    else
        printf "  ${RED}✗${NC} %-35s %s  ${RED}MISSING${NC}\n" "$label" "[${cat}]"
        return 1
    fi
}

echo ""
echo "============================="
echo " Figure 输入文件检查"
echo "============================="
errors=0

echo ""
echo "--- 静态参考数据 (Nature/mine/code/data/) ---"
check "gene_map"       "$GENE_MAP"       "static" || ((errors++))
check "geneset_dir"    "$GENESET_DIR"    "static" || ((errors++))
check "gene_annotation" "$GENE_ANNOTATION" "static" || ((errors++))
check "shet_10bins"    "$SHET_PATH"      "static" || ((errors++))
check "file_id_map"    "$FILE_ID_MAP"    "static" || ((errors++))

echo ""
echo "--- 跑一次的数据 (Nature/mine/code/outputs/) ---"
check "cnmf_regulation" "$CNMF_REGULATION_DIR" "run-once" || ((errors++))
check "limma_logFC_sum" "$LIMMA_PATH"          "run-once" || ((errors++))

echo ""
echo "--- 分 trait 的数据 (run_all/outputs/) ---"
check "posterior_dir"      "$POSTERIOR_DIR"     "per-trait" || ((errors++))
check "program_assoc_dir"  "$PROGRAM_ASSOC_DIR" "per-trait" || ((errors++))

echo ""
if (( errors > 0 )); then
    echo -e "${RED}${errors} 个输入缺失${NC}"
    echo "用 env 覆盖路径，例如: LIMMA_PATH=/other/path bash check.sh"
    exit 1
else
    echo -e "${GREEN}全部输入就绪${NC}"
fi
