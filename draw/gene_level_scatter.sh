#!/usr/bin/env bash
# 生成 Gene-Level 散点图：posterior vs perturb-seq correlation
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="gene_level_scatter"
WORKFLOW="figures-gene-level-scatter"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${GLSCATTER_MEM:-8G}"; CPUS="${GLSCATTER_CPUS:-2}"
PARTITION="${GLSCATTER_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${GLSCATTER_TIME:-${DEFAULT_TIME}}"
TOP_N="${GLSCATTER_TOP_N:-8}"; Y_LIMIT="${GLSCATTER_Y_LIMIT:-8}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*

# 此图不需要基因集
lof_ids=()
while IFS= read -r id; do lof_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$FILE_ID_MAP" | sort -u)
echo "[${FIG_NAME}] $(date)  ${#lof_ids[@]} LoF IDs"

declare -a submit_lines=()
for lof_id in "${lof_ids[@]}"; do
    id_dir="${BATCH_DIR}/${lof_id}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/gene_level_scatter_${lof_id}.sh"

    write_config_header "$cfg"
    cat >> "$cfg" <<YEOF
    gene_level_scatter:
      inputs:
        limma_path: $(yaml_quote "$LIMMA_PATH")
        shet_path: $(yaml_quote "$SHET_PATH")
      parameters:
        lof_ids: [$(yaml_quote "$lof_id")]
        top_n_labels: ${TOP_N}
        y_limit: ${Y_LIMIT}
YEOF

    write_slurm_script "$spath" "gene_level_scatter_${lof_id}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done

write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
