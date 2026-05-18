#!/usr/bin/env bash
# 生成 Gene-Level QQ 图（per trait，自己算 correlation）
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="gene_level_qq"
WORKFLOW="figures-gene-level-qq"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
MEM="${GLQQ_MEM:-32G}"; CPUS="${GLQQ_CPUS:-4}"
PARTITION="${GLQQ_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${GLQQ_TIME:-24:00:00}"
Y_LIMIT="${GLQQ_Y_LIMIT:-10}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*

lof_ids=()
while IFS= read -r id; do lof_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$FILE_ID_MAP" | sort -u)
echo "[${FIG_NAME}] $(date)  ${#lof_ids[@]} LoF IDs"

declare -a submit_lines=()
for lof_id in "${lof_ids[@]}"; do
    id_dir="${BATCH_DIR}/${lof_id}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/gene_level_qq_${lof_id}.sh"

    write_config_header "$cfg"
    cat >> "$cfg" <<YEOF
    gene_level_qq:
      inputs:
        limma_path: $(yaml_quote "$LIMMA_PATH")
        shet_path: $(yaml_quote "$SHET_PATH")
      parameters:
        lof_ids: [$(yaml_quote "$lof_id")]
        y_limit: ${Y_LIMIT}
YEOF

    write_slurm_script "$spath" "gene_level_qq_${lof_id}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done

write_submit "${SCRIPT_DIR}/../scripts/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
