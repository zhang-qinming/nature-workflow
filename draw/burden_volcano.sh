#!/usr/bin/env bash
# 生成 LoF Burden 火山图
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="burden_volcano"
WORKFLOW="figures-burden-volcano"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${BURDEN_VOLCANO_MEM:-8G}"; CPUS="${BURDEN_VOLCANO_CPUS:-2}"
PARTITION="${BURDEN_VOLCANO_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${BURDEN_VOLCANO_TIME:-${DEFAULT_TIME}}"
LABEL_FDR="${BURDEN_VOLCANO_LABEL_FDR:-0.01}"; LINE_FDR="${BURDEN_VOLCANO_LINE_FDR:-0.1}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
echo "[${FIG_NAME}] $(date)  checking inputs..."
check_inputs "$FIG_NAME" \
  "geneset_dir     " "$GENESET_DIR"  "static    " yes \
  "gene_map        " "$GENE_MAP"     "static    " yes \
  "file_id_map     " "$FILE_ID_MAP"  "static    " yes \
  "posterior_dir   " "$POSTERIOR_DIR" "per-trait" yes \
  || exit 1
_build_gs_yaml "$HIGHLIGHT_GENESETS" HL_GS_YAML
_build_gs_yaml "$GENESET_LIST" GS_YAML

lof_ids=()
while IFS= read -r id; do lof_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$FILE_ID_MAP" | sort -u)
echo "[${FIG_NAME}] $(date)  ${#lof_ids[@]} LoF IDs"

declare -a submit_lines=()
for lof_id in "${lof_ids[@]}"; do
    id_dir="${BATCH_DIR}/${lof_id}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/burden_volcano_${lof_id}.sh"

    write_config_header "$cfg"
    cat >> "$cfg" <<YEOF
    burden_volcano:
      parameters:
        lof_ids: [$(yaml_quote "$lof_id")]
        label_fdr_threshold: ${LABEL_FDR}
        line_fdr_threshold: ${LINE_FDR}
        highlight_genesets:
${HL_GS_YAML}        data_genesets:
${GS_YAML}
YEOF

    write_slurm_script "$spath" "burden_volcano_${lof_id}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done

write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
