#!/usr/bin/env bash
# 生成 cNMF 图：program-regulator 散点 + corregulation 散点
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

# === 本图专属配置 ===
FIG_NAME="cnmf"
WORKFLOW="figures-cnmf"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${CNMF_MEM:-8G}"; CPUS="${CNMF_CPUS:-2}"
PARTITION="${CNMF_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${CNMF_TIME:-${DEFAULT_TIME}}"
K="${CNMF_K:-60}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
_build_gs_yaml

# === 收集 LoF IDs ===
lof_ids=()
while IFS= read -r id; do lof_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$FILE_ID_MAP" | sort -u)
echo "[${FIG_NAME}] $(date)  ${#lof_ids[@]} LoF IDs"

declare -a submit_lines=()
for lof_id in "${lof_ids[@]}"; do
    trait_file="$(resolve_trait_file "$lof_id")"
    id_dir="${BATCH_DIR}/${lof_id}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/cnmf_${lof_id}.sh"

    write_config_header "$cfg"
    cat >> "$cfg" <<YEOF
    cnmf:
      inputs:
        program_association_dir: $(yaml_quote "$PROGRAM_ASSOC_DIR")
        cnmf_regulation_dir: $(yaml_quote "$CNMF_REGULATION_DIR")
      parameters:
        k: ${K}
        trait_targets:
          - trait_file: $(yaml_quote "$trait_file")
            trait_id: $(yaml_quote "$lof_id")
YEOF

    write_slurm_script "$spath" "cnmf_${lof_id}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done

write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
echo "  submit: ${BATCH_ROOT}/submit_${FIG_NAME}.sh"
