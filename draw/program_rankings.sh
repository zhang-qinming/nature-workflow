#!/usr/bin/env bash
# 生成 Program 排名条形图
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="program_rankings"
WORKFLOW="figures-program-rankings"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${PRANK_MEM:-4G}"; CPUS="${PRANK_CPUS:-1}"
PARTITION="${PRANK_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${PRANK_TIME:-${DEFAULT_TIME}}"
K="${PRANK_K:-60}"; TOP_N="${PRANK_TOP_N:-12}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*

lof_ids=()
while IFS= read -r id; do lof_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$FILE_ID_MAP" | sort -u)
echo "[${FIG_NAME}] $(date)  ${#lof_ids[@]} LoF IDs"

declare -a submit_lines=()
for lof_id in "${lof_ids[@]}"; do
    trait_file="$(resolve_trait_file "$lof_id")"
    id_dir="${BATCH_DIR}/${lof_id}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/program_rankings_${lof_id}.sh"

    write_config_header "$cfg"
    cat >> "$cfg" <<YEOF
    program_rankings:
      inputs:
        program_association_dir: $(yaml_quote "$PROGRAM_ASSOC_DIR")
      parameters:
        k: ${K}
        top_n: ${TOP_N}
        trait_targets:
          - trait_file: $(yaml_quote "$trait_file")
            trait_id: $(yaml_quote "$lof_id")
YEOF

    write_slurm_script "$spath" "program_rankings_${lof_id}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done

write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
