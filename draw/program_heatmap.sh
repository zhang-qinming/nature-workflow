#!/usr/bin/env bash
# 生成 Program × Trait 热图（汇总图，所有 trait 一起）
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="program_heatmap"
WORKFLOW="figures"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${PHEAT_MEM:-32G}"; CPUS="${PHEAT_CPUS:-4}"
PARTITION="${PHEAT_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${PHEAT_TIME:-24:00:00}"
K="${PHEAT_K:-60}"; METRICS="${PHEAT_METRICS:-program_score,regulator_score}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
cfg="${BATCH_DIR}/pipeline.${FIG_NAME}.yaml"; spath="${BATCH_DIR}/figures_${FIG_NAME}.sh"

# 收集所有 LoF IDs
lof_ids=()
while IFS= read -r id; do lof_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$FILE_ID_MAP" | sort -u)
echo "[${FIG_NAME}] $(date)  ${#lof_ids[@]} LoF IDs"

out_id="$(batch_id "$FIG_NAME" "${lof_ids[@]}")"

write_config_header "$cfg"
cat >> "$cfg" <<YEOF
    program_heatmap:
      inputs:
        program_association_dir: $(yaml_quote "$PROGRAM_ASSOC_DIR")
      parameters:
        k: ${K}
        output_id: $(yaml_quote "$out_id")
        metrics: [$(yaml_quote "$METRICS")]
        trait_targets:
YEOF
for lof_id in "${lof_ids[@]}"; do
    tf="$(resolve_trait_file "$lof_id")"
    cat >> "$cfg" <<YEOF
          - trait_file: $(yaml_quote "$tf")
            trait_id: $(yaml_quote "$lof_id")
YEOF
done

write_slurm_script "$spath" "program_heatmap" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "sbatch \"${spath}\""
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
