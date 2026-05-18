#!/usr/bin/env bash
# 生成 Cross-Trait 相关性热图（汇总图，所有 trait 一起）
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="cross_trait_heatmap"
WORKFLOW="figures"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${CTH_MEM:-32G}"; CPUS="${CTH_CPUS:-4}"
PARTITION="${CTH_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${CTH_TIME:-24:00:00}"
METHOD="${CTH_METHOD:-pearson}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
cfg="${BATCH_DIR}/pipeline.${FIG_NAME}.yaml"; spath="${BATCH_DIR}/figures_${FIG_NAME}.sh"

lof_ids=()
while IFS= read -r id; do lof_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $2!=""{print $2}' "$FILE_ID_MAP" | sort -u)
(( ${#lof_ids[@]} > 1 )) || { echo "Need >1 LoF IDs, got ${#lof_ids[@]}" >&2; exit 1; }
echo "[${FIG_NAME}] $(date)  ${#lof_ids[@]} LoF IDs"

out_id="$(batch_id "$FIG_NAME" "${lof_ids[@]}")"

write_config_header "$cfg"
cat >> "$cfg" <<YEOF
    cross_trait_heatmap:
      parameters:
        output_id: $(yaml_quote "$out_id")
        method: $(yaml_quote "$METHOD")
        lof_ids:
YEOF
for lof_id in "${lof_ids[@]}"; do
    echo "          - $(yaml_quote "$lof_id")" >> "$cfg"
done

write_slurm_script "$spath" "cross_trait_heatmap" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "sbatch \"${spath}\""
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
