#!/usr/bin/env bash
# 生成 cNMF Program Top 基因条形图（不依赖 trait）
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="cnmf_program_top_genes"
WORKFLOW="figures"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${TOPG_MEM:-8G}"; CPUS="${TOPG_CPUS:-2}"
PARTITION="${TOPG_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${TOPG_TIME:-${DEFAULT_TIME}}"
K="${TOPG_K:-60}"; TOP_N="${TOPG_TOP_N:-20}"; PROGRAMS="${TOPG_PROGRAMS:-}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
cfg="${BATCH_DIR}/pipeline.${FIG_NAME}.yaml"; spath="${BATCH_DIR}/figures_${FIG_NAME}.sh"
echo "[${FIG_NAME}] $(date)"

write_config_header "$cfg"
if [[ -n "$PROGRAMS" ]]; then
    cat >> "$cfg" <<YEOF
    cnmf_program_top_genes:
      parameters:
        k: ${K}
        top_n: ${TOP_N}
        programs: [$(yaml_quote "$PROGRAMS")]
YEOF
else
    cat >> "$cfg" <<YEOF
    cnmf_program_top_genes:
      parameters:
        k: ${K}
        top_n: ${TOP_N}
YEOF
fi

write_slurm_script "$spath" "cnmf_program_top_genes" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "sbatch \"${spath}\""
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
