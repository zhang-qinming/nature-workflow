#!/usr/bin/env bash
# 生成 cNMF Program 基因集富集气泡图（不依赖 trait）
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="cnmf_program_enrichment"
WORKFLOW="figures"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${ENRICH_MEM:-8G}"; CPUS="${ENRICH_CPUS:-2}"
PARTITION="${ENRICH_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${ENRICH_TIME:-${DEFAULT_TIME}}"
K="${ENRICH_K:-60}"; TOP_N="${ENRICH_TOP_N:-100}"; PROGRAMS="${ENRICH_PROGRAMS:-}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
cfg="${BATCH_DIR}/pipeline.${FIG_NAME}.yaml"; spath="${BATCH_DIR}/figures_${FIG_NAME}.sh"
echo "[${FIG_NAME}] $(date)"

_build_gs_yaml

write_config_header "$cfg"
if [[ -n "$PROGRAMS" ]]; then
    cat >> "$cfg" <<YEOF
    cnmf_program_enrichment:
      inputs:
        spectra_path: $(yaml_quote "$SPECTRA_PATH")
      parameters:
        k: ${K}
        top_n: ${TOP_N}
        programs: [$(yaml_quote "$PROGRAMS")]
        genesets:
${GS_YAML}
YEOF
else
    cat >> "$cfg" <<YEOF
    cnmf_program_enrichment:
      inputs:
        spectra_path: $(yaml_quote "$SPECTRA_PATH")
      parameters:
        k: ${K}
        top_n: ${TOP_N}
        genesets:
${GS_YAML}
YEOF
fi

write_slurm_script "$spath" "enrich" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "sbatch \"${spath}\""
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
