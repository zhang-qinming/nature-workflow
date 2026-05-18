#!/usr/bin/env bash
# 生成 GWAS Manhattan 全基因组图
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="gwas_manhattan"
WORKFLOW="figures-gwas-manhattan"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${GMAN_MEM:-64G}"; CPUS="${GMAN_CPUS:-4}"
PARTITION="${GMAN_PARTITION:-fat}"; TIME="${GMAN_TIME:-48:00:00}"
FLANK_BP="${GMAN_FLANK_BP:-50000}"; LABEL_P="${GMAN_LABEL_P:-1e-30}"
GENOMEWIDE_P="${GMAN_GENOMEWIDE_P:-5e-8}"

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
_build_gs_yaml

gwas_ids=()
while IFS= read -r id; do gwas_ids+=("$id"); done < <(awk -F'\t' 'NR>1 && $1!=""{print $1}' "$FILE_ID_MAP" | sort -u)
echo "[${FIG_NAME}] $(date)  ${#gwas_ids[@]} GWAS IDs"

declare -a submit_lines=()
for gwas_id in "${gwas_ids[@]}"; do
    id_dir="${BATCH_DIR}/${gwas_id}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/gwas_manhattan_${gwas_id}.sh"

    write_config_header "$cfg"
    cat >> "$cfg" <<YEOF
    gwas_manhattan:
      parameters:
        gwas_ids: [$(yaml_quote "$gwas_id")]
        flank_bp: ${FLANK_BP}
        label_p_threshold: ${LABEL_P}
        genomewide_threshold: ${GENOMEWIDE_P}
        highlight_genesets:
${GS_YAML}
YEOF

    write_slurm_script "$spath" "gwas_manhattan_${gwas_id}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done

write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
