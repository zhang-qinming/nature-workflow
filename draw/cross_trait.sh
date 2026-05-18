#!/usr/bin/env bash
# 生成 Cross-Trait 散点图（成对，需手动指定 pairs）
# 用法: CROSS_TRAIT_PAIRS="GCST90081631:GCST90084301,GCST90081631:GCST90084302" bash gen_cross_trait.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="cross_trait"
WORKFLOW="figures"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${XTRAIT_MEM:-8G}"; CPUS="${XTRAIT_CPUS:-2}"
PARTITION="${XTRAIT_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${XTRAIT_TIME:-${DEFAULT_TIME}}"
TOP_N="${XTRAIT_TOP_N:-12}"

PAIRS="${CROSS_TRAIT_PAIRS:-}"
[[ -n "$PAIRS" ]] || { echo "Usage: CROSS_TRAIT_PAIRS='id1:id2,id3:id4' bash $0" >&2; exit 1; }

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
echo "[${FIG_NAME}] $(date)  pairs=${PAIRS}"

declare -a submit_lines=()
IFS=',' read -r -a pairs_arr <<< "$PAIRS"
for raw in "${pairs_arr[@]}"; do
    raw="$(trim "$raw")"; [[ -z "$raw" ]] && continue
    IFS=':' read -r x y <<< "$raw"; x="$(trim "${x:-}")"; y="$(trim "${y:-}")"
    [[ -n "$x" && -n "$y" ]] || continue
    pid="$(safe_name "${x}__${y}")"
    id_dir="${BATCH_DIR}/${pid}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/figures_${FIG_NAME}_${pid}.sh"

    write_config_header "$cfg"
    cat >> "$cfg" <<YEOF
    cross_trait:
      parameters:
        top_n_labels: ${TOP_N}
        lof_pairs:
          - id_x: $(yaml_quote "$x")
            id_y: $(yaml_quote "$y")
            output_id: $(yaml_quote "$pid")
YEOF

    write_slurm_script "$spath" "cross_trait_${pid}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done

write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
