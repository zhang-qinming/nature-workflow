#!/usr/bin/env bash
# 生成 GWAS Locus Zoom 图（需手动指定坐标 TSV）
# 用法: LOCI_TSV=/path/to/loci.tsv bash gen_gwas_locus_zoom.sh
# TSV 格式 (8列): gwas_id  gwas_file  chrom  start  end  locus_label  flank_bp  output_id
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

FIG_NAME="gwas_locus_zoom"
WORKFLOW="figures"
BATCH_DIR="${BATCH_ROOT}/${FIG_NAME}"
LOG_DIR="${LOGS_ROOT}/${FIG_NAME}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${DEFAULT_PROJECT_ROOT}/outputs}"
MEM="${LZM_MEM:-8G}"; CPUS="${LZM_CPUS:-2}"
PARTITION="${LZM_PARTITION:-${DEFAULT_PARTITION}}"; TIME="${LZM_TIME:-${DEFAULT_TIME}}"
FLANK_BP_DEF="${LZM_FLANK_BP:-250000}"; LABEL_TOP="${LZM_LABEL_TOP:-6}"
GENOMEWIDE_P="${LZM_GENOMEWIDE_P:-5e-8}"

LOCI_TSV="${LOCI_TSV:-}"
[[ -n "$LOCI_TSV" && -f "$LOCI_TSV" ]] || { echo "Usage: LOCI_TSV=/path/to/loci.tsv bash $0" >&2; exit 1; }

mkdir -p "$BATCH_DIR" "$LOG_DIR"
rm -rf "${BATCH_DIR:?}"/*
echo "[${FIG_NAME}] $(date)"

declare -a submit_lines=()
while IFS=$'\t' read -r gwas_id gwas_file chrom start end locus_label flank_bp output_id; do
    [[ -n "${gwas_id:-}" || -n "${gwas_file:-}" ]] || continue
    chrom="$(trim "${chrom:-}")"; start="$(trim "${start:-}")"; end="$(trim "${end:-}")"
    [[ -n "$chrom" && -n "$start" && -n "$end" ]] || continue
    locus_label="${locus_label:-chr${chrom}:${start}-${end}}"
    flank_bp="${flank_bp:-${FLANK_BP_DEF}}"
    output_id="${output_id:-$(safe_name "${gwas_id:-$gwas_file}")__chr${chrom}_${start}_${end}}"
    output_id="$(safe_name "$output_id")"

    id_dir="${BATCH_DIR}/${output_id}"; mkdir -p "$id_dir"
    cfg="${id_dir}/pipeline.${FIG_NAME}.yaml"; spath="${id_dir}/figures_${FIG_NAME}_${output_id}.sh"

    write_config_header "$cfg"
    if [[ -n "${gwas_id:-}" ]]; then
        cat >> "$cfg" <<YEOF
    gwas_locus_zoom:
      parameters:
        genomewide_threshold: ${GENOMEWIDE_P}
        label_top_n: ${LABEL_TOP}
        locus_targets:
          - gwas_id: $(yaml_quote "$gwas_id")
            chrom: $(yaml_quote "$chrom")
            start: ${start}
            end: ${end}
            flank_bp: ${flank_bp}
            locus_label: $(yaml_quote "$locus_label")
            output_id: $(yaml_quote "$output_id")
YEOF
    else
        cat >> "$cfg" <<YEOF
    gwas_locus_zoom:
      parameters:
        genomewide_threshold: ${GENOMEWIDE_P}
        label_top_n: ${LABEL_TOP}
        locus_targets:
          - gwas_file: $(yaml_quote "$gwas_file")
            chrom: $(yaml_quote "$chrom")
            start: ${start}
            end: ${end}
            flank_bp: ${flank_bp}
            locus_label: $(yaml_quote "$locus_label")
            output_id: $(yaml_quote "$output_id")
YEOF
    fi

    write_slurm_script "$spath" "gwas_locus_zoom_${output_id}" "$MEM" "$CPUS" "$PARTITION" "$TIME" "$WORKFLOW" "$cfg"
    submit_lines+=("sbatch \"${spath}\"")
done < <(awk -F'\t' 'NR>1' "$LOCI_TSV")

write_submit "${BATCH_ROOT}/submit_${FIG_NAME}.sh" "${submit_lines[@]}"
echo "[${FIG_NAME}] $(date)  Done: ${BATCH_DIR}"
