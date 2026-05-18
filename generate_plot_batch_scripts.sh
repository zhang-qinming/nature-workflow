#!/usr/bin/env bash

set -euo pipefail

DEFAULT_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_ROOT="${BATCH_ROOT:-${DEFAULT_PROJECT_ROOT}/scripts/plot}"
PLOT_ROOT="${PLOT_ROOT:-${BATCH_ROOT}}"
PLOT_TRAITS_ROOT="${PLOT_TRAITS_ROOT:-${PLOT_ROOT}/traits}"
PLOT_SUPPLEMENTARY_ROOT="${PLOT_SUPPLEMENTARY_ROOT:-${PLOT_ROOT}/supplementary}"
PLOT_TOOLS_ROOT="${PLOT_TOOLS_ROOT:-${PLOT_ROOT}/tools}"
LOGS_ROOT="${LOGS_ROOT:-${PLOT_ROOT}/logs}"
PLOT_CONFIG_NAME="${PLOT_CONFIG_NAME:-pipeline.plot.yaml}"

CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
PLOT_ENV="${PLOT_ENV:-paper-pipeline-plot}"
TOOLS_ENV="${TOOLS_ENV:-paper-pipeline-tools}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

USE_TOOLS_ENV_WRAPPERS="${USE_TOOLS_ENV_WRAPPERS:-1}"
LDSC_REPO_DIR="${LDSC_REPO_DIR:-${HOME}/software/ldsc}"
LDSC_PY="${LDSC_PY:-}"
LIFTOVER_BIN="${LIFTOVER_BIN:-}"

PLOT_TRAIT_NUM="${PLOT_TRAIT_NUM:-86}"
PLOT_TRAIT_NUMS="${PLOT_TRAIT_NUMS:-${PLOT_TRAIT_NUM}}"
PLOT_TRAIT_FILE="${PLOT_TRAIT_FILE:-Backman_2021_${PLOT_TRAIT_NUM}.per_gene_estimates.tsv}"
PLOT_RAW_BURDEN_FILE="${PLOT_RAW_BURDEN_FILE:-Backman_2021_${PLOT_TRAIT_NUM}_M1_001.summary_statistics.csv}"
PLOT_TRAIT_LABEL="${PLOT_TRAIT_LABEL:-MCH}"
PLOT_TRAIT_FILE_TEMPLATE="${PLOT_TRAIT_FILE_TEMPLATE:-Backman_2021_{trait}.per_gene_estimates.tsv}"
PLOT_RAW_BURDEN_FILE_TEMPLATE="${PLOT_RAW_BURDEN_FILE_TEMPLATE:-Backman_2021_{trait}_M1_001.summary_statistics.csv}"
PLOT_TRAIT_LABEL_TEMPLATE="${PLOT_TRAIT_LABEL_TEMPLATE:-Backman_2021_{trait}}"
PLOT_GWAS_TRAIT="${PLOT_GWAS_TRAIT:-30050}"
PLOT_CROSS_TRAIT_TRAITS="${PLOT_CROSS_TRAIT_TRAITS:-86, 88}"
PLOT_CNMF_K="${PLOT_CNMF_K:-60}"
PLOT_CNMF_MODE="${PLOT_CNMF_MODE:-generic}"
PLOT_ENABLE_MULTIPLE_REGRESSION="${PLOT_ENABLE_MULTIPLE_REGRESSION:-inherit}"
PLOT_ENABLE_TRANS_EQTL="${PLOT_ENABLE_TRANS_EQTL:-inherit}"
PLOT_LABEL_PROGRAMS="${PLOT_LABEL_PROGRAMS:-4, 16, 25, 40}"
PLOT_CORREGULATION_PROGRAM_A="${PLOT_CORREGULATION_PROGRAM_A:-P25}"
PLOT_CORREGULATION_PROGRAM_B="${PLOT_CORREGULATION_PROGRAM_B:-P16}"
PLOT_TRANS_EQTL_TOP_N="${PLOT_TRANS_EQTL_TOP_N:-100}"
PLOT_TRANS_EQTL_MATCH_TABLE="${PLOT_TRANS_EQTL_MATCH_TABLE:-}"
PLOT_MCH_POSTERIOR_FILE="${PLOT_MCH_POSTERIOR_FILE:-Backman_2021_86.per_gene_estimates.tsv}"
PLOT_RDW_POSTERIOR_FILE="${PLOT_RDW_POSTERIOR_FILE:-Backman_2021_88.per_gene_estimates.tsv}"
PLOT_IRF_POSTERIOR_FILE="${PLOT_IRF_POSTERIOR_FILE:-Backman_2021_106.per_gene_estimates.tsv}"
PLOT_AUTOPHAGY_GENESET_FILE="${PLOT_AUTOPHAGY_GENESET_FILE:-Autophagosome_genes.txt}"
SUPPLEMENTARY_GROWTH_PROGRAMS="${SUPPLEMENTARY_GROWTH_PROGRAMS:-P16}"

PLOT_GENESET_DIR="${PLOT_GENESET_DIR:-data/geneset}"
PLOT_GWAS_DIR="${PLOT_GWAS_DIR:-data/GWAS}"
PLOT_GENE_MAP_PATH="${PLOT_GENE_MAP_PATH:-data/gencode_v41_gname_gid_ALL_sorted_onlyID}"
PLOT_BACKMAN_CORRESP_PATH="${PLOT_BACKMAN_CORRESP_PATH:-data/Backman_Niele_corresp}"
PLOT_RAW_BURDEN_DIR="${PLOT_RAW_BURDEN_DIR:-data/LoF/raw_burden}"
PLOT_LD_REFERENCE_DIR="${PLOT_LD_REFERENCE_DIR:-data/GWAS/reference}"

PLOT_BURDEN_REG_CORRELATION_DIR="${PLOT_BURDEN_REG_CORRELATION_DIR:-data/Perturbseq/trait_association/K562GW/GeneLevel}"
PLOT_GENE_LEVEL_TRAIT_METADATA="${PLOT_GENE_LEVEL_TRAIT_METADATA:-data/Perturbseq/trait_association/K562GW/GeneLevel/Backman_Niele_corresp_all}"

PLOT_SHET_PATH="${PLOT_SHET_PATH:-data/shet_10bins.txt}"
PLOT_METADATA_PATH="${PLOT_METADATA_PATH:-data/Perturbseq/metadata/gwps_metadata.csv}"
PLOT_TRANS_EQTL_PATH="${PLOT_TRANS_EQTL_PATH:-data/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz}"

PLOT_LDSC_DIR="${PLOT_LDSC_DIR:-data/LDSC}"
PLOT_LDSC_SUMSTATS_DIR="${PLOT_LDSC_SUMSTATS_DIR:-data/LDSC/GWAS}"
PLOT_LDSC_LD_REFERENCE_DIR="${PLOT_LDSC_LD_REFERENCE_DIR:-data/LDSC/ref/1000G_EUR_Phase3_plink}"
# Figure1/1_annotation_LDSC.sh reads both:
#   - ${PLOT_LDSC_BASELINE_DIR}/baseline.<chr>.*
#   - ${PLOT_LDSC_BASELINE_DIR}/<chr>.snp
PLOT_LDSC_BASELINE_DIR="${PLOT_LDSC_BASELINE_DIR:-data/LDSC/ref/baseline_v1.2}"
PLOT_LDSC_WEIGHTS_PREFIX="${PLOT_LDSC_WEIGHTS_PREFIX:-data/LDSC/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.}"
PLOT_LDSC_FRQ_PREFIX="${PLOT_LDSC_FRQ_PREFIX:-data/LDSC/ref/1000G_Phase3_frq/1000G.EUR.QC.}"

PLOT_CHIP_DIR="${PLOT_CHIP_DIR:-data/TF_ChIP}"
PLOT_MULTICELL_CNMF_ROOT="${PLOT_MULTICELL_CNMF_ROOT:-data/Perturbseq/cNMF}"
PLOT_MULTICELL_CNMF_REGULATION_DIR="${PLOT_MULTICELL_CNMF_REGULATION_DIR:-data/Perturbseq/cNMF_regulation}"
PLOT_MULTICELL_TRAIT_ASSOCIATION_DIR="${PLOT_MULTICELL_TRAIT_ASSOCIATION_DIR:-data/Perturbseq/trait_association}"
PLOT_GROWTH_SCREENING_CSV="${PLOT_GROWTH_SCREENING_CSV:-data/Growth_screening_2014.csv}"
PLOT_LIFTOVER_CHAIN="${PLOT_LIFTOVER_CHAIN:-data/hg38ToHg19.over.chain.gz}"

PLOT_PARTITION="${PLOT_PARTITION:-cu,privority,batch01}"
PLOT_TIME="${PLOT_TIME:-30-00:00:00}"

PLOT_GENE_BURDEN_MEM="${PLOT_GENE_BURDEN_MEM:-40G}"
PLOT_GENE_BURDEN_CPUS="${PLOT_GENE_BURDEN_CPUS:-2}"

PLOT_GENE_LEVEL_MEM="${PLOT_GENE_LEVEL_MEM:-8G}"
PLOT_GENE_LEVEL_CPUS="${PLOT_GENE_LEVEL_CPUS:-2}"

PLOT_CNMF_MEM="${PLOT_CNMF_MEM:-100G}"
PLOT_CNMF_CPUS="${PLOT_CNMF_CPUS:-8}"
PLOT_CNMF_TRANS_EQTL_PARALLEL_JOBS="${PLOT_CNMF_TRANS_EQTL_PARALLEL_JOBS:-${PLOT_CNMF_CPUS}}"

PLOT_SUPPLEMENTARY_MEM="${PLOT_SUPPLEMENTARY_MEM:-150G}"
PLOT_SUPPLEMENTARY_CPUS="${PLOT_SUPPLEMENTARY_CPUS:-16}"
PLOT_SUPPLEMENTARY_TIME="${PLOT_SUPPLEMENTARY_TIME:-10-00:00:00}"

mkdir -p "${BATCH_ROOT}"
BATCH_ROOT="$(cd "${BATCH_ROOT}" && pwd)"
PLOT_ROOT="${PLOT_ROOT:-${BATCH_ROOT}}"
mkdir -p "${PLOT_ROOT}"
PLOT_ROOT="$(cd "${PLOT_ROOT}" && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-}"
if [[ -z "${PROJECT_ROOT}" ]]; then
    PROJECT_ROOT="$(cd "${PLOT_ROOT}/../.." && pwd)"
fi
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "PROJECT_ROOT does not exist: ${PROJECT_ROOT}" >&2
    exit 1
fi
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${PROJECT_ROOT}/outputs}"
if [[ "${ARTIFACT_ROOT}" != /* ]]; then
    ARTIFACT_ROOT="$(cd "${PROJECT_ROOT}" && realpath -m "${ARTIFACT_ROOT}")"
fi

for _dir in "${PLOT_TRAITS_ROOT}" "${PLOT_SUPPLEMENTARY_ROOT}" "${PLOT_TOOLS_ROOT}"; do
    if [[ -z "${_dir}" || "${_dir}" == "/" ]]; then
        echo "FATAL: refusing to rm -rf empty or root path" >&2; exit 1
    fi
done
rm -rf "${PLOT_TRAITS_ROOT}" "${PLOT_SUPPLEMENTARY_ROOT}" "${PLOT_TOOLS_ROOT}"
mkdir -p "${PLOT_TRAITS_ROOT}" "${PLOT_SUPPLEMENTARY_ROOT}" "${PLOT_TOOLS_ROOT}" "${LOGS_ROOT}"
PLOT_TRAITS_ROOT="$(cd "${PLOT_TRAITS_ROOT}" && pwd)"
PLOT_SUPPLEMENTARY_ROOT="$(cd "${PLOT_SUPPLEMENTARY_ROOT}" && pwd)"
PLOT_TOOLS_ROOT="$(cd "${PLOT_TOOLS_ROOT}" && pwd)"
LOGS_ROOT="$(cd "${LOGS_ROOT}" && pwd)"

write_tool_wrapper() {
    local script_path="$1"
    local tool_bin="$2"

    cat > "$script_path" <<EOF
#!/usr/bin/env bash
set -euo pipefail

source "${CONDA_SH}"
exec conda run --no-capture-output -n "${TOOLS_ENV}" "${tool_bin}" "\$@"
EOF

    chmod +x "$script_path"
}

write_ldsc_wrapper() {
    local script_path="$1"

    cat > "$script_path" <<EOF
#!/usr/bin/env bash
set -euo pipefail

source "${CONDA_SH}"
exec conda run --no-capture-output -n "${TOOLS_ENV}" python "${LDSC_REPO_DIR}/ldsc.py" "\$@"
EOF

    chmod +x "$script_path"
}

if [[ "${USE_TOOLS_ENV_WRAPPERS}" == "1" ]]; then
    write_ldsc_wrapper "${PLOT_TOOLS_ROOT}/ldsc_py_wrapper.sh"
    write_tool_wrapper "${PLOT_TOOLS_ROOT}/liftOver_wrapper.sh" "liftOver"
    LDSC_PY="${LDSC_PY:-${PLOT_TOOLS_ROOT}/ldsc_py_wrapper.sh}"
    LIFTOVER_BIN="${LIFTOVER_BIN:-${PLOT_TOOLS_ROOT}/liftOver_wrapper.sh}"
else
    LDSC_PY="${LDSC_PY:-${LDSC_REPO_DIR}/ldsc.py}"
    LIFTOVER_BIN="${LIFTOVER_BIN:-liftOver}"
fi

replace_template() {
    local template="$1"
    local trait="$2"
    printf '%s' "${template//\{trait\}/${trait}}"
}

trim() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "$value"
}

yaml_bool() {
    local raw
    raw="$(trim "$1")"
    case "${raw,,}" in
        1|true|yes|on) printf 'true' ;;
        0|false|no|off) printf 'false' ;;
        *)
            echo "Boolean parameter must be one of: 0/1, true/false, yes/no, on/off. Got: ${raw}" >&2
            exit 1
            ;;
    esac
}

yaml_mode() {
    local raw
    raw="$(trim "$1")"
    case "${raw,,}" in
        generic|legacy) printf '%s' "${raw,,}" ;;
        *)
            echo "Mode parameter must be one of: generic, legacy. Got: ${raw}" >&2
            exit 1
            ;;
    esac
}

resolve_path() {
    local raw_path="$1"
    if [[ "$raw_path" = /* ]]; then
        printf '%s' "$raw_path"
    else
        (
            cd "${PROJECT_ROOT}"
            realpath -m "$raw_path"
        )
    fi
}

trait_numbers=()
IFS=',' read -r -a _trait_raw <<< "${PLOT_TRAIT_NUMS}"
for item in "${_trait_raw[@]}"; do
    item="$(trim "$item")"
    [[ -z "$item" ]] && continue
    if [[ "$item" =~ [^0-9] ]]; then
        echo "PLOT_TRAIT_NUMS must contain comma-separated integers" >&2
        exit 1
    fi
    trait_numbers+=("$item")
done

if (( ${#trait_numbers[@]} == 0 )); then
    echo "PLOT_TRAIT_NUMS must contain at least one trait" >&2
    exit 1
fi

trait_file_for() {
    local trait="$1"
    if (( ${#trait_numbers[@]} == 1 )) && [[ "${trait}" == "${PLOT_TRAIT_NUM}" ]]; then
        printf '%s' "${PLOT_TRAIT_FILE}"
    else
        replace_template "${PLOT_TRAIT_FILE_TEMPLATE}" "${trait}"
    fi
}

raw_burden_file_for() {
    local trait="$1"
    if (( ${#trait_numbers[@]} == 1 )) && [[ "${trait}" == "${PLOT_TRAIT_NUM}" ]]; then
        printf '%s' "${PLOT_RAW_BURDEN_FILE}"
    else
        replace_template "${PLOT_RAW_BURDEN_FILE_TEMPLATE}" "${trait}"
    fi
}

trait_label_for() {
    local trait="$1"
    if (( ${#trait_numbers[@]} == 1 )) && [[ "${trait}" == "${PLOT_TRAIT_NUM}" ]]; then
        printf '%s' "${PLOT_TRAIT_LABEL}"
    else
        replace_template "${PLOT_TRAIT_LABEL_TEMPLATE}" "${trait}"
    fi
}

lof_trait_id() {
    local trait_file="$1"
    trait_file="${trait_file%.per_gene_estimates.tsv}"
    if [[ "$trait_file" =~ (GCST[0-9]+) ]]; then
        printf '%s' "${BASH_REMATCH[1]}"
    else
        printf '%s' "${trait_file%%_*}"
    fi
}

plot_safe_name() {
    local trait_num="$1"
    local trait_file="$2"
    trait_file="${trait_file%.per_gene_estimates.tsv}"
    if [[ "$trait_file" =~ (GCST[0-9]+) ]]; then
        printf '%s' "${BASH_REMATCH[1]}"
    else
        printf 'trait%s' "$trait_num"
    fi
}

load_trans_eqtl_match_table() {
    local table="$1"
    local id1 id2 path1 path2 rest

    [[ -z "$table" ]] && return 0
    if [[ ! -f "$table" ]]; then
        echo "WARN: trans-eQTL GWAS match table not found: ${table}; legacy-mode auto matching disabled" >&2
        return 0
    fi

    while IFS=$'\t' read -r id1 id2 path1 path2 rest; do
        [[ -z "${id1:-}" || "$id1" == "id1" ]] && continue
        [[ -z "${id2:-}" || -z "${path1:-}" ]] && continue
        TRANS_EQTL_GWAS_ID_BY_LOF["$id2"]="$id1"
        TRANS_EQTL_GWAS_PATH_BY_LOF["$id2"]="$path1"
    done < "$table"
}

multiple_regression_flag() {
    local mode="$1"
    case "${PLOT_ENABLE_MULTIPLE_REGRESSION,,}" in
        inherit|'')
            [[ "$mode" == "legacy" ]] && printf 'true' || printf 'false'
            ;;
        1|true|yes|on)
            printf 'true'
            ;;
        0|false|no|off)
            printf 'false'
            ;;
        *)
            echo "PLOT_ENABLE_MULTIPLE_REGRESSION must be inherit or a boolean value. Got: ${PLOT_ENABLE_MULTIPLE_REGRESSION}" >&2
            exit 1
            ;;
    esac
}

trans_eqtl_run_flag() {
    local mode="$1"
    local matched="$2"
    case "${PLOT_ENABLE_TRANS_EQTL,,}" in
        inherit|auto|'')
            if [[ "$mode" == "legacy" && "$matched" == "1" ]]; then
                printf 'true'
            else
                printf 'false'
            fi
            ;;
        1|true|yes|on)
            printf 'true'
            ;;
        0|false|no|off)
            printf 'false'
            ;;
        *)
            echo "PLOT_ENABLE_TRANS_EQTL must be inherit/auto or a boolean value. Got: ${PLOT_ENABLE_TRANS_EQTL}" >&2
            exit 1
            ;;
    esac
}

write_slurm_script() {
    local script_path="$1"
    local job_name="$2"
    local mem="$3"
    local cpus="$4"
    local partition="$5"
    local time_limit="$6"
    local workflow_name="$7"
    local config_file="${8:-${PLOT_CONFIG_NAME}}"
    local script_dir
    local config_path

    script_dir="$(cd "$(dirname "$script_path")" && pwd)"
    if [[ "$config_file" = /* ]]; then
        config_path="$config_file"
    else
        config_path="${script_dir}/${config_file}"
    fi

    cat > "$script_path" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${job_name}
#SBATCH --error=${LOGS_ROOT}/${job_name}.error
#SBATCH --output=${LOGS_ROOT}/${job_name}.out
#SBATCH --mem=${mem}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --partition=${partition}
#SBATCH --chdir=${script_dir}
#SBATCH --time=${time_limit}

set -euo pipefail

source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"

cd "${script_dir}"
paper-pipeline run --config "${config_path}" ${workflow_name}
EOF

    chmod +x "$script_path"
}

write_trait_plot_config() {
    local config_path="$1"
    local trait_num="$2"
    local trait_file="$3"
    local raw_burden_file="$4"
    local trait_label="$5"
    local cnmf_mode="$6"
    local run_multiple_regression="$7"
    local run_trans_eqtl_follow_up="$8"
    local trans_eqtl_gwas_trait="$9"
    local trans_eqtl_gwas_file="${10}"

    cat > "$config_path" <<EOF
project_root: ${PROJECT_ROOT}
artifact_root: ${ARTIFACT_ROOT}

executables:
  plot_bash:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${PLOT_ENV}
    - bash
  plot_rscript:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${PLOT_ENV}
    - Rscript
  ldsc_py: ${LDSC_PY}
  liftover: ${LIFTOVER_BIN}

workflows:
  plot:
    gene_burden:
      inputs:
        raw_burden_dir: $(resolve_path "${PLOT_RAW_BURDEN_DIR}")
        gwas_dir: $(resolve_path "${PLOT_GWAS_DIR}")
        geneset_dir: $(resolve_path "${PLOT_GENESET_DIR}")
        gene_map: $(resolve_path "${PLOT_GENE_MAP_PATH}")
        trait_correspondence: $(resolve_path "${PLOT_BACKMAN_CORRESP_PATH}")
        ld_reference_dir: $(resolve_path "${PLOT_LD_REFERENCE_DIR}")
      outputs:
        output_dir: ${ARTIFACT_ROOT}/plots/trait${trait_num}/gene_burden
      parameters:
        gwas_traits: [${PLOT_GWAS_TRAIT}]
        manhattan_gwas_trait: ${PLOT_GWAS_TRAIT}
        volcano_burden_file: ${raw_burden_file}
        lof_enrichment_pairs:
          - burden_trait: ${trait_num}
            gwas_trait: ${PLOT_GWAS_TRAIT}
            trait_name: ${trait_label}
        cross_trait_traits: [${PLOT_CROSS_TRAIT_TRAITS}]

    perturbseq_gene_level:
      inputs:
        raw_burden_dir: $(resolve_path "${PLOT_RAW_BURDEN_DIR}")
        gwas_dir: $(resolve_path "${PLOT_GWAS_DIR}")
        gene_map: $(resolve_path "${PLOT_GENE_MAP_PATH}")
        shet_path: $(resolve_path "${PLOT_SHET_PATH}")
        burden_reg_correlation_dir: $(resolve_path "${PLOT_BURDEN_REG_CORRELATION_DIR}")
        gene_level_trait_metadata: $(resolve_path "${PLOT_GENE_LEVEL_TRAIT_METADATA}")
      outputs:
        output_dir: ${ARTIFACT_ROOT}/plots/trait${trait_num}/perturbseq_gene_level
      parameters:
        trait_files: [${trait_file}]
        enrichment_gwas_trait: ${PLOT_GWAS_TRAIT}
        enrichment_comp: closest
        enrichment_target_gene: HBA1
        enrichment_lof_burden_file: ${raw_burden_file}
        render_trait_file: ${trait_file}

    perturbseq_cnmf:
      inputs:
        gwas_dir: $(resolve_path "${PLOT_GWAS_DIR}")
        gene_map: $(resolve_path "${PLOT_GENE_MAP_PATH}")
        shet_path: $(resolve_path "${PLOT_SHET_PATH}")
        metadata: $(resolve_path "${PLOT_METADATA_PATH}")
        trans_eqtl_path: $(resolve_path "${PLOT_TRANS_EQTL_PATH}")
      outputs:
        output_dir: ${ARTIFACT_ROOT}/plots/trait${trait_num}/perturbseq_cnmf
      parameters:
        mode: ${cnmf_mode}
        k: ${PLOT_CNMF_K}
        trait_files: [${trait_file}]
        run_multiple_regression: ${run_multiple_regression}
        run_trans_eqtl_follow_up: ${run_trans_eqtl_follow_up}
        plot_label_programs: [${PLOT_LABEL_PROGRAMS}]
        corregulation_pairs:
          - program_a: ${PLOT_CORREGULATION_PROGRAM_A}
            program_b: ${PLOT_CORREGULATION_PROGRAM_B}
EOF

    if [[ "$run_multiple_regression" == "true" ]]; then
        cat >> "$config_path" <<EOF
        posterior_trait_files:
          MCH: ${PLOT_MCH_POSTERIOR_FILE}
          RDW: ${PLOT_RDW_POSTERIOR_FILE}
          IRF: ${PLOT_IRF_POSTERIOR_FILE}
EOF
    fi

    if [[ "$run_trans_eqtl_follow_up" == "true" ]]; then
        cat >> "$config_path" <<EOF
        trans_eqtl_top_n: ${PLOT_TRANS_EQTL_TOP_N}
        trans_eqtl_parallel_jobs: ${PLOT_CNMF_TRANS_EQTL_PARALLEL_JOBS}
        trans_eqtl_output_label: ${trait_label}
        trans_eqtl_regulator_file: regulators_enrichment_K${PLOT_CNMF_K}_${trait_file}
        trans_eqtl_gwas_trait: ${trans_eqtl_gwas_trait}
        trans_eqtl_gwas_file: ${trans_eqtl_gwas_file}
EOF
    fi
}

write_supplementary_plot_config() {
    local config_path="$1"

    cat > "$config_path" <<EOF
project_root: ${PROJECT_ROOT}
artifact_root: ${ARTIFACT_ROOT}

executables:
  plot_bash:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${PLOT_ENV}
    - bash
  plot_rscript:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${PLOT_ENV}
    - Rscript
  ldsc_py: ${LDSC_PY}
  liftover: ${LIFTOVER_BIN}

workflows:
  plot:
    supplementary:
      inputs:
        ldsc_dir: $(resolve_path "${PLOT_LDSC_DIR}")
        ldsc_ld_reference_dir: $(resolve_path "${PLOT_LDSC_LD_REFERENCE_DIR}")
        ldsc_baseline_dir: $(resolve_path "${PLOT_LDSC_BASELINE_DIR}")
        ldsc_weights_dir: $(resolve_path "${PLOT_LDSC_WEIGHTS_PREFIX}")
        ldsc_frq_dir: $(resolve_path "${PLOT_LDSC_FRQ_PREFIX}")
        ldsc_sumstats_dir: $(resolve_path "${PLOT_LDSC_SUMSTATS_DIR}")
        chip_dir: $(resolve_path "${PLOT_CHIP_DIR}")
        multicell_cnmf_root: $(resolve_path "${PLOT_MULTICELL_CNMF_ROOT}")
        multicell_cnmf_regulation_dir: $(resolve_path "${PLOT_MULTICELL_CNMF_REGULATION_DIR}")
        multicell_trait_association_dir: $(resolve_path "${PLOT_MULTICELL_TRAIT_ASSOCIATION_DIR}")
        gwas_dir: $(resolve_path "${PLOT_GWAS_DIR}")
        growth_screening_csv: $(resolve_path "${PLOT_GROWTH_SCREENING_CSV}")
        gene_map: $(resolve_path "${PLOT_GENE_MAP_PATH}")
        shet_path: $(resolve_path "${PLOT_SHET_PATH}")
        metadata: $(resolve_path "${PLOT_METADATA_PATH}")
        geneset_dir: $(resolve_path "${PLOT_GENESET_DIR}")
        liftOver_chain: $(resolve_path "${PLOT_LIFTOVER_CHAIN}")
      outputs:
        output_dir: ${ARTIFACT_ROOT}/plots/supplementary
      parameters:
        k: ${PLOT_CNMF_K}
        ldsc_annotation_beds: []
        ldsc_sumstats: []
        genetic_correlation_pairs: []
        run_chip_collection: false
        growth_screening_programs: [${SUPPLEMENTARY_GROWTH_PROGRAMS}]
        posterior_trait_files:
          MCH: ${PLOT_MCH_POSTERIOR_FILE}
          RDW: ${PLOT_RDW_POSTERIOR_FILE}
          IRF: ${PLOT_IRF_POSTERIOR_FILE}
        autophagy_geneset_file: ${PLOT_AUTOPHAGY_GENESET_FILE}
        figure5_permutation_jobs:
          - program_n: 5
            regulator_n: 3
            program_top_def: 200
            lof_thresh: 0.1
            trait: MCH
        leave_one_out_jobs:
          - program_n: 5
            regulator_n: 3
            program_top_def: 200
            trait: MCH
        effect_direction_traits: [MCH, RDW, IRF]
        effect_direction_program_top_defs: [100, 200]
EOF
}

write_trait_submitter() {
    local script_path="$1"

    cat > "$script_path" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

gene_burden_jobid="$(sbatch --parsable plot_gene_burden.sh)"
echo "plot-gene-burden jobid: ${gene_burden_jobid}"

gene_level_jobid="$(sbatch --parsable plot_perturbseq_gene_level.sh)"
echo "plot-perturbseq-gene-level jobid: ${gene_level_jobid}"

cnmf_jobid="$(sbatch --parsable plot_perturbseq_cnmf.sh)"
echo "plot-perturbseq-cnmf jobid: ${cnmf_jobid}"
EOF

    chmod +x "$script_path"
}

write_supplementary_submitter() {
    local script_path="$1"

    cat > "$script_path" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

jobid="$(sbatch --parsable plot_supplementary.sh)"
echo "plot-supplementary jobid: ${jobid}"
EOF

    chmod +x "$script_path"
}

write_root_submitter() {
    local script_path="$1"

    cat > "$script_path" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

shopt -s nullglob
for trait_dir in traits/*; do
    if [[ -x "${trait_dir}/submit_trait_plots.sh" ]]; then
        echo "[${trait_dir}]"
        "${trait_dir}/submit_trait_plots.sh"
    fi
done
shopt -u nullglob

if [[ "${INCLUDE_SUPPLEMENTARY:-0}" == "1" ]]; then
    supplementary/submit_supplementary.sh
else
    echo "Skipping plot-supplementary; set INCLUDE_SUPPLEMENTARY=1 to submit it."
fi
EOF

    chmod +x "$script_path"
}

if [[ -n "$PLOT_TRANS_EQTL_MATCH_TABLE" && "$PLOT_TRANS_EQTL_MATCH_TABLE" != /* ]]; then
    PLOT_TRANS_EQTL_MATCH_TABLE="$(resolve_path "$PLOT_TRANS_EQTL_MATCH_TABLE")"
fi
declare -A TRANS_EQTL_GWAS_ID_BY_LOF=()
declare -A TRANS_EQTL_GWAS_PATH_BY_LOF=()
load_trans_eqtl_match_table "$PLOT_TRANS_EQTL_MATCH_TABLE"

declare -A SEEN_PLOT_SAFE_NAMES=()
cnmf_mode="$(yaml_mode "${PLOT_CNMF_MODE}")"
for trait_num in "${trait_numbers[@]}"; do
    trait_file="$(trait_file_for "${trait_num}")"
    raw_burden_file="$(raw_burden_file_for "${trait_num}")"
    trait_label="$(trait_label_for "${trait_num}")"
    plot_safe="$(plot_safe_name "${trait_num}" "$trait_file")"
    if [[ -n "${SEEN_PLOT_SAFE_NAMES[$plot_safe]:-}" ]]; then
        echo "Duplicate normalized trait name '${plot_safe}'. Adjust trait templates to avoid folder collisions." >&2
        exit 1
    fi
    SEEN_PLOT_SAFE_NAMES["$plot_safe"]=1
    trait_dir="${PLOT_TRAITS_ROOT}/${plot_safe}"
    mkdir -p "${trait_dir}"
    lof_id="$(lof_trait_id "$trait_file")"
    matched_gwas="0"
    trans_eqtl_gwas_trait=""
    trans_eqtl_gwas_file=""

    if [[ -n "${TRANS_EQTL_GWAS_PATH_BY_LOF[$lof_id]:-}" ]]; then
        matched_gwas="1"
        trans_eqtl_gwas_trait="${TRANS_EQTL_GWAS_ID_BY_LOF[$lof_id]}"
        trans_eqtl_gwas_file="${TRANS_EQTL_GWAS_PATH_BY_LOF[$lof_id]}"
    fi
    run_multiple_regression="$(multiple_regression_flag "$cnmf_mode")"
    run_trans_eqtl_follow_up="$(trans_eqtl_run_flag "$cnmf_mode" "$matched_gwas")"

    if [[ "$run_trans_eqtl_follow_up" == "true" ]]; then
        if [[ "$matched_gwas" != "1" ]]; then
            echo "trans-eQTL requested for ${trait_file}, but no GWAS match found for ${lof_id} in ${PLOT_TRANS_EQTL_MATCH_TABLE}" >&2
            exit 1
        fi
        if [[ "$trans_eqtl_gwas_file" = /* ]]; then
            [[ -f "$trans_eqtl_gwas_file" ]] || {
                echo "Missing trans-eQTL GWAS file for ${trait_file}: ${trans_eqtl_gwas_file}" >&2
                exit 1
            }
        else
            [[ -f "$(resolve_path "${PLOT_GWAS_DIR}/${trans_eqtl_gwas_file}")" ]] || {
                echo "Missing trans-eQTL GWAS file for ${trait_file}: ${PLOT_GWAS_DIR}/${trans_eqtl_gwas_file}" >&2
                exit 1
            }
        fi
    elif [[ "${PLOT_ENABLE_TRANS_EQTL,,}" =~ ^(inherit|auto)?$ && "$cnmf_mode" == "legacy" && "$matched_gwas" == "0" ]]; then
        echo "WARN: trans-eQTL disabled for ${trait_file}: no GWAS match for ${lof_id}" >&2
    fi

    write_trait_plot_config "${trait_dir}/${PLOT_CONFIG_NAME}" "${trait_num}" "${trait_file}" "${raw_burden_file}" "${trait_label}" "${cnmf_mode}" "${run_multiple_regression}" "${run_trans_eqtl_follow_up}" "${trans_eqtl_gwas_trait}" "${trans_eqtl_gwas_file}"

    write_slurm_script \
        "${trait_dir}/plot_gene_burden.sh" \
        "gene_burden_${plot_safe}" \
        "${PLOT_GENE_BURDEN_MEM}" \
        "${PLOT_GENE_BURDEN_CPUS}" \
        "${PLOT_PARTITION}" \
        "${PLOT_TIME}" \
        "plot-gene-burden"

    write_slurm_script \
        "${trait_dir}/plot_perturbseq_gene_level.sh" \
        "gene_level_${plot_safe}" \
        "${PLOT_GENE_LEVEL_MEM}" \
        "${PLOT_GENE_LEVEL_CPUS}" \
        "${PLOT_PARTITION}" \
        "${PLOT_TIME}" \
        "plot-perturbseq-gene-level"

    write_slurm_script \
        "${trait_dir}/plot_perturbseq_cnmf.sh" \
        "plot_cnmf_${plot_safe}" \
        "${PLOT_CNMF_MEM}" \
        "${PLOT_CNMF_CPUS}" \
        "${PLOT_PARTITION}" \
        "${PLOT_TIME}" \
        "plot-perturbseq-cnmf"

    write_trait_submitter "${trait_dir}/submit_trait_plots.sh"
done

write_supplementary_plot_config "${PLOT_SUPPLEMENTARY_ROOT}/${PLOT_CONFIG_NAME}"

write_slurm_script \
    "${PLOT_SUPPLEMENTARY_ROOT}/plot_supplementary.sh" \
    "plot_supplementary" \
    "${PLOT_SUPPLEMENTARY_MEM}" \
    "${PLOT_SUPPLEMENTARY_CPUS}" \
    "${PLOT_PARTITION}" \
    "${PLOT_SUPPLEMENTARY_TIME}" \
    "plot-supplementary"
write_supplementary_submitter "${PLOT_SUPPLEMENTARY_ROOT}/submit_supplementary.sh"

write_root_submitter "${PLOT_ROOT}/submit_plot_workflows.sh"

echo "Generated plot scripts under:"
echo "  traits: ${PLOT_TRAITS_ROOT}/*"
echo "  supplementary: ${PLOT_SUPPLEMENTARY_ROOT}"
echo "  tools:  ${PLOT_TOOLS_ROOT}"
echo "  logs:   ${LOGS_ROOT}"
echo "  cnmf mode: ${cnmf_mode}"
echo "  trans-eQTL: ${PLOT_ENABLE_TRANS_EQTL} (match table: ${PLOT_TRANS_EQTL_MATCH_TABLE:-none})"
