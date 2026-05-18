#!/usr/bin/env bash

set -euo pipefail

START_TRAIT="${START_TRAIT:-1}"
END_TRAIT="${END_TRAIT:-150}"
DEFAULT_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_ROOT="${BATCH_ROOT:-${DEFAULT_PROJECT_ROOT}/scripts/run_once}"
SHARED_ROOT="${SHARED_ROOT:-${BATCH_ROOT}/shared}"
LOGS_ROOT="${LOGS_ROOT:-${BATCH_ROOT}/logs}"
SHARED_CONFIG_NAME="${SHARED_CONFIG_NAME:-pipeline.shared.yaml}"

TRAIT_PREFIX="${TRAIT_PREFIX:-Backman_2021_}"
SHARED_TRAIT_FILE="${SHARED_TRAIT_FILE:-${TRAIT_PREFIX}${START_TRAIT}.per_gene_estimates.tsv}"

CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
PERTURBSEQ_ENV="${PERTURBSEQ_ENV:-paper-pipeline-perturbseq}"
R_ENV="${R_ENV:-paper-pipeline-r}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

GENE_LEVEL_PREP_MEM="${GENE_LEVEL_PREP_MEM:-700G}"
GENE_LEVEL_PREP_CPUS="${GENE_LEVEL_PREP_CPUS:-20}"
GENE_LEVEL_PREP_PARTITION="${GENE_LEVEL_PREP_PARTITION:-fat}"
GENE_LEVEL_PREP_TIME="${GENE_LEVEL_PREP_TIME:-7-00:00:00}"

GENE_LEVEL_LIMMA_MEM="${GENE_LEVEL_LIMMA_MEM:-50G}"
GENE_LEVEL_LIMMA_CPUS="${GENE_LEVEL_LIMMA_CPUS:-10}"
GENE_LEVEL_LIMMA_PARTITION="${GENE_LEVEL_LIMMA_PARTITION:-cu,privority,batch01}"
GENE_LEVEL_LIMMA_TIME="${GENE_LEVEL_LIMMA_TIME:-7-00:00:00}"
GENE_LEVEL_TOTAL_CHUNKS="${GENE_LEVEL_TOTAL_CHUNKS:-198}"
GENE_LEVEL_CHUNKS_PER_JOB="${GENE_LEVEL_CHUNKS_PER_JOB:-5}"

GENE_LEVEL_SUMMARY_MEM="${GENE_LEVEL_SUMMARY_MEM:-40G}"
GENE_LEVEL_SUMMARY_CPUS="${GENE_LEVEL_SUMMARY_CPUS:-4}"
GENE_LEVEL_SUMMARY_PARTITION="${GENE_LEVEL_SUMMARY_PARTITION:-cu,privority,batch01}"
GENE_LEVEL_SUMMARY_TIME="${GENE_LEVEL_SUMMARY_TIME:-7-00:00:00}"

ESSENTIAL_PREP_MEM="${ESSENTIAL_PREP_MEM:-450G}"
ESSENTIAL_PREP_CPUS="${ESSENTIAL_PREP_CPUS:-20}"
ESSENTIAL_PREP_PARTITION="${ESSENTIAL_PREP_PARTITION:-fat}"
ESSENTIAL_PREP_TIME="${ESSENTIAL_PREP_TIME:-7-00:00:00}"
ESSENTIAL_K="${ESSENTIAL_K:-60}"
ESSENTIAL_TOTAL_WORKERS="${ESSENTIAL_TOTAL_WORKERS:-20}"
ESSENTIAL_WORKERS_PER_JOB="${ESSENTIAL_WORKERS_PER_JOB:-2}"
ESSENTIAL_PROGRAMS_PER_JOB="${ESSENTIAL_PROGRAMS_PER_JOB:-10}"

ESSENTIAL_FACTORIZE_MEM="${ESSENTIAL_FACTORIZE_MEM:-40G}"
ESSENTIAL_FACTORIZE_CPUS="${ESSENTIAL_FACTORIZE_CPUS:-10}"
ESSENTIAL_FACTORIZE_PARTITION="${ESSENTIAL_FACTORIZE_PARTITION:-cu,privority,batch01}"
ESSENTIAL_FACTORIZE_TIME="${ESSENTIAL_FACTORIZE_TIME:-10-00:00:00}"

ESSENTIAL_POST_MEM="${ESSENTIAL_POST_MEM:-100G}"
ESSENTIAL_POST_CPUS="${ESSENTIAL_POST_CPUS:-40}"
ESSENTIAL_POST_PARTITION="${ESSENTIAL_POST_PARTITION:-cu,privority,batch01}"
ESSENTIAL_POST_TIME="${ESSENTIAL_POST_TIME:-7-00:00:00}"

GENOMEWIDE_PREP_MEM="${GENOMEWIDE_PREP_MEM:-450G}"
GENOMEWIDE_PREP_CPUS="${GENOMEWIDE_PREP_CPUS:-10}"
GENOMEWIDE_PREP_PARTITION="${GENOMEWIDE_PREP_PARTITION:-fat}"
GENOMEWIDE_PREP_TIME="${GENOMEWIDE_PREP_TIME:-7-00:00:00}"
GENOMEWIDE_ASSOCIATION_K="${GENOMEWIDE_ASSOCIATION_K:-60}"
GENOMEWIDE_TOTAL_WORKERS="${GENOMEWIDE_TOTAL_WORKERS:-100}"
GENOMEWIDE_WORKERS_PER_JOB="${GENOMEWIDE_WORKERS_PER_JOB:-2}"
GENOMEWIDE_PROGRAMS_PER_JOB="${GENOMEWIDE_PROGRAMS_PER_JOB:-1}"

GENOMEWIDE_FACTORIZE_MEM="${GENOMEWIDE_FACTORIZE_MEM:-80G}"
GENOMEWIDE_FACTORIZE_CPUS="${GENOMEWIDE_FACTORIZE_CPUS:-16}"
GENOMEWIDE_FACTORIZE_PARTITION="${GENOMEWIDE_FACTORIZE_PARTITION:-cu,privority,batch01}"
GENOMEWIDE_FACTORIZE_TIME="${GENOMEWIDE_FACTORIZE_TIME:-10-00:00:00}"

GENOMEWIDE_POST_MEM="${GENOMEWIDE_POST_MEM:-60G}"
GENOMEWIDE_POST_CPUS="${GENOMEWIDE_POST_CPUS:-15}"
GENOMEWIDE_POST_PARTITION="${GENOMEWIDE_POST_PARTITION:-cu,privority,batch01}"
GENOMEWIDE_POST_TIME="${GENOMEWIDE_POST_TIME:-7-00:00:00}"

if [[ "$START_TRAIT" =~ [^0-9] ]] || [[ "$END_TRAIT" =~ [^0-9] ]]; then
    echo "START_TRAIT and END_TRAIT must be integers" >&2
    exit 1
fi

if (( START_TRAIT > END_TRAIT )); then
    echo "START_TRAIT must be <= END_TRAIT" >&2
    exit 1
fi

for value_name in \
    GENE_LEVEL_TOTAL_CHUNKS \
    GENE_LEVEL_CHUNKS_PER_JOB \
    ESSENTIAL_K \
    ESSENTIAL_TOTAL_WORKERS \
    ESSENTIAL_WORKERS_PER_JOB \
    ESSENTIAL_PROGRAMS_PER_JOB \
    GENOMEWIDE_ASSOCIATION_K \
    GENOMEWIDE_TOTAL_WORKERS \
    GENOMEWIDE_WORKERS_PER_JOB \
    GENOMEWIDE_PROGRAMS_PER_JOB; do
    value="${!value_name}"
    if [[ "$value" =~ [^0-9] ]] || (( value < 1 )); then
        echo "${value_name} must be a positive integer" >&2
        exit 1
    fi
done

mkdir -p "${BATCH_ROOT}"
BATCH_ROOT="$(cd "${BATCH_ROOT}" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-}"
if [[ -z "${PROJECT_ROOT}" ]]; then
    PROJECT_ROOT="$(cd "${BATCH_ROOT}/../.." && pwd)"
fi
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${PROJECT_ROOT}/outputs}"
if [[ "${ARTIFACT_ROOT}" != /* ]]; then
    ARTIFACT_ROOT="$(cd "${PROJECT_ROOT}" && realpath -m "${ARTIFACT_ROOT}")"
fi
mkdir -p "${SHARED_ROOT}" "${LOGS_ROOT}"
SHARED_ROOT="$(cd "${SHARED_ROOT}" && pwd)"
LOGS_ROOT="$(cd "${LOGS_ROOT}" && pwd)"
if [[ -z "${SHARED_ROOT}" || "${SHARED_ROOT}" == "/" ]]; then
    echo "FATAL: refusing to rm -rf empty or root path" >&2; exit 1
fi
rm -rf "${SHARED_ROOT}"
mkdir -p "${SHARED_ROOT}"

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

write_slurm_script() {
    local script_path="$1"
    local job_name="$2"
    local mem="$3"
    local cpus="$4"
    local partition="$5"
    local time_limit="$6"
    local workflow_name="$7"
    local gres="${8:-}"
    local config_file="${9:-${SHARED_CONFIG_NAME}}"
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
EOF

    if [[ -n "$gres" ]]; then
        cat >> "$script_path" <<EOF
#SBATCH --gres=${gres}
EOF
    fi

cat >> "$script_path" <<EOF
#SBATCH --time=${time_limit}

set -euo pipefail

source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"

cd "${script_dir}"
paper-pipeline run --config "${config_path}" ${workflow_name}
EOF

    chmod +x "$script_path"
}

indices_yaml() {
    local start="$1"
    local per_job="$2"
    local total="$3"
    local zero_based="${4:-0}"
    local values=()
    local end=$(( start + per_job - 1 ))
    if (( end > total )); then
        end=$total
    fi
    for (( idx=start; idx<=end; idx++ )); do
        if (( zero_based == 1 )); then
            values+=("$(( idx - 1 ))")
        else
            values+=("${idx}")
        fi
    done
    local joined
    printf -v joined '%s, ' "${values[@]}"
    joined="${joined%, }"
    printf '%s' "${joined}"
}

write_shared_config() {
    local config_path="$1"
    local essential_worker_indices="${2:-}"
    local genomewide_worker_indices="${3:-}"
    local gene_level_chunk_indices="${4:-}"
    local essential_program_indices="${5:-}"
    local genomewide_program_indices="${6:-}"

    cat > "$config_path" <<EOF
project_root: ${PROJECT_ROOT}
artifact_root: ${ARTIFACT_ROOT}

executables:
  perturbseq_runner:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${PERTURBSEQ_ENV}
  rscript:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${R_ENV}
    - Rscript

workflows:
  perturbseq:
    gene_level:
      inputs:
        h5ad: $(resolve_path "data/K562_gwps_raw_singlecell_01.h5ad")
        metadata: $(resolve_path "data/Perturbseq/metadata/gwps_metadata.csv")
      parameters:
        filter_min_genes: 500
        filter_min_cells: 500
        chunk_size: 50
        total_chunks: ${GENE_LEVEL_TOTAL_CHUNKS}
EOF

    if [[ -n "${gene_level_chunk_indices}" ]]; then
        cat >> "$config_path" <<EOF
        chunk_indices: [${gene_level_chunk_indices}]
EOF
    fi

    cat >> "$config_path" <<EOF

    cnmf_essential:
      cell_types: [K562_essential_raw_singlecell_01]
      trait_files: [${SHARED_TRAIT_FILE}]
      inputs:
        h5ad_dir: $(resolve_path "data")
        metadata_dir: $(resolve_path "data/Perturbseq/metadata")
        gene_map: $(resolve_path "data/gencode_v41_gname_gid_ALL_sorted_onlyID")
        shet_path: $(resolve_path "data/shet_10bins.txt")
      parameters:
        k: ${ESSENTIAL_K}
        total_workers: ${ESSENTIAL_TOTAL_WORKERS}
EOF

    if [[ -n "${essential_worker_indices}" ]]; then
        cat >> "$config_path" <<EOF
        worker_indices: [${essential_worker_indices}]
EOF
    fi

    if [[ -n "${essential_program_indices}" ]]; then
        cat >> "$config_path" <<EOF
        program_indices: [${essential_program_indices}]
EOF
    fi

    cat >> "$config_path" <<EOF
        n_iter: 100
        seed: 14
        numgenes: 2000
        consensus_density_threshold: 0.4
        mito_column: mitopercent
        mito_threshold: 0.3
        min_genes_per_cell: 100
        min_counts_per_cell: 100
        min_cells_per_gene: 10
        random_iterations: 100000
        top_n_program_genes: 100

    cnmf_genomewide:
      trait_files: [${SHARED_TRAIT_FILE}]
      inputs:
        h5ad: $(resolve_path "data/K562_gwps_raw_singlecell_01.h5ad")
        metadata: $(resolve_path "data/Perturbseq/metadata/gwps_metadata.csv")
        gene_map: $(resolve_path "data/gencode_v41_gname_gid_ALL_sorted_onlyID")
        shet_path: $(resolve_path "data/shet_10bins.txt")
      parameters:
        ks: [30, 60, 90, 120]
        total_workers: ${GENOMEWIDE_TOTAL_WORKERS}
EOF

    if [[ -n "${genomewide_worker_indices}" ]]; then
        cat >> "$config_path" <<EOF
        worker_indices: [${genomewide_worker_indices}]
EOF
    fi

    cat >> "$config_path" <<EOF
        n_iter: 100
        seed: 14
        numgenes: 2000
        aggregate_name: cNMF_all
        per_k_name_template: cNMF_K{k}
        consensus_ks: [${GENOMEWIDE_ASSOCIATION_K}]
        consensus_density_threshold: 0.5
        association_k: ${GENOMEWIDE_ASSOCIATION_K}
        random_iterations: 100000
        top_n_program_genes: 100
EOF

    if [[ -n "${genomewide_program_indices}" ]]; then
        cat >> "$config_path" <<EOF
        program_indices: [${genomewide_program_indices}]
EOF
    fi
}

write_dependency_submitter() {
    local script_path="$1"
    local prepare_script="$2"
    local factorize_glob="$3"
    local postbase_script="$4"

    cat > "$script_path" <<EOF
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
cd "\$SCRIPT_DIR/shared"

prepare_jobid="\$(sbatch --parsable "${prepare_script}")"
echo "prepare jobid: \${prepare_jobid}"

factorize_jobids=()
shopt -s nullglob
for script in ${factorize_glob}; do
    jobid="\$(sbatch --parsable --dependency=afterok:\${prepare_jobid} "\${script}")"
    factorize_jobids+=("\${jobid}")
    echo "factorize jobid: \${jobid} (\${script})"
done
shopt -u nullglob

if (( \${#factorize_jobids[@]} == 0 )); then
    echo "No factorize scripts matched: ${factorize_glob}" >&2
    exit 1
fi

dependency="\$(IFS=:; echo "\${factorize_jobids[*]}")"
postbase_jobid="\$(sbatch --parsable --dependency=afterok:\${dependency} "${postbase_script}")"
echo "postbase jobid: \${postbase_jobid}"
EOF

    chmod +x "$script_path"
}

write_split_postbase_dependency_submitter() {
    local script_path="$1"
    local prepare_script="$2"
    local factorize_glob="$3"
    local postbase_base_script="$4"
    local regulation_glob="$5"

    cat > "$script_path" <<EOF
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
cd "\$SCRIPT_DIR/shared"

prepare_jobid="\$(sbatch --parsable "${prepare_script}")"
echo "prepare jobid: \${prepare_jobid}"

factorize_jobids=()
shopt -s nullglob
for script in ${factorize_glob}; do
    jobid="\$(sbatch --parsable --dependency=afterok:\${prepare_jobid} "\${script}")"
    factorize_jobids+=("\${jobid}")
    echo "factorize jobid: \${jobid} (\${script})"
done
shopt -u nullglob

if (( \${#factorize_jobids[@]} == 0 )); then
    echo "No factorize scripts matched: ${factorize_glob}" >&2
    exit 1
fi

factorize_dependency="\$(IFS=:; echo "\${factorize_jobids[*]}")"
postbase_base_jobid="\$(sbatch --parsable --dependency=afterok:\${factorize_dependency} "${postbase_base_script}")"
echo "postbase base jobid: \${postbase_base_jobid}"

regulation_jobids=()
shopt -s nullglob
for script in ${regulation_glob}; do
    jobid="\$(sbatch --parsable --dependency=afterok:\${postbase_base_jobid} "\${script}")"
    regulation_jobids+=("\${jobid}")
    echo "postbase regulation jobid: \${jobid} (\${script})"
done
shopt -u nullglob

if (( \${#regulation_jobids[@]} == 0 )); then
    echo "No postbase regulation scripts matched: ${regulation_glob}" >&2
    exit 1
fi
EOF

    chmod +x "$script_path"
}

write_shared_config "${SHARED_ROOT}/${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_gene_level_prepare.sh" \
    "pertseq_genelevel_prepare" \
    "${GENE_LEVEL_PREP_MEM}" \
    "${GENE_LEVEL_PREP_CPUS}" \
    "${GENE_LEVEL_PREP_PARTITION}" \
    "${GENE_LEVEL_PREP_TIME}" \
    "perturbseq-gene-level-prepare" \
    "" \
    "${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_gene_level_summarize.sh" \
    "pertseq_genelevel_summarize" \
    "${GENE_LEVEL_SUMMARY_MEM}" \
    "${GENE_LEVEL_SUMMARY_CPUS}" \
    "${GENE_LEVEL_SUMMARY_PARTITION}" \
    "${GENE_LEVEL_SUMMARY_TIME}" \
    "perturbseq-gene-level-summarize" \
    "" \
    "${SHARED_CONFIG_NAME}"

for (( chunk_start=1; chunk_start<=GENE_LEVEL_TOTAL_CHUNKS; chunk_start+=GENE_LEVEL_CHUNKS_PER_JOB )); do
    chunk_indices="$(indices_yaml "${chunk_start}" "${GENE_LEVEL_CHUNKS_PER_JOB}" "${GENE_LEVEL_TOTAL_CHUNKS}" 0)"
    chunk_last=$(( chunk_start + GENE_LEVEL_CHUNKS_PER_JOB - 1 ))
    if (( chunk_last > GENE_LEVEL_TOTAL_CHUNKS )); then
        chunk_last=${GENE_LEVEL_TOTAL_CHUNKS}
    fi
    gene_level_cfg="pipeline.gene_level.limma_${chunk_start}_${chunk_last}.yaml"
    gene_level_script="perturbseq_gene_level_limma_${chunk_start}_${chunk_last}.sh"
    write_shared_config "${SHARED_ROOT}/${gene_level_cfg}" "" "" "${chunk_indices}"
    write_slurm_script \
        "${SHARED_ROOT}/${gene_level_script}" \
        "pertseq_genelevel_limma_${chunk_start}_${chunk_last}" \
        "${GENE_LEVEL_LIMMA_MEM}" \
        "${GENE_LEVEL_LIMMA_CPUS}" \
        "${GENE_LEVEL_LIMMA_PARTITION}" \
        "${GENE_LEVEL_LIMMA_TIME}" \
        "perturbseq-gene-level-limma" \
        "" \
        "${gene_level_cfg}"
done

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_essential_prepare.sh" \
    "pertseq_cnmf_essential_prepare" \
    "${ESSENTIAL_PREP_MEM}" \
    "${ESSENTIAL_PREP_CPUS}" \
    "${ESSENTIAL_PREP_PARTITION}" \
    "${ESSENTIAL_PREP_TIME}" \
    "perturbseq-cnmf-essential-prepare" \
    "" \
    "${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_essential_postbase.sh" \
    "pertseq_cnmf_essential_postbase" \
    "${ESSENTIAL_POST_MEM}" \
    "${ESSENTIAL_POST_CPUS}" \
    "${ESSENTIAL_POST_PARTITION}" \
    "${ESSENTIAL_POST_TIME}" \
    "perturbseq-cnmf-essential-postbase" \
    "" \
    "${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_essential_postbase_base.sh" \
    "pertseq_cnmf_essential_postbase_base" \
    "${ESSENTIAL_POST_MEM}" \
    "${ESSENTIAL_POST_CPUS}" \
    "${ESSENTIAL_POST_PARTITION}" \
    "${ESSENTIAL_POST_TIME}" \
    "perturbseq-cnmf-essential-postbase-base" \
    "" \
    "${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_genomewide_prepare.sh" \
    "pertseq_cnmf_genomewide_prepare" \
    "${GENOMEWIDE_PREP_MEM}" \
    "${GENOMEWIDE_PREP_CPUS}" \
    "${GENOMEWIDE_PREP_PARTITION}" \
    "${GENOMEWIDE_PREP_TIME}" \
    "perturbseq-cnmf-genomewide-prepare" \
    "" \
    "${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_genomewide_postbase.sh" \
    "pertseq_cnmf_genomewide_postbase" \
    "${GENOMEWIDE_POST_MEM}" \
    "${GENOMEWIDE_POST_CPUS}" \
    "${GENOMEWIDE_POST_PARTITION}" \
    "${GENOMEWIDE_POST_TIME}" \
    "perturbseq-cnmf-genomewide-postbase" \
    "" \
    "${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_genomewide_postbase_base.sh" \
    "pertseq_cnmf_genomewide_postbase_base" \
    "${GENOMEWIDE_POST_MEM}" \
    "${GENOMEWIDE_POST_CPUS}" \
    "${GENOMEWIDE_POST_PARTITION}" \
    "${GENOMEWIDE_POST_TIME}" \
    "perturbseq-cnmf-genomewide-postbase-base" \
    "" \
    "${SHARED_CONFIG_NAME}"

for (( worker_start=1; worker_start<=ESSENTIAL_TOTAL_WORKERS; worker_start+=ESSENTIAL_WORKERS_PER_JOB )); do
    essential_indices="$(indices_yaml "${worker_start}" "${ESSENTIAL_WORKERS_PER_JOB}" "${ESSENTIAL_TOTAL_WORKERS}" 1)"
    essential_last=$(( worker_start + ESSENTIAL_WORKERS_PER_JOB - 1 ))
    if (( essential_last > ESSENTIAL_TOTAL_WORKERS )); then
        essential_last=${ESSENTIAL_TOTAL_WORKERS}
    fi
    essential_cfg="pipeline.essential.factorize_${worker_start}_${essential_last}.yaml"
    essential_script="perturbseq_cnmf_essential_factorize_${worker_start}_${essential_last}.sh"
    write_shared_config "${SHARED_ROOT}/${essential_cfg}" "${essential_indices}" "" ""
    write_slurm_script \
        "${SHARED_ROOT}/${essential_script}" \
        "pertseq_cnmf_essential_factorize_${worker_start}_${essential_last}" \
        "${ESSENTIAL_FACTORIZE_MEM}" \
        "${ESSENTIAL_FACTORIZE_CPUS}" \
        "${ESSENTIAL_FACTORIZE_PARTITION}" \
        "${ESSENTIAL_FACTORIZE_TIME}" \
        "perturbseq-cnmf-essential-factorize" \
        "" \
        "${essential_cfg}"
done

for (( program_start=1; program_start<=ESSENTIAL_K; program_start+=ESSENTIAL_PROGRAMS_PER_JOB )); do
    essential_program_indices="$(indices_yaml "${program_start}" "${ESSENTIAL_PROGRAMS_PER_JOB}" "${ESSENTIAL_K}" 0)"
    program_last=$(( program_start + ESSENTIAL_PROGRAMS_PER_JOB - 1 ))
    if (( program_last > ESSENTIAL_K )); then
        program_last=${ESSENTIAL_K}
    fi
    essential_cfg="pipeline.essential.postbase_regulation_${program_start}_${program_last}.yaml"
    essential_script="perturbseq_cnmf_essential_postbase_regulation_${program_start}_${program_last}.sh"
    write_shared_config "${SHARED_ROOT}/${essential_cfg}" "" "" "" "${essential_program_indices}" ""
    write_slurm_script \
        "${SHARED_ROOT}/${essential_script}" \
        "pertseq_cnmf_essential_postreg_${program_start}_${program_last}" \
        "${ESSENTIAL_POST_MEM}" \
        "${ESSENTIAL_POST_CPUS}" \
        "${ESSENTIAL_POST_PARTITION}" \
        "${ESSENTIAL_POST_TIME}" \
        "perturbseq-cnmf-essential-postbase-regulation" \
        "" \
        "${essential_cfg}"
done

for (( worker_start=1; worker_start<=GENOMEWIDE_TOTAL_WORKERS; worker_start+=GENOMEWIDE_WORKERS_PER_JOB )); do
    genomewide_indices="$(indices_yaml "${worker_start}" "${GENOMEWIDE_WORKERS_PER_JOB}" "${GENOMEWIDE_TOTAL_WORKERS}" 1)"
    genomewide_last=$(( worker_start + GENOMEWIDE_WORKERS_PER_JOB - 1 ))
    if (( genomewide_last > GENOMEWIDE_TOTAL_WORKERS )); then
        genomewide_last=${GENOMEWIDE_TOTAL_WORKERS}
    fi
    genomewide_cfg="pipeline.genomewide.factorize_${worker_start}_${genomewide_last}.yaml"
    genomewide_script="perturbseq_cnmf_genomewide_factorize_${worker_start}_${genomewide_last}.sh"
    write_shared_config "${SHARED_ROOT}/${genomewide_cfg}" "" "${genomewide_indices}" ""
    write_slurm_script \
        "${SHARED_ROOT}/${genomewide_script}" \
        "pertseq_cnmf_genomewide_factorize_${worker_start}_${genomewide_last}" \
        "${GENOMEWIDE_FACTORIZE_MEM}" \
        "${GENOMEWIDE_FACTORIZE_CPUS}" \
        "${GENOMEWIDE_FACTORIZE_PARTITION}" \
        "${GENOMEWIDE_FACTORIZE_TIME}" \
        "perturbseq-cnmf-genomewide-factorize" \
        "" \
        "${genomewide_cfg}"
done

for (( program_start=1; program_start<=GENOMEWIDE_ASSOCIATION_K; program_start+=GENOMEWIDE_PROGRAMS_PER_JOB )); do
    genomewide_program_indices="$(indices_yaml "${program_start}" "${GENOMEWIDE_PROGRAMS_PER_JOB}" "${GENOMEWIDE_ASSOCIATION_K}" 0)"
    program_last=$(( program_start + GENOMEWIDE_PROGRAMS_PER_JOB - 1 ))
    if (( program_last > GENOMEWIDE_ASSOCIATION_K )); then
        program_last=${GENOMEWIDE_ASSOCIATION_K}
    fi
    genomewide_cfg="pipeline.genomewide.postbase_regulation_${program_start}_${program_last}.yaml"
    genomewide_script="perturbseq_cnmf_genomewide_postbase_regulation_${program_start}_${program_last}.sh"
    write_shared_config "${SHARED_ROOT}/${genomewide_cfg}" "" "" "" "" "${genomewide_program_indices}"
    write_slurm_script \
        "${SHARED_ROOT}/${genomewide_script}" \
        "pertseq_cnmf_gw_postreg_${program_start}_${program_last}" \
        "${GENOMEWIDE_POST_MEM}" \
        "${GENOMEWIDE_POST_CPUS}" \
        "${GENOMEWIDE_POST_PARTITION}" \
        "${GENOMEWIDE_POST_TIME}" \
        "perturbseq-cnmf-genomewide-postbase-regulation" \
        "" \
        "${genomewide_cfg}"
done

write_split_postbase_dependency_submitter \
    "${BATCH_ROOT}/submit_essential_chain.sh" \
    "perturbseq_cnmf_essential_prepare.sh" \
    "perturbseq_cnmf_essential_factorize_*.sh" \
    "perturbseq_cnmf_essential_postbase_base.sh" \
    "perturbseq_cnmf_essential_postbase_regulation_*.sh"

write_split_postbase_dependency_submitter \
    "${BATCH_ROOT}/submit_genomewide_chain.sh" \
    "perturbseq_cnmf_genomewide_prepare.sh" \
    "perturbseq_cnmf_genomewide_factorize_*.sh" \
    "perturbseq_cnmf_genomewide_postbase_base.sh" \
    "perturbseq_cnmf_genomewide_postbase_regulation_*.sh"

cat > "${BATCH_ROOT}/submit_gene_level_chain.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/shared"

prepare_jobid="$(sbatch --parsable perturbseq_gene_level_prepare.sh)"
echo "prepare jobid: ${prepare_jobid}"

limma_jobids=()
shopt -s nullglob
for script in perturbseq_gene_level_limma_*.sh; do
    jobid="$(sbatch --parsable --dependency=afterok:${prepare_jobid} "${script}")"
    limma_jobids+=("${jobid}")
    echo "limma jobid: ${jobid} (${script})"
done
shopt -u nullglob

if (( ${#limma_jobids[@]} == 0 )); then
    echo "No gene-level limma scripts found" >&2
    exit 1
fi

dependency="$(IFS=:; echo "${limma_jobids[*]}")"
summarize_jobid="$(sbatch --parsable --dependency=afterok:${dependency} perturbseq_gene_level_summarize.sh)"
echo "summarize jobid: ${summarize_jobid}"
EOF

chmod +x "${BATCH_ROOT}/submit_gene_level_chain.sh"

echo "Generated shared run-once scripts under: ${BATCH_ROOT}"
