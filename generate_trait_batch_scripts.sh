#!/usr/bin/env bash

set -euo pipefail

START_TRAIT="${START_TRAIT:-1}"
END_TRAIT="${END_TRAIT:-150}"
DEFAULT_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_ROOT="${BATCH_ROOT:-${DEFAULT_PROJECT_ROOT}/scripts/traits}"
CONFIG_NAME="${CONFIG_NAME:-pipeline.yaml}"
LOGS_ROOT="${LOGS_ROOT:-${BATCH_ROOT}/logs}"

TRAIT_PREFIX="${TRAIT_PREFIX:-Backman_2021_}"

CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
GENEBAYES_ENV="${GENEBAYES_ENV:-}"
GENEBAYES_GPU_ENV_IS_SET="${GENEBAYES_GPU_ENV+x}"
GENEBAYES_CPU_ENV="${GENEBAYES_CPU_ENV:-${GENEBAYES_ENV:-paper-pipeline-genebayes}}"
GENEBAYES_GPU_ENV="${GENEBAYES_GPU_ENV:-paper-pipeline-genebayes-gpu}"
PERTURBSEQ_ENV="${PERTURBSEQ_ENV:-paper-pipeline-perturbseq}"
R_ENV="${R_ENV:-paper-pipeline-r}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

GENEBAYES_MEM="${GENEBAYES_MEM:-200G}"
GENEBAYES_CPUS="${GENEBAYES_CPUS:-50}"
GENEBAYES_PARTITION="${GENEBAYES_PARTITION:-A800,gpu}"
GENEBAYES_GRES="${GENEBAYES_GRES:-gpu:1}"
GENEBAYES_TIME="${GENEBAYES_TIME:-10-00:00:00}"
GENEBAYES_CPU_MEM="${GENEBAYES_CPU_MEM:-60G}"
GENEBAYES_CPU_CPUS="${GENEBAYES_CPU_CPUS:-16}"
GENEBAYES_CPU_PARTITION="${GENEBAYES_CPU_PARTITION:-cu,privority,batch01,fat}"
GENEBAYES_CPU_TIME="${GENEBAYES_CPU_TIME:-7-00:00:00}"
GENEBAYES_GPU_MEM="${GENEBAYES_GPU_MEM:-${GENEBAYES_MEM}}"
GENEBAYES_GPU_CPUS="${GENEBAYES_GPU_CPUS:-${GENEBAYES_CPUS}}"
GENEBAYES_GPU_PARTITION="${GENEBAYES_GPU_PARTITION:-${GENEBAYES_PARTITION}}"
GENEBAYES_GPU_GRES="${GENEBAYES_GPU_GRES:-${GENEBAYES_GRES}}"
GENEBAYES_GPU_TIME="${GENEBAYES_GPU_TIME:-${GENEBAYES_TIME}}"
GENEBAYES_TRAITS_PER_GPU_JOB="${GENEBAYES_TRAITS_PER_GPU_JOB:-12}"
GENEBAYES_THREADS_PER_TRAIT="${GENEBAYES_THREADS_PER_TRAIT:-}"
GENEBAYES_PACKED_ROOT="${GENEBAYES_PACKED_ROOT:-${BATCH_ROOT}/genebayes_packed}"
GENEBAYES_PARALLEL_A800="${GENEBAYES_PARALLEL_A800:-3}"
GENEBAYES_PARALLEL_GPU="${GENEBAYES_PARALLEL_GPU:-2}"

ESSENTIAL_ASSOC_MEM="${ESSENTIAL_ASSOC_MEM:-10G}"
ESSENTIAL_ASSOC_CPUS="${ESSENTIAL_ASSOC_CPUS:-2}"
ESSENTIAL_ASSOC_PARTITION="${ESSENTIAL_ASSOC_PARTITION:-cu,privority,batch01}"
ESSENTIAL_ASSOC_TIME="${ESSENTIAL_ASSOC_TIME:-7-00:00:00}"

GENOMEWIDE_ASSOC_MEM="${GENOMEWIDE_ASSOC_MEM:-30G}"
GENOMEWIDE_ASSOC_CPUS="${GENOMEWIDE_ASSOC_CPUS:-2}"
GENOMEWIDE_ASSOC_PARTITION="${GENOMEWIDE_ASSOC_PARTITION:-cu,privority,batch01,fat}"
GENOMEWIDE_ASSOC_TIME="${GENOMEWIDE_ASSOC_TIME:-7-00:00:00}"

if [[ "$START_TRAIT" =~ [^0-9] ]] || [[ "$END_TRAIT" =~ [^0-9] ]]; then
    echo "START_TRAIT and END_TRAIT must be integers" >&2
    exit 1
fi

if (( START_TRAIT > END_TRAIT )); then
    echo "START_TRAIT must be <= END_TRAIT" >&2
    exit 1
fi

if [[ -n "${GENEBAYES_ENV}" && -z "${GENEBAYES_GPU_ENV_IS_SET}" ]]; then
    echo "GENEBAYES_ENV now only aliases the CPU environment; set GENEBAYES_GPU_ENV explicitly if needed." >&2
fi

for value_name in GENEBAYES_TRAITS_PER_GPU_JOB GENEBAYES_PARALLEL_A800 GENEBAYES_PARALLEL_GPU; do
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
GENEBAYES_PACKED_ROOT="${GENEBAYES_PACKED_ROOT:-${BATCH_ROOT}/genebayes_packed}"
mkdir -p "${GENEBAYES_PACKED_ROOT}"
GENEBAYES_PACKED_ROOT="$(cd "${GENEBAYES_PACKED_ROOT}" && pwd)"
mkdir -p "${LOGS_ROOT}"
LOGS_ROOT="$(cd "${LOGS_ROOT}" && pwd)"

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
    local config_file="${9:-${CONFIG_NAME}}"
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

write_trait_config() {
    local config_path="$1"
    local trait_num="$2"
    local trait_file="$3"

    cat > "$config_path" <<EOF
project_root: ${PROJECT_ROOT}
artifact_root: ${ARTIFACT_ROOT}

executables:
  genebayes_runner:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${GENEBAYES_GPU_ENV}
  genebayes_cpu_runner:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${GENEBAYES_CPU_ENV}
  genebayes_gpu_runner:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${GENEBAYES_GPU_ENV}
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
  genebayes:
    traits: [${trait_num}]
    trait_prefix: ${TRAIT_PREFIX}
    inputs:
      s_het: $(resolve_path "data/LoF/GeneBayes_geneFeatures/s_het.tsv")
      reduced_embeddings: $(resolve_path "data/LoF/GeneBayes_geneFeatures/reduced_embeddings.tsv")
      celltype_ntpm: $(resolve_path "data/LoF/GeneBayes_geneFeatures/celltype_nTPM.tsv")
      geneformer_mixture: $(resolve_path "data/LoF/GeneBayes_geneFeatures/GeneFormer_cellclassifier_mixture2.tsv")
      traits_dir: $(resolve_path "data/LoF/raw_burden")
      train_gene_list: $(resolve_path "data/LoF/raw_burden/10086.csv")
      val_gene_list: $(resolve_path "data/LoF/raw_burden/10086.csv")
    parameters:
      response_file_template: "{trait}.summary_statistics.csv"
      prefer_gpu: true
      integration_lb: -3.99
      integration_ub: 3.99
      batch_size: 100
      n_integration_pts: 50000
      lr: 0.3
      total_iterations: 1000
      early_stopping_iter: 10
      n_trees_per_iteration: 2
      reg_alpha: 2
      reg_lambda: 2
      subsample: 0.8
      min_child_weight: 3
      max_depth: 3

  perturbseq:
    cnmf_essential:
      cell_types: [K562_essential_raw_singlecell_01]
      trait_files: [${trait_file}]
      inputs:
        h5ad_dir: $(resolve_path "data")
        metadata_dir: $(resolve_path "data/Perturbseq/metadata")
        gene_map: $(resolve_path "data/gencode_v41_gname_gid_ALL_sorted_onlyID")
        shet_path: $(resolve_path "data/shet_10bins.txt")
      parameters:
        k: 60
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
      trait_files: [${trait_file}]
      inputs:
        h5ad: $(resolve_path "data/K562_gwps_raw_singlecell_01.h5ad")
        metadata: $(resolve_path "data/Perturbseq/metadata/gwps_metadata.csv")
        gene_map: $(resolve_path "data/gencode_v41_gname_gid_ALL_sorted_onlyID")
        shet_path: $(resolve_path "data/shet_10bins.txt")
      parameters:
        ks: [30, 60, 90, 120]
        n_iter: 100
        seed: 14
        numgenes: 2000
        aggregate_name: cNMF_all
        per_k_name_template: cNMF_K{k}
        consensus_ks: [60]
        consensus_density_threshold: 0.5
        association_k: 60
        random_iterations: 100000
        top_n_program_genes: 100
EOF
}

write_trait_submitter() {
    local script_path="$1"

    cat > "$script_path" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

genebayes_cpu_jobid="$(sbatch --parsable genebayes_cpu.sh)"
echo "genebayes cpu jobid: ${genebayes_cpu_jobid}"

genebayes_gpu_jobid="$(sbatch --parsable --dependency=afterok:${genebayes_cpu_jobid} genebayes_gpu.sh)"
echo "genebayes gpu jobid: ${genebayes_gpu_jobid}"

essential_dependency="${genebayes_gpu_jobid}"
if [[ -n "${ESSENTIAL_POSTBASE_JOBID:-}" ]]; then
    essential_dependency="${essential_dependency}:${ESSENTIAL_POSTBASE_JOBID}"
else
    echo "ESSENTIAL_POSTBASE_JOBID is not set; assuming essential postbase has already finished."
fi

genomewide_dependency="${genebayes_gpu_jobid}"
if [[ -n "${GENOMEWIDE_POSTBASE_JOBID:-}" ]]; then
    genomewide_dependency="${genomewide_dependency}:${GENOMEWIDE_POSTBASE_JOBID}"
else
    echo "GENOMEWIDE_POSTBASE_JOBID is not set; assuming genomewide postbase has already finished."
fi

essential_assoc_jobid="$(sbatch --parsable --dependency=afterok:${essential_dependency} perturbseq_cnmf_essential_association.sh)"
echo "essential association jobid: ${essential_assoc_jobid}"

genomewide_assoc_jobid="$(sbatch --parsable --dependency=afterok:${genomewide_dependency} perturbseq_cnmf_genomewide_association.sh)"
echo "genomewide association jobid: ${genomewide_assoc_jobid}"
EOF

    chmod +x "$script_path"
}

write_genebayes_packed_script() {
    local script_path="$1"
    local job_name="$2"
    local mem="$3"
    local cpus="$4"
    local partition="$5"
    local time_limit="$6"
    local gres="${7:-}"
    shift 7
    local traits=("$@")
    local script_dir

    script_dir="$(cd "$(dirname "$script_path")" && pwd)"

    cat > "$script_path" <<EOF
#!/usr/bin/env bash
# Run multiple GeneBayes GPU posterior jobs inside one GPU allocation.
# The CPU feature stage must already have completed for every trait config.
# A800 defaults to ${GENEBAYES_PARALLEL_A800} concurrent traits.
# Other GPU nodes default to ${GENEBAYES_PARALLEL_GPU} concurrent traits.
#SBATCH --job-name=${job_name}
#SBATCH --error=${script_dir}/${job_name}.error
#SBATCH --output=${script_dir}/${job_name}.out
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

detect_max_parallel() {
    if [[ "\${SLURM_JOB_PARTITION:-}" == *A800* ]]; then
        echo "${GENEBAYES_PARALLEL_A800}"
        return
    fi
    if command -v nvidia-smi >/dev/null 2>&1; then
        gpu_name="\$(nvidia-smi --query-gpu=name --format=csv,noheader | head -n 1 | tr -d '\r')"
        if [[ "\${gpu_name}" == *A800* ]]; then
            echo "${GENEBAYES_PARALLEL_A800}"
            return
        fi
    fi
    echo "${GENEBAYES_PARALLEL_GPU}"
}

TRAITS=(
EOF

    for trait in "${traits[@]}"; do
        cat >> "$script_path" <<EOF
    "${trait}"
EOF
    done

    cat >> "$script_path" <<EOF
)

CONFIGS=(
EOF

    for trait in "${traits[@]}"; do
        cat >> "$script_path" <<EOF
    "../trait${trait}/${CONFIG_NAME}"
EOF
    done

    cat >> "$script_path" <<EOF
)

MAX_PARALLEL="\$(detect_max_parallel)"
if (( MAX_PARALLEL > \${#TRAITS[@]} )); then
    MAX_PARALLEL=\${#TRAITS[@]}
fi
if (( MAX_PARALLEL < 1 )); then
    MAX_PARALLEL=1
fi

THREADS_PER_TRAIT="${GENEBAYES_THREADS_PER_TRAIT}"
if [[ -z "\${THREADS_PER_TRAIT}" ]]; then
    cpu_budget="\${SLURM_CPUS_PER_TASK:-${cpus}}"
    THREADS_PER_TRAIT=\$(( cpu_budget / MAX_PARALLEL ))
    if (( THREADS_PER_TRAIT < 1 )); then
        THREADS_PER_TRAIT=1
    fi
fi

export OMP_NUM_THREADS="\${THREADS_PER_TRAIT}"
export OPENBLAS_NUM_THREADS="\${THREADS_PER_TRAIT}"
export MKL_NUM_THREADS="\${THREADS_PER_TRAIT}"
export NUMEXPR_NUM_THREADS="\${THREADS_PER_TRAIT}"

echo "Detected max parallel trait jobs on this GPU allocation: \${MAX_PARALLEL}"
echo "Threads per trait: \${THREADS_PER_TRAIT}"

status=0
start_idx=0
total_traits="\${#CONFIGS[@]}"

while (( start_idx < total_traits )); do
    end_idx=\$(( start_idx + MAX_PARALLEL - 1 ))
    if (( end_idx >= total_traits )); then
        end_idx=\$(( total_traits - 1 ))
    fi

    pids=()
    names=()
    for (( idx=start_idx; idx<=end_idx; idx++ )); do
        cfg="\${CONFIGS[\$idx]}"
        trait="\${TRAITS[\$idx]}"
        log_prefix="genebayes_trait\${trait}"
        echo "[\$(date '+%F %T')] START genebayes gpu trait \${trait}"
        (
            set -euo pipefail
            paper-pipeline run --config "\${cfg}" genebayes-gpu
        ) > "\${log_prefix}.log" 2>&1 &
        pids+=("\$!")
        names+=("\${trait}")
    done

    for idx in "\${!pids[@]}"; do
        if wait "\${pids[\$idx]}"; then
            echo "[\$(date '+%F %T')] END   genebayes gpu trait \${names[\$idx]}"
        else
            echo "[\$(date '+%F %T')] FAIL  genebayes gpu trait \${names[\$idx]}" >&2
            status=1
        fi
    done

    start_idx=\$(( end_idx + 1 ))
done

exit "\${status}"
EOF

    chmod +x "$script_path"
}

trait_numbers=()

for num in $(seq "$START_TRAIT" "$END_TRAIT"); do
    batch_dir="${BATCH_ROOT}/trait${num}"
    config_path="${batch_dir}/${CONFIG_NAME}"
    trait_file="${TRAIT_PREFIX}${num}.per_gene_estimates.tsv"

    if [[ -z "${batch_dir}" || "${batch_dir}" == "/" ]]; then
        echo "FATAL: refusing to rm -rf empty or root path" >&2; exit 1
    fi
    rm -rf "$batch_dir"
    mkdir -p "$batch_dir"
    write_trait_config "${config_path}" "${num}" "${trait_file}"

    write_slurm_script \
        "${batch_dir}/genebayes_cpu.sh" \
        "genebayes_cpu_${num}" \
        "${GENEBAYES_CPU_MEM}" \
        "${GENEBAYES_CPU_CPUS}" \
        "${GENEBAYES_CPU_PARTITION}" \
        "${GENEBAYES_CPU_TIME}" \
        "genebayes-cpu"

    write_slurm_script \
        "${batch_dir}/genebayes_gpu.sh" \
        "genebayes_gpu_${num}" \
        "${GENEBAYES_GPU_MEM}" \
        "${GENEBAYES_GPU_CPUS}" \
        "${GENEBAYES_GPU_PARTITION}" \
        "${GENEBAYES_GPU_TIME}" \
        "genebayes-gpu" \
        "${GENEBAYES_GPU_GRES}"

    write_slurm_script \
        "${batch_dir}/perturbseq_cnmf_essential_association.sh" \
        "pertseq_cnmf_essential_association_${num}" \
        "${ESSENTIAL_ASSOC_MEM}" \
        "${ESSENTIAL_ASSOC_CPUS}" \
        "${ESSENTIAL_ASSOC_PARTITION}" \
        "${ESSENTIAL_ASSOC_TIME}" \
        "perturbseq-cnmf-essential-association"

    write_slurm_script \
        "${batch_dir}/perturbseq_cnmf_genomewide_association.sh" \
        "pertseq_cnmf_genomewide_association_${num}" \
        "${GENOMEWIDE_ASSOC_MEM}" \
        "${GENOMEWIDE_ASSOC_CPUS}" \
        "${GENOMEWIDE_ASSOC_PARTITION}" \
        "${GENOMEWIDE_ASSOC_TIME}" \
        "perturbseq-cnmf-genomewide-association"

    write_trait_submitter "${batch_dir}/submit_trait_chain.sh"

    trait_numbers+=("${num}")
done

if (( GENEBAYES_TRAITS_PER_GPU_JOB > 1 )); then
    if [[ -z "${GENEBAYES_PACKED_ROOT}" || "${GENEBAYES_PACKED_ROOT}" == "/" ]]; then
        echo "FATAL: refusing to rm -rf empty or root path" >&2; exit 1
    fi
    rm -rf "${GENEBAYES_PACKED_ROOT}"
    mkdir -p "${GENEBAYES_PACKED_ROOT}"
    group_index=1
    total_traits=${#trait_numbers[@]}

    for (( start_idx=0; start_idx<total_traits; start_idx+=GENEBAYES_TRAITS_PER_GPU_JOB )); do
        group=("${trait_numbers[@]:start_idx:GENEBAYES_TRAITS_PER_GPU_JOB}")
        first_trait="${group[0]}"
        last_trait="${group[${#group[@]}-1]}"
        script_path="${GENEBAYES_PACKED_ROOT}/genebayes_pack_${first_trait}_${last_trait}.sh"
        job_name=$(printf "gb_pack_%03d" "${group_index}")

        write_genebayes_packed_script \
            "${script_path}" \
            "${job_name}" \
            "${GENEBAYES_GPU_MEM}" \
            "${GENEBAYES_GPU_CPUS}" \
            "${GENEBAYES_GPU_PARTITION}" \
            "${GENEBAYES_GPU_TIME}" \
            "${GENEBAYES_GPU_GRES}" \
            "${group[@]}"

        group_index=$(( group_index + 1 ))
    done
fi

echo "Generated per-trait scripts under: ${BATCH_ROOT}/trait*"
if (( GENEBAYES_TRAITS_PER_GPU_JOB > 1 )); then
    echo "Generated packed GeneBayes GPU scripts under: ${GENEBAYES_PACKED_ROOT}"
fi
