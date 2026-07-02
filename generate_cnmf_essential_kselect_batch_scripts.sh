#!/usr/bin/env bash

set -euo pipefail

DEFAULT_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_ROOT="${BATCH_ROOT:-${DEFAULT_PROJECT_ROOT}/scripts/cnmf_essential_kselect}"
SHARED_ROOT="${SHARED_ROOT:-${BATCH_ROOT}/shared}"
LOGS_ROOT="${LOGS_ROOT:-${BATCH_ROOT}/logs}"
SHARED_CONFIG_NAME="${SHARED_CONFIG_NAME:-pipeline.shared.yaml}"

CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
PERTURBSEQ_ENV="${PERTURBSEQ_ENV:-paper-pipeline-perturbseq}"
R_ENV="${R_ENV:-paper-pipeline-r}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"

KSELECT_CELL="${KSELECT_CELL:-K562_essential_raw_singlecell_01}"
KSELECT_KS="${KSELECT_KS:-30, 40, 50, 60, 70, 80, 90, 120}"
KSELECT_TOTAL_WORKERS="${KSELECT_TOTAL_WORKERS:-20}"
KSELECT_WORKERS_PER_JOB="${KSELECT_WORKERS_PER_JOB:-2}"
KSELECT_N_ITER="${KSELECT_N_ITER:-100}"
KSELECT_SEED="${KSELECT_SEED:-14}"
KSELECT_NUMGENES="${KSELECT_NUMGENES:-2000}"
KSELECT_OUTPUT_LABEL="${KSELECT_OUTPUT_LABEL:-cnmf_essential_kselect}"

KSELECT_PREP_MEM="${KSELECT_PREP_MEM:-450G}"
KSELECT_PREP_CPUS="${KSELECT_PREP_CPUS:-10}"
KSELECT_PREP_PARTITION="${KSELECT_PREP_PARTITION:-fat}"
KSELECT_PREP_TIME="${KSELECT_PREP_TIME:-7-00:00:00}"

KSELECT_FACTORIZE_MEM="${KSELECT_FACTORIZE_MEM:-80G}"
KSELECT_FACTORIZE_CPUS="${KSELECT_FACTORIZE_CPUS:-16}"
KSELECT_FACTORIZE_PARTITION="${KSELECT_FACTORIZE_PARTITION:-cu,privority,batch01}"
KSELECT_FACTORIZE_TIME="${KSELECT_FACTORIZE_TIME:-10-00:00:00}"

KSELECT_POST_MEM="${KSELECT_POST_MEM:-60G}"
KSELECT_POST_CPUS="${KSELECT_POST_CPUS:-15}"
KSELECT_POST_PARTITION="${KSELECT_POST_PARTITION:-cu,privority,batch01}"
KSELECT_POST_TIME="${KSELECT_POST_TIME:-7-00:00:00}"

for value_name in KSELECT_TOTAL_WORKERS KSELECT_WORKERS_PER_JOB KSELECT_N_ITER KSELECT_NUMGENES; do
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
elif [[ "${PROJECT_ROOT}" != /* ]]; then
    PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
fi
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${PROJECT_ROOT}/outputs}"
if [[ "${ARTIFACT_ROOT}" != /* ]]; then
    ARTIFACT_ROOT="$(cd "${PROJECT_ROOT}" && realpath -m "${ARTIFACT_ROOT}")"
fi

mkdir -p "${SHARED_ROOT}" "${LOGS_ROOT}"
SHARED_ROOT="$(cd "${SHARED_ROOT}" && pwd)"
LOGS_ROOT="$(cd "${LOGS_ROOT}" && pwd)"
if [[ -z "${SHARED_ROOT}" || "${SHARED_ROOT}" == "/" ]]; then
    echo "FATAL: refusing to rm -rf empty or root path" >&2
    exit 1
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

indices_yaml() {
    local start="$1"
    local per_job="$2"
    local total="$3"
    local values=()
    local end=$(( start + per_job - 1 ))
    if (( end > total )); then
        end=$total
    fi
    for (( idx=start; idx<=end; idx++ )); do
        values+=("$(( idx - 1 ))")
    done
    local joined
    printf -v joined '%s, ' "${values[@]}"
    joined="${joined%, }"
    printf '%s' "${joined}"
}

write_slurm_script() {
    local script_path="$1"
    local job_name="$2"
    local mem="$3"
    local cpus="$4"
    local partition="$5"
    local time_limit="$6"
    local workflow_name="$7"
    local config_file="${8:-${SHARED_CONFIG_NAME}}"
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

write_kselect_config() {
    local config_path="$1"
    local worker_indices="${2:-}"

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
    cnmf_essential_kselect:
      inputs:
        h5ad: ${ARTIFACT_ROOT}/perturbseq/cnmf_essential/filtered_data/${KSELECT_CELL}.h5ad
      outputs:
        cnmf_root: ${ARTIFACT_ROOT}/perturbseq/${KSELECT_OUTPUT_LABEL}/cNMF
      parameters:
        ks: [${KSELECT_KS}]
        total_workers: ${KSELECT_TOTAL_WORKERS}
EOF

    if [[ -n "${worker_indices}" ]]; then
        cat >> "$config_path" <<EOF
        worker_indices: [${worker_indices}]
EOF
    fi

    cat >> "$config_path" <<EOF
        n_iter: ${KSELECT_N_ITER}
        seed: ${KSELECT_SEED}
        numgenes: ${KSELECT_NUMGENES}
        aggregate_name: cNMF_all
        per_k_name_template: cNMF_K{k}
EOF
}

write_submitter() {
    local script_path="$1"

    cat > "$script_path" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/shared"

prepare_jobid="$(sbatch --parsable perturbseq_cnmf_essential_kselect_prepare.sh)"
echo "prepare jobid: ${prepare_jobid}"

factorize_jobids=()
shopt -s nullglob
for script in perturbseq_cnmf_essential_kselect_factorize_*.sh; do
    jobid="$(sbatch --parsable --dependency=afterok:${prepare_jobid} "${script}")"
    factorize_jobids+=("${jobid}")
    echo "factorize jobid: ${jobid} (${script})"
done
shopt -u nullglob

if (( ${#factorize_jobids[@]} == 0 )); then
    echo "No factorize scripts found" >&2
    exit 1
fi

dependency="$(IFS=:; echo "${factorize_jobids[*]}")"
postbase_jobid="$(sbatch --parsable --dependency=afterok:${dependency} perturbseq_cnmf_essential_kselect_postbase_base.sh)"
echo "postbase base jobid: ${postbase_jobid}"
EOF

    chmod +x "$script_path"
}

write_kselect_config "${SHARED_ROOT}/${SHARED_CONFIG_NAME}"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_essential_kselect_prepare.sh" \
    "pertseq_cnmf_ess_kselect_prep" \
    "${KSELECT_PREP_MEM}" \
    "${KSELECT_PREP_CPUS}" \
    "${KSELECT_PREP_PARTITION}" \
    "${KSELECT_PREP_TIME}" \
    "perturbseq-cnmf-essential-kselect-prepare"

write_slurm_script \
    "${SHARED_ROOT}/perturbseq_cnmf_essential_kselect_postbase_base.sh" \
    "pertseq_cnmf_ess_kselect_post" \
    "${KSELECT_POST_MEM}" \
    "${KSELECT_POST_CPUS}" \
    "${KSELECT_POST_PARTITION}" \
    "${KSELECT_POST_TIME}" \
    "perturbseq-cnmf-essential-kselect-postbase-base"

for (( worker_start=1; worker_start<=KSELECT_TOTAL_WORKERS; worker_start+=KSELECT_WORKERS_PER_JOB )); do
    worker_indices="$(indices_yaml "${worker_start}" "${KSELECT_WORKERS_PER_JOB}" "${KSELECT_TOTAL_WORKERS}")"
    worker_last=$(( worker_start + KSELECT_WORKERS_PER_JOB - 1 ))
    if (( worker_last > KSELECT_TOTAL_WORKERS )); then
        worker_last=${KSELECT_TOTAL_WORKERS}
    fi
    worker_cfg="pipeline.essential_kselect.factorize_${worker_start}_${worker_last}.yaml"
    worker_script="perturbseq_cnmf_essential_kselect_factorize_${worker_start}_${worker_last}.sh"
    write_kselect_config "${SHARED_ROOT}/${worker_cfg}" "${worker_indices}"
    write_slurm_script \
        "${SHARED_ROOT}/${worker_script}" \
        "pertseq_cnmf_ess_ks_fact_${worker_start}_${worker_last}" \
        "${KSELECT_FACTORIZE_MEM}" \
        "${KSELECT_FACTORIZE_CPUS}" \
        "${KSELECT_FACTORIZE_PARTITION}" \
        "${KSELECT_FACTORIZE_TIME}" \
        "perturbseq-cnmf-essential-kselect-factorize" \
        "${worker_cfg}"
done

write_submitter "${BATCH_ROOT}/submit_essential_kselect_chain.sh"

cat <<EOF
Generated essential cNMF K-selection scripts under: ${BATCH_ROOT}
Submit with:
  cd ${BATCH_ROOT}
  ./submit_essential_kselect_chain.sh

Expected plot:
  ${ARTIFACT_ROOT}/perturbseq/${KSELECT_OUTPUT_LABEL}/cNMF/cNMF_all/cNMF_all.k_selection.png
EOF
