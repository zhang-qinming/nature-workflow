#!/usr/bin/env bash

set -euo pipefail

DEFAULT_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_ROOT="${BATCH_ROOT:-${DEFAULT_PROJECT_ROOT}/scripts/kolf}"
SHARED_ROOT="${SHARED_ROOT:-${BATCH_ROOT}/shared}"
JOBS_ROOT="${JOBS_ROOT:-${BATCH_ROOT}/jobs}"
LOGS_ROOT="${LOGS_ROOT:-${BATCH_ROOT}/logs}"
SHARED_CONFIG_NAME="${SHARED_CONFIG_NAME:-pipeline.shared.yaml}"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
if [[ "${PROJECT_ROOT}" != /* ]]; then
    PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
fi
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${PROJECT_ROOT}/outputs}"
if [[ "${ARTIFACT_ROOT}" != /* ]]; then
    ARTIFACT_ROOT="$(cd "${PROJECT_ROOT}" && realpath -m "${ARTIFACT_ROOT}")"
fi

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
PERTURBSEQ_ENV="${PERTURBSEQ_ENV:-paper-pipeline-perturbseq}"
R_ENV="${R_ENV:-paper-pipeline-r}"

KOLF_ROOT="${KOLF_ROOT:-${PROJECT_ROOT}/data/KOLF}"
KOLF_PREPARED_ROOT="${KOLF_PREPARED_ROOT:-${KOLF_ROOT}/pipeline_ready}"
KOLF_INPUT_H5AD="${KOLF_INPUT_H5AD:-${KOLF_ROOT}/KOLF_Pan_Genome_QC_Filtered.h5ad}"
KOLF_CELL_METADATA="${KOLF_CELL_METADATA:-${KOLF_ROOT}/cell_metadata.csv}"
KOLF_GENE_METADATA="${KOLF_GENE_METADATA:-${KOLF_ROOT}/gene_metadata.csv}"
KOLF_PREFIX="${KOLF_PREFIX:-KOLF}"
STANDARDIZE_SCRIPT="${STANDARDIZE_SCRIPT:-${PROJECT_ROOT}/data/KOLF/standardize_kolf_for_pipeline.py}"

KOLF_PIPELINE_H5AD="${KOLF_PIPELINE_H5AD:-${KOLF_PREPARED_ROOT}/${KOLF_PREFIX}.pipeline.h5ad}"
KOLF_PIPELINE_METADATA="${KOLF_PIPELINE_METADATA:-${KOLF_PREPARED_ROOT}/${KOLF_PREFIX}.cell_metadata.pipeline.csv}"
KOLF_GENE_MAP="${KOLF_GENE_MAP:-${KOLF_PREPARED_ROOT}/${KOLF_PREFIX}.gene_map.tsv}"

KOLF_OUTPUT_LABEL="${KOLF_OUTPUT_LABEL:-kolf}"
CNMF_ROOT="${CNMF_ROOT:-${ARTIFACT_ROOT}/perturbseq/${KOLF_OUTPUT_LABEL}/cNMF}"
REGULATION_DIR="${REGULATION_DIR:-${ARTIFACT_ROOT}/perturbseq/${KOLF_OUTPUT_LABEL}/cNMF_regulation}"
REG_SCRIPT="${REG_SCRIPT:-${PROJECT_ROOT}/src/paper_pipeline/dataprep/perturbseq/r_scripts/cnmf_regulatory_effects.R}"

K_VALUES="${K_VALUES:-30 60 90 120 150}"
KSELECT_KS="${KSELECT_KS:-${K_VALUES// /, }}"
KSELECT_TOTAL_WORKERS="${KSELECT_TOTAL_WORKERS:-100}"
KSELECT_WORKERS_PER_JOB="${KSELECT_WORKERS_PER_JOB:-1}"
KSELECT_N_ITER="${KSELECT_N_ITER:-100}"
KSELECT_SEED="${KSELECT_SEED:-14}"
KSELECT_NUMGENES="${KSELECT_NUMGENES:-2000}"
DENSITY_THRESHOLD="${DENSITY_THRESHOLD:-0.5}"
DENSITY_LABEL="${DENSITY_LABEL:-${DENSITY_THRESHOLD//./_}}"
PROGRAMS_PER_JOB="${PROGRAMS_PER_JOB:-5}"

STANDARDIZE_MEM="${STANDARDIZE_MEM:-950G}"
STANDARDIZE_CPUS="${STANDARDIZE_CPUS:-56}"
STANDARDIZE_PARTITION="${STANDARDIZE_PARTITION:-fat}"
STANDARDIZE_TIME="${STANDARDIZE_TIME:-7-00:00:00}"

KSELECT_PREP_MEM="${KSELECT_PREP_MEM:-950G}"
KSELECT_PREP_CPUS="${KSELECT_PREP_CPUS:-56}"
KSELECT_PREP_PARTITION="${KSELECT_PREP_PARTITION:-fat}"
KSELECT_PREP_TIME="${KSELECT_PREP_TIME:-7-00:00:00}"

KSELECT_FACTORIZE_MEM="${KSELECT_FACTORIZE_MEM:-40G}"
KSELECT_FACTORIZE_CPUS="${KSELECT_FACTORIZE_CPUS:-9}"
KSELECT_FACTORIZE_PARTITION="${KSELECT_FACTORIZE_PARTITION:-cu,privority,batch01}"
KSELECT_FACTORIZE_TIME="${KSELECT_FACTORIZE_TIME:-30-00:00:00}"

KSELECT_POST_MEM="${KSELECT_POST_MEM:-950G}"
KSELECT_POST_CPUS="${KSELECT_POST_CPUS:-56}"
KSELECT_POST_PARTITION="${KSELECT_POST_PARTITION:-fat}"
KSELECT_POST_TIME="${KSELECT_POST_TIME:-7-00:00:00}"

CONSENSUS_MEM="${CONSENSUS_MEM:-900G}"
CONSENSUS_CPUS="${CONSENSUS_CPUS:-50}"
CONSENSUS_PARTITION="${CONSENSUS_PARTITION:-fat}"
CONSENSUS_TIME="${CONSENSUS_TIME:-7-00:00:00}"

REGULATION_MEM="${REGULATION_MEM:-60G}"
REGULATION_CPUS="${REGULATION_CPUS:-12}"
REGULATION_PARTITION="${REGULATION_PARTITION:-cu,privority,batch01}"
REGULATION_TIME="${REGULATION_TIME:-7-00:00:00}"

for value_name in KSELECT_TOTAL_WORKERS KSELECT_WORKERS_PER_JOB KSELECT_N_ITER KSELECT_NUMGENES PROGRAMS_PER_JOB; do
    value="${!value_name}"
    if [[ "$value" =~ [^0-9] ]] || (( value < 1 )); then
        echo "${value_name} must be a positive integer" >&2
        exit 1
    fi
done

mkdir -p "${BATCH_ROOT}" "${SHARED_ROOT}" "${JOBS_ROOT}" "${LOGS_ROOT}"
BATCH_ROOT="$(cd "${BATCH_ROOT}" && pwd)"
SHARED_ROOT="$(cd "${SHARED_ROOT}" && pwd)"
JOBS_ROOT="$(cd "${JOBS_ROOT}" && pwd)"
LOGS_ROOT="$(cd "${LOGS_ROOT}" && pwd)"

for path_name in SHARED_ROOT JOBS_ROOT; do
    path_value="${!path_name}"
    if [[ -z "${path_value}" || "${path_value}" == "/" ]]; then
        echo "FATAL: refusing to remove empty or root ${path_name}" >&2
        exit 1
    fi
done
rm -rf "${SHARED_ROOT}" "${JOBS_ROOT}"
mkdir -p "${SHARED_ROOT}" "${JOBS_ROOT}"

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

write_pipeline_job() {
    local script_path="$1"
    local job_name="$2"
    local mem="$3"
    local cpus="$4"
    local partition="$5"
    local time_limit="$6"
    local workflow_name="$7"
    local config_file="${8:-${SHARED_CONFIG_NAME}}"
    local config_path

    if [[ "${config_file}" = /* ]]; then
        config_path="${config_file}"
    else
        config_path="${SHARED_ROOT}/${config_file}"
    fi

    cat > "${script_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${LOGS_ROOT}/${job_name}_%j.out
#SBATCH --error=${LOGS_ROOT}/${job_name}_%j.err
#SBATCH --mem=${mem}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --partition=${partition}
#SBATCH --time=${time_limit}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"

paper-pipeline run --config "${config_path}" ${workflow_name}
EOF
    chmod +x "${script_path}"
}

write_kselect_config() {
    local config_path="$1"
    local worker_indices="${2:-}"

    cat > "${config_path}" <<EOF
project_root: ${PROJECT_ROOT}
artifact_root: ${ARTIFACT_ROOT}

executables:
  perturbseq_runner:
    - conda
    - run
    - --no-capture-output
    - -n
    - ${PERTURBSEQ_ENV}

workflows:
  perturbseq:
    cnmf_essential_kselect:
      inputs:
        h5ad: ${KOLF_PIPELINE_H5AD}
      outputs:
        cnmf_root: ${CNMF_ROOT}
      parameters:
        ks: [${KSELECT_KS}]
        total_workers: ${KSELECT_TOTAL_WORKERS}
EOF

    if [[ -n "${worker_indices}" ]]; then
        cat >> "${config_path}" <<EOF
        worker_indices: [${worker_indices}]
EOF
    fi

    cat >> "${config_path}" <<EOF
        n_iter: ${KSELECT_N_ITER}
        seed: ${KSELECT_SEED}
        numgenes: ${KSELECT_NUMGENES}
        aggregate_name: cNMF_all
        per_k_name_template: cNMF_K{k}
EOF
}

write_standardize_job() {
    local script_path="${SHARED_ROOT}/kolf_standardize.sh"
    cat > "${script_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=kolf_standardize
#SBATCH --output=${LOGS_ROOT}/kolf_standardize_%j.out
#SBATCH --error=${LOGS_ROOT}/kolf_standardize_%j.err
#SBATCH --mem=${STANDARDIZE_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${STANDARDIZE_CPUS}
#SBATCH --partition=${STANDARDIZE_PARTITION}
#SBATCH --time=${STANDARDIZE_TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

source "${CONDA_SH}"
conda activate "${PERTURBSEQ_ENV}"

python "${STANDARDIZE_SCRIPT}" \\
  --input-h5ad "${KOLF_INPUT_H5AD}" \\
  --cell-metadata "${KOLF_CELL_METADATA}" \\
  --gene-metadata "${KOLF_GENE_METADATA}" \\
  --output-dir "${KOLF_PREPARED_ROOT}" \\
  --prefix "${KOLF_PREFIX}" \\
  --write-h5ad
EOF
    chmod +x "${script_path}"
}

write_consensus_job() {
    local k="$1"
    local job_path="${JOBS_ROOT}/consensus_K${k}.sh"
    cat > "${job_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=kolf_cons_K${k}
#SBATCH --output=${LOGS_ROOT}/consensus_K${k}_%j.out
#SBATCH --error=${LOGS_ROOT}/consensus_K${k}_%j.err
#SBATCH --mem=${CONSENSUS_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CONSENSUS_CPUS}
#SBATCH --partition=${CONSENSUS_PARTITION}
#SBATCH --time=${CONSENSUS_TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

source "${CONDA_SH}"

conda run --no-capture-output -n "${PERTURBSEQ_ENV}" \\
  cnmf consensus \\
  --output-dir "${CNMF_ROOT}" \\
  --name cNMF_all \\
  --components "${k}" \\
  --local-density-threshold "${DENSITY_THRESHOLD}" \\
  --show-clustering
EOF
    chmod +x "${job_path}"
}

write_regulation_job() {
    local k="$1"
    local start="$2"
    local end="$3"
    local job_path="${JOBS_ROOT}/regulation_K${k}_program_${start}_${end}.sh"
    cat > "${job_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=kolf_reg_K${k}_${start}_${end}
#SBATCH --output=${LOGS_ROOT}/regulation_K${k}_${start}_${end}_%j.out
#SBATCH --error=${LOGS_ROOT}/regulation_K${k}_${start}_${end}_%j.err
#SBATCH --mem=${REGULATION_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${REGULATION_CPUS}
#SBATCH --partition=${REGULATION_PARTITION}
#SBATCH --time=${REGULATION_TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

source "${CONDA_SH}"
mkdir -p "${REGULATION_DIR}"

USAGES_PATH="${CNMF_ROOT}/cNMF_all/cNMF_all.usages.k_${k}.dt_${DENSITY_LABEL}.consensus.txt"
if [[ ! -s "\${USAGES_PATH}" ]]; then
    echo "Missing consensus usages file: \${USAGES_PATH}" >&2
    exit 1
fi

for program_index in \$(seq ${start} ${end}); do
  conda run --no-capture-output -n "${R_ENV}" Rscript \\
    "${REG_SCRIPT}" \\
    "${KOLF_PIPELINE_METADATA}" \\
    "\${USAGES_PATH}" \\
    "\${program_index}" \\
    "${REGULATION_DIR}/K${k}_program\${program_index}_perturb_effects.txt"
done
EOF
    chmod +x "${job_path}"
}

write_kselect_config "${SHARED_ROOT}/${SHARED_CONFIG_NAME}"
write_standardize_job

write_pipeline_job \
    "${SHARED_ROOT}/kolf_cnmf_kselect_prepare.sh" \
    "kolf_kselect_prepare" \
    "${KSELECT_PREP_MEM}" \
    "${KSELECT_PREP_CPUS}" \
    "${KSELECT_PREP_PARTITION}" \
    "${KSELECT_PREP_TIME}" \
    "perturbseq-cnmf-essential-kselect-prepare"

write_pipeline_job \
    "${SHARED_ROOT}/kolf_cnmf_kselect_postbase_base.sh" \
    "kolf_kselect_postbase" \
    "${KSELECT_POST_MEM}" \
    "${KSELECT_POST_CPUS}" \
    "${KSELECT_POST_PARTITION}" \
    "${KSELECT_POST_TIME}" \
    "perturbseq-cnmf-essential-kselect-postbase-base"

for (( worker_start=1; worker_start<=KSELECT_TOTAL_WORKERS; worker_start+=KSELECT_WORKERS_PER_JOB )); do
    worker_indices="$(indices_yaml "${worker_start}" "${KSELECT_WORKERS_PER_JOB}" "${KSELECT_TOTAL_WORKERS}")"
    worker_last=$(( worker_start + KSELECT_WORKERS_PER_JOB - 1 ))
    if (( worker_last > KSELECT_TOTAL_WORKERS )); then
        worker_last="${KSELECT_TOTAL_WORKERS}"
    fi
    worker_cfg="pipeline.kolf_kselect.factorize_${worker_start}_${worker_last}.yaml"
    write_kselect_config "${SHARED_ROOT}/${worker_cfg}" "${worker_indices}"
    write_pipeline_job \
        "${SHARED_ROOT}/kolf_cnmf_kselect_factorize_${worker_start}_${worker_last}.sh" \
        "kolf_ks_fact_${worker_start}_${worker_last}" \
        "${KSELECT_FACTORIZE_MEM}" \
        "${KSELECT_FACTORIZE_CPUS}" \
        "${KSELECT_FACTORIZE_PARTITION}" \
        "${KSELECT_FACTORIZE_TIME}" \
        "perturbseq-cnmf-essential-kselect-factorize" \
        "${worker_cfg}"
done

for k in ${K_VALUES}; do
    if [[ "$k" =~ [^0-9] ]] || (( k < 1 )); then
        echo "Invalid K value: ${k}" >&2
        exit 1
    fi
    write_consensus_job "${k}"
    for (( start=1; start<=k; start+=PROGRAMS_PER_JOB )); do
        end=$(( start + PROGRAMS_PER_JOB - 1 ))
        if (( end > k )); then
            end="${k}"
        fi
        write_regulation_job "${k}" "${start}" "${end}"
    done
done

cat > "${BATCH_ROOT}/submit_full_chain.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_DIR="${SCRIPT_DIR}/shared"
JOBS_DIR="${SCRIPT_DIR}/jobs"

standardize_jobid="$(sbatch --parsable "${SHARED_DIR}/kolf_standardize.sh")"
echo "standardize jobid: ${standardize_jobid}"

prepare_jobid="$(sbatch --parsable --dependency=afterok:${standardize_jobid} "${SHARED_DIR}/kolf_cnmf_kselect_prepare.sh")"
echo "prepare jobid: ${prepare_jobid}"

factorize_jobids=()
shopt -s nullglob
for script in "${SHARED_DIR}"/kolf_cnmf_kselect_factorize_*.sh; do
    jobid="$(sbatch --parsable --dependency=afterok:${prepare_jobid} "${script}")"
    factorize_jobids+=("${jobid}")
    echo "factorize jobid: ${jobid} ($(basename "${script}"))"
done

if (( ${#factorize_jobids[@]} == 0 )); then
    echo "No factorize scripts found" >&2
    exit 1
fi

factorize_dependency="$(IFS=:; echo "${factorize_jobids[*]}")"
postbase_jobid="$(sbatch --parsable --dependency=afterok:${factorize_dependency} "${SHARED_DIR}/kolf_cnmf_kselect_postbase_base.sh")"
echo "postbase jobid: ${postbase_jobid}"

for consensus_script in "${JOBS_DIR}"/consensus_K*.sh; do
    base="$(basename "${consensus_script}")"
    k="${base#consensus_K}"
    k="${k%.sh}"
    consensus_jobid="$(sbatch --parsable --dependency=afterok:${postbase_jobid} "${consensus_script}")"
    echo "consensus K=${k}: ${consensus_jobid}"

    for regulation_script in "${JOBS_DIR}"/regulation_K${k}_program_*.sh; do
        regulation_jobid="$(sbatch --parsable --dependency=afterok:${consensus_jobid} "${regulation_script}")"
        echo "  regulation $(basename "${regulation_script}"): ${regulation_jobid}"
    done
done
shopt -u nullglob
EOF
chmod +x "${BATCH_ROOT}/submit_full_chain.sh"

cat > "${BATCH_ROOT}/submit_cnmf_regulation_chain.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_DIR="${SCRIPT_DIR}/shared"
JOBS_DIR="${SCRIPT_DIR}/jobs"

prepare_jobid="$(sbatch --parsable "${SHARED_DIR}/kolf_cnmf_kselect_prepare.sh")"
echo "prepare jobid: ${prepare_jobid}"

factorize_jobids=()
shopt -s nullglob
for script in "${SHARED_DIR}"/kolf_cnmf_kselect_factorize_*.sh; do
    jobid="$(sbatch --parsable --dependency=afterok:${prepare_jobid} "${script}")"
    factorize_jobids+=("${jobid}")
    echo "factorize jobid: ${jobid} ($(basename "${script}"))"
done

if (( ${#factorize_jobids[@]} == 0 )); then
    echo "No factorize scripts found" >&2
    exit 1
fi

factorize_dependency="$(IFS=:; echo "${factorize_jobids[*]}")"
postbase_jobid="$(sbatch --parsable --dependency=afterok:${factorize_dependency} "${SHARED_DIR}/kolf_cnmf_kselect_postbase_base.sh")"
echo "postbase jobid: ${postbase_jobid}"

for consensus_script in "${JOBS_DIR}"/consensus_K*.sh; do
    base="$(basename "${consensus_script}")"
    k="${base#consensus_K}"
    k="${k%.sh}"
    consensus_jobid="$(sbatch --parsable --dependency=afterok:${postbase_jobid} "${consensus_script}")"
    echo "consensus K=${k}: ${consensus_jobid}"

    for regulation_script in "${JOBS_DIR}"/regulation_K${k}_program_*.sh; do
        regulation_jobid="$(sbatch --parsable --dependency=afterok:${consensus_jobid} "${regulation_script}")"
        echo "  regulation $(basename "${regulation_script}"): ${regulation_jobid}"
    done
done
shopt -u nullglob
EOF
chmod +x "${BATCH_ROOT}/submit_cnmf_regulation_chain.sh"

cat <<EOF
Generated KOLF batch scripts under: ${BATCH_ROOT}
K values: ${K_VALUES}
Shared config:
  ${SHARED_ROOT}/${SHARED_CONFIG_NAME}
Expected standardized h5ad:
  ${KOLF_PIPELINE_H5AD}
Expected K-selection plot:
  ${CNMF_ROOT}/cNMF_all/cNMF_all.k_selection.png
Expected regulation directory:
  ${REGULATION_DIR}

Submit full chain with:
  cd ${BATCH_ROOT}
  ./submit_full_chain.sh

If ${KOLF_PIPELINE_H5AD} already exists, skip standardization with:
  cd ${BATCH_ROOT}
  ./submit_cnmf_regulation_chain.sh
EOF
