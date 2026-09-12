#!/usr/bin/env bash

set -euo pipefail

# Generate local Slurm wrappers only; this script never submits jobs.
DEFAULT_PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BATCH_ROOT="${BATCH_ROOT:-${DEFAULT_PROJECT_ROOT}/scripts/cd4}"
SHARED_ROOT="${SHARED_ROOT:-${BATCH_ROOT}/shared}"
JOBS_ROOT="${JOBS_ROOT:-${BATCH_ROOT}/jobs}"
LOGS_ROOT="${LOGS_ROOT:-${BATCH_ROOT}/logs}"
SHARED_CONFIG_NAME="${SHARED_CONFIG_NAME:-pipeline.shared.yaml}"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
if [[ "${PROJECT_ROOT}" != /* ]]; then PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"; fi
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${PROJECT_ROOT}/outputs}"
if [[ "${ARTIFACT_ROOT}" != /* ]]; then ARTIFACT_ROOT="$(cd "${PROJECT_ROOT}" && realpath -m "${ARTIFACT_ROOT}")"; fi

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
PERTURBSEQ_ENV="${PERTURBSEQ_ENV:-paper-pipeline-perturbseq}"
R_ENV="${R_ENV:-paper-pipeline-r}"

CD4_ROOT="${CD4_ROOT:-${PROJECT_ROOT}/data/CD4}"
CD4_INPUT_H5AD="${CD4_INPUT_H5AD:-${CD4_ROOT}/GWCD4_Rest.filtered_merged.h5ad}"
CD4_PREPARED_ROOT="${CD4_PREPARED_ROOT:-${CD4_ROOT}/pipeline_ready}"
CD4_PREFIX="${CD4_PREFIX:-CD4}"
STANDARDIZE_SCRIPT="${STANDARDIZE_SCRIPT:-${CD4_ROOT}/prepare_cd4_for_pipeline.py}"
CD4_METADATA="${CD4_METADATA:-${CD4_PREPARED_ROOT}/${CD4_PREFIX}.cell_metadata.pipeline.csv}"
CD4_GENE_MAP="${CD4_GENE_MAP:-${CD4_PREPARED_ROOT}/${CD4_PREFIX}.gene_map.tsv}"
CD4_REPORT="${CD4_REPORT:-${CD4_PREPARED_ROOT}/${CD4_PREFIX}.prepare_report.json}"

OUTPUT_LABEL="${OUTPUT_LABEL:-cd4}"
CNMF_ROOT="${CNMF_ROOT:-${ARTIFACT_ROOT}/perturbseq/${OUTPUT_LABEL}/cNMF}"
REGULATION_DIR="${REGULATION_DIR:-${ARTIFACT_ROOT}/perturbseq/${OUTPUT_LABEL}/cNMF_regulation}"
REG_SCRIPT="${REG_SCRIPT:-${PROJECT_ROOT}/src/paper_pipeline/dataprep/perturbseq/r_scripts/cnmf_regulatory_effects.R}"

K_VALUES="${K_VALUES:-30 60 90 120 150}"
KSELECT_KS="${KSELECT_KS:-${K_VALUES// /, }}"
TOTAL_WORKERS="${KSELECT_TOTAL_WORKERS:-100}"
WORKERS_PER_JOB="${KSELECT_WORKERS_PER_JOB:-1}"
N_ITER="${KSELECT_N_ITER:-100}"
SEED="${KSELECT_SEED:-14}"
NUMGENES="${KSELECT_NUMGENES:-2000}"
DENSITY_THRESHOLD="${DENSITY_THRESHOLD:-0.5}"
DENSITY_LABEL="${DENSITY_LABEL:-${DENSITY_THRESHOLD//./_}}"
PROGRAMS_PER_JOB="${PROGRAMS_PER_JOB:-5}"

STD_MEM="${STANDARDIZE_MEM:-128G}"
STD_CPUS="${STANDARDIZE_CPUS:-8}"
STD_PARTITION="${STANDARDIZE_PARTITION:-fat}"
STD_TIME="${STANDARDIZE_TIME:-2-00:00:00}"
PREP_MEM="${KSELECT_PREP_MEM:-950G}"
PREP_CPUS="${KSELECT_PREP_CPUS:-56}"
PREP_PARTITION="${KSELECT_PREP_PARTITION:-fat}"
PREP_TIME="${KSELECT_PREP_TIME:-14-00:00:00}"
FACT_MEM="${KSELECT_FACTORIZE_MEM:-96G}"
FACT_CPUS="${KSELECT_FACTORIZE_CPUS:-16}"
FACT_PARTITION="${KSELECT_FACTORIZE_PARTITION:-cu,privority,batch01}"
FACT_TIME="${KSELECT_FACTORIZE_TIME:-30-00:00:00}"
POST_MEM="${KSELECT_POST_MEM:-950G}"
POST_CPUS="${KSELECT_POST_CPUS:-56}"
POST_PARTITION="${KSELECT_POST_PARTITION:-fat}"
POST_TIME="${KSELECT_POST_TIME:-14-00:00:00}"
CONS_MEM="${CONSENSUS_MEM:-950G}"
CONS_CPUS="${CONSENSUS_CPUS:-56}"
CONS_PARTITION="${CONSENSUS_PARTITION:-fat}"
CONS_TIME="${CONSENSUS_TIME:-14-00:00:00}"
REG_MEM="${REGULATION_MEM:-80G}"
REG_CPUS="${REGULATION_CPUS:-12}"
REG_PARTITION="${REGULATION_PARTITION:-cu,privority,batch01}"
REG_TIME="${REGULATION_TIME:-14-00:00:00}"

for name in TOTAL_WORKERS WORKERS_PER_JOB N_ITER NUMGENES PROGRAMS_PER_JOB; do
    value="${!name}"
    if [[ "${value}" =~ [^0-9] ]] || (( value < 1 )); then
        echo "${name} must be a positive integer" >&2; exit 1
    fi
done
if (( TOTAL_WORKERS % WORKERS_PER_JOB != 0 )); then
    echo "KSELECT_WORKERS_PER_JOB must divide KSELECT_TOTAL_WORKERS" >&2; exit 1
fi
for k in ${K_VALUES}; do
    if [[ "${k}" =~ [^0-9] ]] || (( k < 1 )); then echo "Invalid K: ${k}" >&2; exit 1; fi
done

mkdir -p "${BATCH_ROOT}" "${SHARED_ROOT}" "${JOBS_ROOT}" "${LOGS_ROOT}"
BATCH_ROOT="$(cd "${BATCH_ROOT}" && pwd)"
SHARED_ROOT="$(cd "${SHARED_ROOT}" && pwd)"
JOBS_ROOT="$(cd "${JOBS_ROOT}" && pwd)"
LOGS_ROOT="$(cd "${LOGS_ROOT}" && pwd)"
for path_name in BATCH_ROOT PROJECT_ROOT CD4_ROOT; do
    path_value="${!path_name}"
    if [[ "${path_value}" == "/" ]]; then
        echo "FATAL: refusing to use root as ${path_name}" >&2
        exit 1
    fi
done
if [[ "${BATCH_ROOT}" == "${PROJECT_ROOT}" || "${BATCH_ROOT}" == "${CD4_ROOT}" ]]; then
    echo "FATAL: BATCH_ROOT must be a dedicated batch-script directory" >&2
    echo "       BATCH_ROOT=${BATCH_ROOT}" >&2
    exit 1
fi
for name in SHARED_ROOT JOBS_ROOT; do
    value="${!name}"
    if [[ -z "${value}" || "${value}" == "/" || "${value}" == "${PROJECT_ROOT}" || "${value}" == "${CD4_ROOT}" ]]; then
        echo "Refusing to remove ${name}: ${value}" >&2
        exit 1
    fi
done
rm -rf "${SHARED_ROOT}" "${JOBS_ROOT}"
mkdir -p "${SHARED_ROOT}" "${JOBS_ROOT}"

indices_yaml() {
    local start="$1" per_job="$2" total="$3" end values joined
    end=$((start + per_job - 1)); (( end > total )) && end="$total"
    values=()
    for ((i=start; i<=end; i++)); do values+=("$((i - 1))"); done
    printf -v joined '%s, ' "${values[@]}"
    printf '%s' "${joined%, }"
}

write_pipeline_job() {
    local path="$1" name="$2" mem="$3" cpus="$4" partition="$5" time="$6" workflow="$7" cfg="${8:-${SHARED_CONFIG_NAME}}" cfg_path
    [[ "${cfg}" = /* ]] && cfg_path="${cfg}" || cfg_path="${SHARED_ROOT}/${cfg}"
    cat > "${path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${name}
#SBATCH --output=${LOGS_ROOT}/${name}_%j.out
#SBATCH --error=${LOGS_ROOT}/${name}_%j.err
#SBATCH --mem=${mem}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --partition=${partition}
#SBATCH --time=${time}
#SBATCH --chdir=${PROJECT_ROOT}
set -euo pipefail
source "${CONDA_SH}"
conda activate "${CONTROL_ENV}"
paper-pipeline run --config "${cfg_path}" ${workflow}
EOF
    chmod +x "${path}"
}

write_config() {
    local path="$1" workers="${2:-}"
    cat > "${path}" <<EOF
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
        h5ad: ${CD4_INPUT_H5AD}
      outputs:
        cnmf_root: ${CNMF_ROOT}
      parameters:
        ks: [${KSELECT_KS}]
        total_workers: ${TOTAL_WORKERS}
EOF
    if [[ -n "${workers}" ]]; then printf '        worker_indices: [%s]\n' "${workers}" >> "${path}"; fi
    cat >> "${path}" <<EOF
        n_iter: ${N_ITER}
        seed: ${SEED}
        numgenes: ${NUMGENES}
        aggregate_name: cNMF_all
        per_k_name_template: cNMF_K{k}
EOF
}

cat > "${SHARED_ROOT}/cd4_prepare_metadata.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=cd4_prepare_metadata
#SBATCH --output=${LOGS_ROOT}/cd4_prepare_metadata_%j.out
#SBATCH --error=${LOGS_ROOT}/cd4_prepare_metadata_%j.err
#SBATCH --mem=${STD_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${STD_CPUS}
#SBATCH --partition=${STD_PARTITION}
#SBATCH --time=${STD_TIME}
#SBATCH --chdir=${PROJECT_ROOT}
set -euo pipefail
source "${CONDA_SH}"
conda activate "${PERTURBSEQ_ENV}"
if [[ -s "${CD4_METADATA}" && -s "${CD4_GENE_MAP}" && -s "${CD4_REPORT}" ]]; then
    echo "CD4 pipeline metadata already exists; skipping conversion."
    exit 0
fi
python "${STANDARDIZE_SCRIPT}" --input-h5ad "${CD4_INPUT_H5AD}" --output-dir "${CD4_PREPARED_ROOT}" --prefix "${CD4_PREFIX}" --force
EOF
chmod +x "${SHARED_ROOT}/cd4_prepare_metadata.sh"

write_config "${SHARED_ROOT}/${SHARED_CONFIG_NAME}"
write_pipeline_job "${SHARED_ROOT}/cd4_cnmf_kselect_prepare.sh" cd4_kselect_prepare "${PREP_MEM}" "${PREP_CPUS}" "${PREP_PARTITION}" "${PREP_TIME}" perturbseq-cnmf-essential-kselect-prepare
write_pipeline_job "${SHARED_ROOT}/cd4_cnmf_kselect_postbase_base.sh" cd4_kselect_postbase "${POST_MEM}" "${POST_CPUS}" "${POST_PARTITION}" "${POST_TIME}" perturbseq-cnmf-essential-kselect-postbase-base

for ((start=1; start<=TOTAL_WORKERS; start+=WORKERS_PER_JOB)); do
    workers="$(indices_yaml "${start}" "${WORKERS_PER_JOB}" "${TOTAL_WORKERS}")"
    last=$((start + WORKERS_PER_JOB - 1)); (( last > TOTAL_WORKERS )) && last="$TOTAL_WORKERS"
    cfg="pipeline.cd4_kselect.factorize_${start}_${last}.yaml"
    write_config "${SHARED_ROOT}/${cfg}" "${workers}"
    write_pipeline_job "${SHARED_ROOT}/cd4_cnmf_kselect_factorize_${start}_${last}.sh" "cd4_ks_fact_${start}_${last}" "${FACT_MEM}" "${FACT_CPUS}" "${FACT_PARTITION}" "${FACT_TIME}" perturbseq-cnmf-essential-kselect-factorize "${cfg}"
done

for k in ${K_VALUES}; do
    cat > "${JOBS_ROOT}/consensus_K${k}.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=cd4_cons_K${k}
#SBATCH --output=${LOGS_ROOT}/consensus_K${k}_%j.out
#SBATCH --error=${LOGS_ROOT}/consensus_K${k}_%j.err
#SBATCH --mem=${CONS_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CONS_CPUS}
#SBATCH --partition=${CONS_PARTITION}
#SBATCH --time=${CONS_TIME}
#SBATCH --chdir=${PROJECT_ROOT}
set -euo pipefail
source "${CONDA_SH}"
conda run --no-capture-output -n "${PERTURBSEQ_ENV}" cnmf consensus --output-dir "${CNMF_ROOT}" --name cNMF_all --components "${k}" --local-density-threshold "${DENSITY_THRESHOLD}" --show-clustering
EOF
    chmod +x "${JOBS_ROOT}/consensus_K${k}.sh"

    for ((p=1; p<=k; p+=PROGRAMS_PER_JOB)); do
        end=$((p + PROGRAMS_PER_JOB - 1)); (( end > k )) && end="$k"
        cat > "${JOBS_ROOT}/regulation_K${k}_program_${p}_${end}.sh" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=cd4_reg_K${k}_${p}_${end}
#SBATCH --output=${LOGS_ROOT}/regulation_K${k}_${p}_${end}_%j.out
#SBATCH --error=${LOGS_ROOT}/regulation_K${k}_${p}_${end}_%j.err
#SBATCH --mem=${REG_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${REG_CPUS}
#SBATCH --partition=${REG_PARTITION}
#SBATCH --time=${REG_TIME}
#SBATCH --chdir=${PROJECT_ROOT}
set -euo pipefail
source "${CONDA_SH}"
mkdir -p "${REGULATION_DIR}"
USAGES_PATH="${CNMF_ROOT}/cNMF_all/cNMF_all.usages.k_${k}.dt_${DENSITY_LABEL}.consensus.txt"
[[ -s "\${USAGES_PATH}" ]] || { echo "Missing usages: \${USAGES_PATH}" >&2; exit 1; }
[[ -s "${CD4_METADATA}" ]] || { echo "Missing metadata: ${CD4_METADATA}" >&2; exit 1; }
OUTPUT_PREFIX="${REGULATION_DIR}/K${k}_program"
for program_index in \$(seq ${p} ${end}); do
  conda run --no-capture-output -n "${R_ENV}" Rscript "${REG_SCRIPT}" "${CD4_METADATA}" "\${USAGES_PATH}" "\${program_index}" "\${OUTPUT_PREFIX}\${program_index}_perturb_effects.txt"
done
EOF
        chmod +x "${JOBS_ROOT}/regulation_K${k}_program_${p}_${end}.sh"
    done
done

{
printf '%s\n' '#!/usr/bin/env bash'
printf 'CD4_INPUT_H5AD=%q\n' "${CD4_INPUT_H5AD}"
printf 'STANDARDIZE_SCRIPT=%q\n' "${STANDARDIZE_SCRIPT}"
cat <<'EOF'
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SHARED_DIR="$SCRIPT_DIR/shared"
JOBS_DIR="$SCRIPT_DIR/jobs"
[[ -s "${CD4_INPUT_H5AD}" ]] || { echo "Missing CD4 h5ad: ${CD4_INPUT_H5AD}" >&2; exit 1; }
[[ -s "${STANDARDIZE_SCRIPT}" ]] || { echo "Missing CD4 conversion script: ${STANDARDIZE_SCRIPT}" >&2; exit 1; }
metadata_jobid="$(sbatch --parsable "$SHARED_DIR/cd4_prepare_metadata.sh")"
echo "metadata preparation jobid: $metadata_jobid"
prepare_jobid="$(sbatch --parsable --dependency=afterok:$metadata_jobid "$SHARED_DIR/cd4_cnmf_kselect_prepare.sh")"
echo "cNMF prepare jobid: $prepare_jobid"
factorize_dependency=""
found_factorize=0
shopt -s nullglob
for script in "$SHARED_DIR"/cd4_cnmf_kselect_factorize_*.sh; do
    jobid="$(sbatch --parsable --dependency=afterok:$prepare_jobid "$script")"
    factorize_dependency="${factorize_dependency:+$factorize_dependency:}$jobid"
    found_factorize=1
    echo "factorize jobid: $jobid ($(basename "$script"))"
done
(( found_factorize == 1 )) || { echo "No factorize scripts found" >&2; exit 1; }
postbase_jobid="$(sbatch --parsable --dependency=afterok:$factorize_dependency "$SHARED_DIR/cd4_cnmf_kselect_postbase_base.sh")"
echo "postbase jobid: $postbase_jobid"
for consensus_script in "$JOBS_DIR"/consensus_K*.sh; do
    base="$(basename "$consensus_script")"
    k="$(printf '%s' "$base" | sed -e 's/^consensus_K//' -e 's/\\.sh$//')"
    consensus_jobid="$(sbatch --parsable --dependency=afterok:$postbase_jobid "$consensus_script")"
    echo "consensus K=$k: $consensus_jobid"
    for regulation_script in "$JOBS_DIR"/regulation_K${k}_program_*.sh; do
        regulation_jobid="$(sbatch --parsable --dependency=afterok:$consensus_jobid "$regulation_script")"
        echo "regulation $(basename "$regulation_script"): $regulation_jobid"
    done
done
shopt -u nullglob
EOF
} > "${BATCH_ROOT}/submit_full_chain.sh"
chmod +x "${BATCH_ROOT}/submit_full_chain.sh"

{
printf '%s\n' '#!/usr/bin/env bash'
printf 'CD4_METADATA=%q\n' "${CD4_METADATA}"
printf 'CD4_GENE_MAP=%q\n' "${CD4_GENE_MAP}"
printf 'CD4_REPORT=%q\n' "${CD4_REPORT}"
cat <<'EOF'
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SHARED_DIR="$SCRIPT_DIR/shared"
JOBS_DIR="$SCRIPT_DIR/jobs"
for required_file in "${CD4_METADATA}" "${CD4_GENE_MAP}" "${CD4_REPORT}"; do
    [[ -s "${required_file}" ]] || { echo "Missing CD4 pipeline-ready file: ${required_file}" >&2; exit 1; }
done
prepare_jobid="$(sbatch --parsable "$SHARED_DIR/cd4_cnmf_kselect_prepare.sh")"
echo "cNMF prepare jobid: $prepare_jobid"
factorize_dependency=""
found_factorize=0
shopt -s nullglob
for script in "$SHARED_DIR"/cd4_cnmf_kselect_factorize_*.sh; do
    jobid="$(sbatch --parsable --dependency=afterok:$prepare_jobid "$script")"
    factorize_dependency="${factorize_dependency:+$factorize_dependency:}$jobid"
    found_factorize=1
    echo "factorize jobid: $jobid ($(basename "$script"))"
done
(( found_factorize == 1 )) || { echo "No factorize scripts found" >&2; exit 1; }
postbase_jobid="$(sbatch --parsable --dependency=afterok:$factorize_dependency "$SHARED_DIR/cd4_cnmf_kselect_postbase_base.sh")"
echo "postbase jobid: $postbase_jobid"
for consensus_script in "$JOBS_DIR"/consensus_K*.sh; do
    base="$(basename "$consensus_script")"
    k="$(printf '%s' "$base" | sed -e 's/^consensus_K//' -e 's/\\.sh$//')"
    consensus_jobid="$(sbatch --parsable --dependency=afterok:$postbase_jobid "$consensus_script")"
    echo "consensus K=$k: $consensus_jobid"
    for regulation_script in "$JOBS_DIR"/regulation_K${k}_program_*.sh; do
        regulation_jobid="$(sbatch --parsable --dependency=afterok:$consensus_jobid "$regulation_script")"
        echo "regulation $(basename "$regulation_script"): $regulation_jobid"
    done
done
shopt -u nullglob
EOF
} > "${BATCH_ROOT}/submit_cnmf_regulation_chain.sh"
chmod +x "${BATCH_ROOT}/submit_cnmf_regulation_chain.sh"

cat <<EOF
Generated CD4 batch scripts under: ${BATCH_ROOT}
K values: ${K_VALUES}
Total cNMF workers: ${TOTAL_WORKERS}
Programs per regulation job: ${PROGRAMS_PER_JOB}
Input h5ad used directly by cNMF: ${CD4_INPUT_H5AD}
Pipeline metadata: ${CD4_METADATA}
Gene map: ${CD4_GENE_MAP}
Preparation report: ${CD4_REPORT}
Expected K-selection plot: ${CNMF_ROOT}/cNMF_all/cNMF_all.k_selection.png
Expected regulation directory: ${REGULATION_DIR}
Submit full chain:
  cd ${BATCH_ROOT}
  ./submit_full_chain.sh
Skip conversion when pipeline_ready files exist:
  cd ${BATCH_ROOT}
  ./submit_cnmf_regulation_chain.sh
EOF
