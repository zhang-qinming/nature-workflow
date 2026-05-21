#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cnmf}"
JOB_NAME="${JOB_NAME:-cnmf}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cnmf}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"

STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
LOGS_DIR="${LOGS_DIR:-${RUN_ROOT}/logs}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"
GENERATED_DIR="${GENERATED_DIR:-${BATCH_ROOT}/jobs}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"

WORKER_SCRIPT="${WORKER_SCRIPT:-${SCRIPT_DIR}/run_one.sh}"
SBATCH_MEM="${SBATCH_MEM:-8G}"
SBATCH_CPUS="${SBATCH_CPUS:-2}"
SBATCH_PARTITION="${SBATCH_PARTITION:-cu,fat,batch01,privority}"
SBATCH_TIME="${SBATCH_TIME:-12:00:00}"

CNMF_K="${CNMF_K:-60}"
CNMF_PLOT_LABEL_PROGRAMS="${CNMF_PLOT_LABEL_PROGRAMS:-4,16,25,40}"
CNMF_CORREGULATION_PAIRS="${CNMF_CORREGULATION_PAIRS:-}"

mkdir -p "${BATCH_ROOT}" "${GENERATED_DIR}" "${STATUS_DIR}" "${LOGS_DIR}" "${FAILURE_DIR}"

if [ ! -f "${FILE_ID_MAP}" ]; then
    echo "FILE_ID_MAP not found: ${FILE_ID_MAP}" >&2
    exit 1
fi

if [ ! -f "${WORKER_SCRIPT}" ]; then
    echo "WORKER_SCRIPT not found: ${WORKER_SCRIPT}" >&2
    exit 1
fi

mapfile -t ALL_IDS < <(awk -F '\t' 'NR > 1 { sub(/\r$/, "", $2); print $2 }' "${FILE_ID_MAP}")

if [ "${#ALL_IDS[@]}" -eq 0 ]; then
    echo "No LoF IDs found in ${FILE_ID_MAP}" >&2
    exit 1
fi

rm -f "${GENERATED_DIR}"/*.sbatch

{
    printf 'lof_id\tscript_path\n'
    for lof_id in "${ALL_IDS[@]}"; do
        script_path="${GENERATED_DIR}/${lof_id}.sbatch"
        cat > "${script_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOGS_DIR}/${lof_id}_%j.out
#SBATCH --error=${LOGS_DIR}/${lof_id}_%j.err
#SBATCH --mem=${SBATCH_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SBATCH_CPUS}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

export LOF_ID=${lof_id}
export PROJECT_ROOT=${PROJECT_ROOT}
export FILE_ID_MAP=${FILE_ID_MAP}
export TASK_NAME=${TASK_NAME}
export JOB_NAME=${JOB_NAME}
export OUTPUT_ROOT=${OUTPUT_ROOT}
export OUTPUT_DIR=${OUTPUT_DIR}
export RUN_ROOT=${RUN_ROOT}
export STATUS_DIR=${STATUS_DIR}
export LOGS_DIR=${LOGS_DIR}
export FAILURE_DIR=${FAILURE_DIR}
export CNMF_K=${CNMF_K}
export CNMF_PLOT_LABEL_PROGRAMS=${CNMF_PLOT_LABEL_PROGRAMS}
export CNMF_CORREGULATION_PAIRS=${CNMF_CORREGULATION_PAIRS}

bash ${WORKER_SCRIPT}
EOF
        chmod +x "${script_path}"
        printf '%s\t%s\n' "${lof_id}" "${script_path}"
    done
} > "${MANIFEST_PATH}"

echo "Generated ${#ALL_IDS[@]} sbatch scripts under: ${GENERATED_DIR}"
echo "Manifest: ${MANIFEST_PATH}"
echo "Job name: ${JOB_NAME}"
echo "Next: bash test/cnmf/2_submit.sh"
