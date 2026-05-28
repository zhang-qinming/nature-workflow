#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cross_trait_scatter}"
JOB_NAME="${JOB_NAME:-ctscat}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cross_trait}"
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

# Prefer explicit pairs: CROSS_TRAIT_PAIRS=id_x:id_y,id_a:id_b
# CROSS_TRAIT_IDS can be used for small ad hoc all-vs-all batches.
CROSS_TRAIT_PAIRS="${CROSS_TRAIT_PAIRS:-}"
CROSS_TRAIT_IDS="${CROSS_TRAIT_IDS:-}"
CROSS_TRAIT_MAX_PAIRS="${CROSS_TRAIT_MAX_PAIRS:-100}"
CROSS_TRAIT_TOP_N_LABELS="${CROSS_TRAIT_TOP_N_LABELS:-12}"
CROSS_TRAIT_HIGHLIGHT_GENES="${CROSS_TRAIT_HIGHLIGHT_GENES:-}"

mkdir -p "${BATCH_ROOT}" "${GENERATED_DIR}" "${STATUS_DIR}" "${LOGS_DIR}" "${FAILURE_DIR}"

if [ ! -f "${FILE_ID_MAP}" ]; then
    echo "FILE_ID_MAP not found: ${FILE_ID_MAP}" >&2
    exit 1
fi

if [ ! -f "${WORKER_SCRIPT}" ]; then
    echo "WORKER_SCRIPT not found: ${WORKER_SCRIPT}" >&2
    exit 1
fi

trim() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "$value"
}

safe_name() {
    local raw="$1"
    raw="${raw//[^A-Za-z0-9_.-]/_}"
    printf '%s' "$raw"
}

pairs=()
if [ -n "${CROSS_TRAIT_PAIRS}" ]; then
    IFS=',' read -r -a RAW_PAIRS <<< "${CROSS_TRAIT_PAIRS}"
    for raw_pair in "${RAW_PAIRS[@]}"; do
        raw_pair="$(trim "${raw_pair}")"
        [ -z "${raw_pair}" ] && continue
        IFS=':' read -r id_x id_y <<< "${raw_pair}"
        id_x="$(trim "${id_x:-}")"
        id_y="$(trim "${id_y:-}")"
        if [ -z "${id_x}" ] || [ -z "${id_y}" ]; then
            echo "Invalid CROSS_TRAIT_PAIRS entry: ${raw_pair}" >&2
            exit 1
        fi
        pairs+=("${id_x}:${id_y}")
    done
elif [ -n "${CROSS_TRAIT_IDS}" ]; then
    IFS=',' read -r -a RAW_IDS <<< "${CROSS_TRAIT_IDS}"
    ids=()
    for raw_id in "${RAW_IDS[@]}"; do
        raw_id="$(trim "${raw_id}")"
        [ -z "${raw_id}" ] && continue
        ids+=("${raw_id}")
    done
    if [ "${#ids[@]}" -lt 2 ]; then
        echo "CROSS_TRAIT_IDS requires at least two IDs." >&2
        exit 1
    fi
    for ((i = 0; i < ${#ids[@]}; i++)); do
        for ((j = i + 1; j < ${#ids[@]}; j++)); do
            pairs+=("${ids[$i]}:${ids[$j]}")
        done
    done
else
    echo "Set CROSS_TRAIT_PAIRS=id_x:id_y,id_a:id_b to generate scatter jobs." >&2
    echo "For small all-vs-all batches, set CROSS_TRAIT_IDS=id1,id2,id3." >&2
    exit 1
fi

if [ "${#pairs[@]}" -eq 0 ]; then
    echo "No cross-trait pairs were resolved." >&2
    exit 1
fi

if [ "${#pairs[@]}" -gt "${CROSS_TRAIT_MAX_PAIRS}" ]; then
    echo "Refusing to generate ${#pairs[@]} scatter jobs; CROSS_TRAIT_MAX_PAIRS=${CROSS_TRAIT_MAX_PAIRS}." >&2
    echo "Increase CROSS_TRAIT_MAX_PAIRS if this is intentional." >&2
    exit 1
fi

rm -f "${GENERATED_DIR}"/*.sbatch

{
    printf 'output_id\tid_x\tid_y\tscript_path\n'
    for pair in "${pairs[@]}"; do
        IFS=':' read -r id_x id_y <<< "${pair}"
        output_id="$(safe_name "${id_x}__${id_y}")"
        script_path="${GENERATED_DIR}/${output_id}.sbatch"
        cat > "${script_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOGS_DIR}/${output_id}_%j.out
#SBATCH --error=${LOGS_DIR}/${output_id}_%j.err
#SBATCH --mem=${SBATCH_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SBATCH_CPUS}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

export OUTPUT_ID=${output_id}
export ID_X=${id_x}
export ID_Y=${id_y}
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
export CROSS_TRAIT_TOP_N_LABELS=${CROSS_TRAIT_TOP_N_LABELS}
export CROSS_TRAIT_HIGHLIGHT_GENES=${CROSS_TRAIT_HIGHLIGHT_GENES}

bash ${WORKER_SCRIPT}
EOF
        chmod +x "${script_path}"
        printf '%s\t%s\t%s\t%s\n' "${output_id}" "${id_x}" "${id_y}" "${script_path}"
    done
} > "${MANIFEST_PATH}"

echo "Generated ${#pairs[@]} cross-trait scatter sbatch scripts under: ${GENERATED_DIR}"
echo "Manifest: ${MANIFEST_PATH}"
echo "Job name: ${JOB_NAME}"
echo "Next: bash test/cross_trait_scatter/2_submit.sh"
