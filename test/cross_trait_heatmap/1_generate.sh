#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cross_trait_heatmap}"
JOB_NAME="${JOB_NAME:-ctheat}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cross_trait_heatmap}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"

STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
LOGS_DIR="${LOGS_DIR:-${RUN_ROOT}/logs}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"
GENERATED_DIR="${GENERATED_DIR:-${BATCH_ROOT}/jobs}"
MANIFEST_PATH="${MANIFEST_PATH:-${BATCH_ROOT}/manifest.tsv}"

WORKER_SCRIPT="${WORKER_SCRIPT:-${SCRIPT_DIR}/run_one.sh}"
SBATCH_MEM="${SBATCH_MEM:-32G}"
SBATCH_CPUS="${SBATCH_CPUS:-4}"
SBATCH_PARTITION="${SBATCH_PARTITION:-cu,fat,batch01,privority}"
SBATCH_TIME="${SBATCH_TIME:-24:00:00}"

CROSS_TRAIT_IDS="${CROSS_TRAIT_IDS:-}"
CROSS_TRAIT_METHOD="${CROSS_TRAIT_METHOD:-pearson}"
CROSS_TRAIT_OUTPUT_ID="${CROSS_TRAIT_OUTPUT_ID:-}"
CROSS_TRAIT_RENDER_PLOT="${CROSS_TRAIT_RENDER_PLOT:-0}"

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

IDS=()
if [ -n "${CROSS_TRAIT_IDS}" ]; then
    IFS=',' read -r -a RAW_IDS <<< "${CROSS_TRAIT_IDS}"
    for raw_id in "${RAW_IDS[@]}"; do
        raw_id="$(trim "${raw_id}")"
        [ -z "${raw_id}" ] && continue
        IDS+=("${raw_id}")
    done
else
    mapfile -t IDS < <(awk -F '\t' 'NR > 1 { sub(/\r$/, "", $2); print $2 }' "${FILE_ID_MAP}")
fi

if [ "${#IDS[@]}" -lt 2 ]; then
    echo "cross_trait_heatmap requires at least two LoF IDs." >&2
    echo "Set CROSS_TRAIT_IDS=id1,id2,... or provide at least two rows in FILE_ID_MAP." >&2
    exit 1
fi

if [ "${CROSS_TRAIT_RENDER_PLOT}" != "0" ] && [ "${CROSS_TRAIT_RENDER_PLOT,,}" != "false" ] && [ "${#IDS[@]}" -gt 200 ]; then
    echo "Refusing to render a static heatmap for ${#IDS[@]} traits." >&2
    echo "Use CROSS_TRAIT_RENDER_PLOT=0 for frontend TSV output, or set a smaller CROSS_TRAIT_IDS subset." >&2
    exit 1
fi

if [ -z "${CROSS_TRAIT_OUTPUT_ID}" ]; then
    first_id="${IDS[0]}"
    last_id="${IDS[$((${#IDS[@]} - 1))]}"
    CROSS_TRAIT_OUTPUT_ID="$(safe_name "cross_trait_heatmap_n${#IDS[@]}__${first_id}__${last_id}")"
else
    CROSS_TRAIT_OUTPUT_ID="$(safe_name "${CROSS_TRAIT_OUTPUT_ID}")"
fi

LOF_IDS_CSV="$(IFS=','; printf '%s' "${IDS[*]}")"
script_path="${GENERATED_DIR}/${CROSS_TRAIT_OUTPUT_ID}.sbatch"

rm -f "${GENERATED_DIR}"/*.sbatch

cat > "${script_path}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${LOGS_DIR}/${CROSS_TRAIT_OUTPUT_ID}_%j.out
#SBATCH --error=${LOGS_DIR}/${CROSS_TRAIT_OUTPUT_ID}_%j.err
#SBATCH --mem=${SBATCH_MEM}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SBATCH_CPUS}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --chdir=${PROJECT_ROOT}

set -euo pipefail

export OUTPUT_ID=${CROSS_TRAIT_OUTPUT_ID}
export LOF_IDS=${LOF_IDS_CSV}
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
export CROSS_TRAIT_METHOD=${CROSS_TRAIT_METHOD}
export CROSS_TRAIT_RENDER_PLOT=${CROSS_TRAIT_RENDER_PLOT}

bash ${WORKER_SCRIPT}
EOF
chmod +x "${script_path}"

{
    printf 'output_id\tlof_ids\tmethod\tscript_path\n'
    printf '%s\t%s\t%s\t%s\n' "${CROSS_TRAIT_OUTPUT_ID}" "${LOF_IDS_CSV}" "${CROSS_TRAIT_METHOD}" "${script_path}"
} > "${MANIFEST_PATH}"

echo "Generated cross-trait heatmap sbatch script: ${script_path}"
echo "Manifest: ${MANIFEST_PATH}"
echo "Output ID: ${CROSS_TRAIT_OUTPUT_ID}"
echo "LoF IDs: ${#IDS[@]}"
echo "Render plot: ${CROSS_TRAIT_RENDER_PLOT}"
echo "Note: cross_trait_heatmap is group-level; the full trait set generates one matrix TSV job for frontend rendering."
echo "Next: bash test/cross_trait_heatmap/2_submit.sh"
