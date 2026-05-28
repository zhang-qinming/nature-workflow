#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_ID="${OUTPUT_ID:?OUTPUT_ID not set}"
LOF_IDS="${LOF_IDS:?LOF_IDS not set}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BASE_CONFIG="${BASE_CONFIG:-${SCRIPT_DIR}/config.base.yaml}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cross_trait_heatmap}"
JOB_NAME="${JOB_NAME:-ctheat}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cross_trait_heatmap}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"

CROSS_TRAIT_METHOD="${CROSS_TRAIT_METHOD:-pearson}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"

SUCCESS=0
TEMP_CONFIG=""

mkdir -p "${STATUS_DIR}" "${FAILURE_DIR}"

cleanup() {
    if [ -n "${TEMP_CONFIG}" ] && [ -f "${TEMP_CONFIG}" ]; then
        rm -f "${TEMP_CONFIG}"
    fi
}

on_exit() {
    exit_code=$?
    if [ "${SUCCESS}" -ne 1 ] && [ "${exit_code}" -ne 0 ]; then
        {
            printf 'time=%s\n' "$(date -Iseconds)"
            printf 'output_id=%s\n' "${OUTPUT_ID}"
            printf 'lof_ids=%s\n' "${LOF_IDS}"
            printf 'task_name=%s\n' "${TASK_NAME}"
            printf 'job_name=%s\n' "${JOB_NAME}"
            printf 'job_id=%s\n' "${SLURM_JOB_ID:-local}"
            printf 'host=%s\n' "$(hostname)"
            printf 'exit_code=%s\n' "${exit_code}"
        } > "${STATUS_DIR}/${OUTPUT_ID}.failed"
        cp "${STATUS_DIR}/${OUTPUT_ID}.failed" "${FAILURE_DIR}/${OUTPUT_ID}.failed"
        rm -f "${STATUS_DIR}/${OUTPUT_ID}.running"
    fi
    cleanup
}

trap on_exit EXIT

if [ -f "${CONDA_SH}" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
    conda activate "${CONTROL_ENV}"
fi

TEMP_CONFIG="$(mktemp "/tmp/${TASK_NAME}_${OUTPUT_ID}_${SLURM_JOB_ID:-$$}.XXXXXX.yaml")"

python3 - "${OUTPUT_ID}" "${LOF_IDS}" "${PROJECT_ROOT}" "${BASE_CONFIG}" "${TEMP_CONFIG}" "${OUTPUT_ROOT}" "${OUTPUT_DIR}" "${FILE_ID_MAP}" "${CROSS_TRAIT_METHOD}" <<'PYEOF'
import sys
from pathlib import Path

import yaml

output_id = sys.argv[1]
lof_ids = [item.strip() for item in sys.argv[2].split(",") if item.strip()]
proj_root = Path(sys.argv[3])
base_path = Path(sys.argv[4])
temp_path = Path(sys.argv[5])
output_root = Path(sys.argv[6])
output_dir = Path(sys.argv[7])
file_id_map = Path(sys.argv[8])
method = sys.argv[9]

if len(lof_ids) < 2:
    raise SystemExit("cross_trait_heatmap requires at least two LoF IDs")

with open(base_path, encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

config["project_root"] = str(proj_root)
config["artifact_root"] = str(output_root)

shared = config.get("workflows", {}).get("figures", {}).get("shared_inputs", {})
for key in ("file_id_map", "posterior_dir", "gene_map"):
    value = shared.get(key)
    if value and isinstance(value, str) and not value.startswith("/"):
        shared[key] = str((base_path.parent / value).resolve())

shared["file_id_map"] = str(file_id_map.resolve())
config["workflows"]["figures"]["shared_inputs"] = shared
config["workflows"]["figures"]["cross_trait_heatmap"] = {
    "outputs": {"output_dir": str(output_dir)},
    "parameters": {
        "lof_ids": lof_ids,
        "output_id": output_id,
        "method": method,
    },
}

with open(temp_path, "w", encoding="utf-8") as handle:
    yaml.dump(config, handle, default_flow_style=False, sort_keys=False, allow_unicode=True)
PYEOF

paper-pipeline run --config "${TEMP_CONFIG}" figures-cross-trait-heatmap

printf 'time=%s\njob_id=%s\n' "$(date -Iseconds)" "${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${OUTPUT_ID}.ok"
rm -f "${STATUS_DIR}/${OUTPUT_ID}.failed" "${STATUS_DIR}/${OUTPUT_ID}.running" "${FAILURE_DIR}/${OUTPUT_ID}.failed"
SUCCESS=1
