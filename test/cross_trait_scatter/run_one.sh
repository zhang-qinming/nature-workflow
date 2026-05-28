#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_ID="${OUTPUT_ID:?OUTPUT_ID not set}"
ID_X="${ID_X:?ID_X not set}"
ID_Y="${ID_Y:?ID_Y not set}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BASE_CONFIG="${BASE_CONFIG:-${SCRIPT_DIR}/config.base.yaml}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cross_trait_scatter}"
JOB_NAME="${JOB_NAME:-ctscat}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cross_trait}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"

CROSS_TRAIT_TOP_N_LABELS="${CROSS_TRAIT_TOP_N_LABELS:-12}"
CROSS_TRAIT_HIGHLIGHT_GENES="${CROSS_TRAIT_HIGHLIGHT_GENES:-}"

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
            printf 'id_x=%s\n' "${ID_X}"
            printf 'id_y=%s\n' "${ID_Y}"
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

python3 - "${OUTPUT_ID}" "${ID_X}" "${ID_Y}" "${PROJECT_ROOT}" "${BASE_CONFIG}" "${TEMP_CONFIG}" "${OUTPUT_ROOT}" "${OUTPUT_DIR}" "${FILE_ID_MAP}" "${CROSS_TRAIT_TOP_N_LABELS}" "${CROSS_TRAIT_HIGHLIGHT_GENES}" <<'PYEOF'
import sys
from pathlib import Path

import yaml

output_id = sys.argv[1]
id_x = sys.argv[2]
id_y = sys.argv[3]
proj_root = Path(sys.argv[4])
base_path = Path(sys.argv[5])
temp_path = Path(sys.argv[6])
output_root = Path(sys.argv[7])
output_dir = Path(sys.argv[8])
file_id_map = Path(sys.argv[9])
top_n_labels = int(sys.argv[10])
highlight_genes = [item.strip() for item in sys.argv[11].split(",") if item.strip()]

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
config["workflows"]["figures"]["cross_trait"] = {
    "outputs": {"output_dir": str(output_dir)},
    "parameters": {
        "top_n_labels": top_n_labels,
        "highlight_genes": highlight_genes,
        "lof_pairs": [
            {
                "id_x": id_x,
                "id_y": id_y,
                "output_id": output_id,
            }
        ],
    },
}

with open(temp_path, "w", encoding="utf-8") as handle:
    yaml.dump(config, handle, default_flow_style=False, sort_keys=False, allow_unicode=True)
PYEOF

paper-pipeline run --config "${TEMP_CONFIG}" figures-cross-trait

printf 'time=%s\njob_id=%s\n' "$(date -Iseconds)" "${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${OUTPUT_ID}.ok"
rm -f "${STATUS_DIR}/${OUTPUT_ID}.failed" "${STATUS_DIR}/${OUTPUT_ID}.running" "${FAILURE_DIR}/${OUTPUT_ID}.failed"
SUCCESS=1
