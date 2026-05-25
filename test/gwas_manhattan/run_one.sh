#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GWAS_ID="${GWAS_ID:?GWAS_ID not set}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BASE_CONFIG="${BASE_CONFIG:-${SCRIPT_DIR}/config.base.yaml}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-gwas_manhattan}"
JOB_NAME="${JOB_NAME:-gman}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/gwas_manhattan}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"

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
            printf 'gwas_id=%s\n' "${GWAS_ID}"
            printf 'task_name=%s\n' "${TASK_NAME}"
            printf 'job_name=%s\n' "${JOB_NAME}"
            printf 'job_id=%s\n' "${SLURM_JOB_ID:-local}"
            printf 'host=%s\n' "$(hostname)"
            printf 'exit_code=%s\n' "${exit_code}"
        } > "${STATUS_DIR}/${GWAS_ID}.failed"
        cp "${STATUS_DIR}/${GWAS_ID}.failed" "${FAILURE_DIR}/${GWAS_ID}.failed"
        rm -f "${STATUS_DIR}/${GWAS_ID}.running"
    fi
    cleanup
}

trap on_exit EXIT

if [ -f "${CONDA_SH}" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
    conda activate "${CONTROL_ENV}"
fi

TEMP_CONFIG="$(mktemp "/tmp/${TASK_NAME}_${GWAS_ID}_${SLURM_JOB_ID:-$$}.XXXXXX.yaml")"

python3 - "${GWAS_ID}" "${PROJECT_ROOT}" "${BASE_CONFIG}" "${TEMP_CONFIG}" "${OUTPUT_ROOT}" "${OUTPUT_DIR}" "${FILE_ID_MAP}" <<'PYEOF'
import sys
import os
from pathlib import Path

import yaml

gwas_id = sys.argv[1]
proj_root = Path(sys.argv[2])
base_path = Path(sys.argv[3])
temp_path = Path(sys.argv[4])
output_root = Path(sys.argv[5])
output_dir = Path(sys.argv[6])
file_id_map = Path(sys.argv[7])

with open(base_path, encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

config["project_root"] = str(proj_root)
config["artifact_root"] = str(output_root)

shared = config.get("workflows", {}).get("figures", {}).get("shared_inputs", {})
for key in ("file_id_map", "gene_map", "gene_annotation", "geneset_dir", "spectra_path"):
    value = shared.get(key)
    if value and isinstance(value, str) and not value.startswith("/"):
        shared[key] = str((base_path.parent / value).resolve())

shared["file_id_map"] = str(file_id_map.resolve())
config["workflows"]["figures"]["shared_inputs"] = shared

shared_parameters = config.get("workflows", {}).get("figures", {}).get("shared_parameters", {})
if not isinstance(shared_parameters, dict):
    shared_parameters = {}

def env_or_default(name, default, caster):
    raw = os.environ.get(name)
    if raw is None or raw == "":
        return default
    return caster(raw)

sampling_trigger_rows = env_or_default(
    "GWAS_SAMPLING_TRIGGER_ROWS",
    int(shared_parameters.get("sampling_trigger_rows", 100000)),
    int,
)
sampling_base_points = env_or_default(
    "GWAS_SAMPLING_BASE_POINTS",
    int(shared_parameters.get("sampling_base_points", 50000)),
    int,
)
sampling_fraction = env_or_default(
    "GWAS_SAMPLING_FRACTION",
    float(shared_parameters.get("sampling_fraction", 0.01)),
    float,
)
sampling_max_points = env_or_default(
    "GWAS_SAMPLING_MAX_POINTS",
    int(shared_parameters.get("sampling_max_points", 300000)),
    int,
)
sampling_seed = env_or_default(
    "GWAS_SAMPLING_SEED",
    int(shared_parameters.get("sampling_seed", 1)),
    int,
)

config["workflows"]["figures"]["gwas_manhattan"] = {
    "outputs": {"output_dir": str(output_dir)},
    "parameters": {
        "gwas_ids": [gwas_id],
        "flank_bp": 50000,
        "label_p_threshold": 1e-30,
        "genomewide_threshold": 5e-8,
        "highlight_genesets": [
            "HALLMARK_HEME_METABOLISM",
            "Hematopoiesisgenes",
            "mitotic_cell_cycle",
            "positive_macromolecule_synthesis",
        ],
    },
}

parameters = config["workflows"]["figures"]["gwas_manhattan"]["parameters"]
parameters["sampling_trigger_rows"] = sampling_trigger_rows
parameters["sampling_base_points"] = sampling_base_points
parameters["sampling_fraction"] = sampling_fraction
parameters["sampling_max_points"] = sampling_max_points
parameters["sampling_seed"] = sampling_seed

with open(temp_path, "w", encoding="utf-8") as handle:
    yaml.dump(config, handle, default_flow_style=False, sort_keys=False, allow_unicode=True)
PYEOF

paper-pipeline run --config "${TEMP_CONFIG}" figures-gwas-manhattan

printf 'time=%s\njob_id=%s\n' "$(date -Iseconds)" "${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${GWAS_ID}.ok"
rm -f "${STATUS_DIR}/${GWAS_ID}.failed" "${STATUS_DIR}/${GWAS_ID}.running" "${FAILURE_DIR}/${GWAS_ID}.failed"
SUCCESS=1
