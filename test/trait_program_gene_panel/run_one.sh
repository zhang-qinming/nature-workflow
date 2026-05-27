#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOF_ID="${LOF_ID:?LOF_ID not set}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BASE_CONFIG="${BASE_CONFIG:-${SCRIPT_DIR}/config.base.yaml}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-trait_program_gene_panel}"
JOB_NAME="${JOB_NAME:-tpgp}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/trait_program_gene_panel}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"

TPGP_K="${TPGP_K:-60}"
TPGP_PROGRAM_N="${TPGP_PROGRAM_N:-5}"
TPGP_REGULATOR_N="${TPGP_REGULATOR_N:-3}"
TPGP_MAX_GENES_PER_SIDE="${TPGP_MAX_GENES_PER_SIDE:-8}"
TPGP_HIT_ABS_GAMMA_THRESHOLD="${TPGP_HIT_ABS_GAMMA_THRESHOLD:-0.1}"
TPGP_LOADING_TOP_N="${TPGP_LOADING_TOP_N:-200}"
TPGP_REGULATOR_FDR_THRESHOLD="${TPGP_REGULATOR_FDR_THRESHOLD:-0.05}"
TPGP_RANDOM_ITERATIONS="${TPGP_RANDOM_ITERATIONS:-10000}"
TPGP_RENDER_PLOT="${TPGP_RENDER_PLOT:-1}"

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
            printf 'lof_id=%s\n' "${LOF_ID}"
            printf 'task_name=%s\n' "${TASK_NAME}"
            printf 'job_name=%s\n' "${JOB_NAME}"
            printf 'job_id=%s\n' "${SLURM_JOB_ID:-local}"
            printf 'host=%s\n' "$(hostname)"
            printf 'exit_code=%s\n' "${exit_code}"
        } > "${STATUS_DIR}/${LOF_ID}.failed"
        cp "${STATUS_DIR}/${LOF_ID}.failed" "${FAILURE_DIR}/${LOF_ID}.failed"
        rm -f "${STATUS_DIR}/${LOF_ID}.running"
    fi
    cleanup
}

trap on_exit EXIT

if [ -f "${CONDA_SH}" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
    conda activate "${CONTROL_ENV}"
fi

TEMP_CONFIG="$(mktemp "/tmp/${TASK_NAME}_${LOF_ID}_${SLURM_JOB_ID:-$$}.XXXXXX.yaml")"

python3 - "${LOF_ID}" "${PROJECT_ROOT}" "${BASE_CONFIG}" "${TEMP_CONFIG}" "${OUTPUT_ROOT}" "${OUTPUT_DIR}" "${FILE_ID_MAP}" "${TPGP_K}" "${TPGP_PROGRAM_N}" "${TPGP_REGULATOR_N}" "${TPGP_MAX_GENES_PER_SIDE}" "${TPGP_HIT_ABS_GAMMA_THRESHOLD}" "${TPGP_LOADING_TOP_N}" "${TPGP_REGULATOR_FDR_THRESHOLD}" "${TPGP_RANDOM_ITERATIONS}" "${TPGP_RENDER_PLOT}" <<'PYEOF'
import sys
from pathlib import Path

import yaml

lof_id = sys.argv[1]
proj_root = Path(sys.argv[2])
base_path = Path(sys.argv[3])
temp_path = Path(sys.argv[4])
output_root = Path(sys.argv[5])
output_dir = Path(sys.argv[6])
file_id_map = Path(sys.argv[7])
k = int(sys.argv[8])
program_n = int(sys.argv[9])
regulator_n = int(sys.argv[10])
max_genes_per_side = int(sys.argv[11])
hit_abs_gamma_threshold = float(sys.argv[12])
loading_top_n = int(sys.argv[13])
regulator_fdr_threshold = float(sys.argv[14])
random_iterations = int(sys.argv[15])
render_plot = sys.argv[16].strip().lower() in {"1", "true", "yes", "on"}

with open(base_path, encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

config["project_root"] = str(proj_root)
config["artifact_root"] = str(output_root)

shared = config.get("workflows", {}).get("figures", {}).get("shared_inputs", {})
for key in (
    "file_id_map",
    "posterior_dir",
    "gene_map",
    "spectra_path",
    "regulation_dir",
    "program_association_dir",
):
    value = shared.get(key)
    if value and isinstance(value, str) and not value.startswith("/"):
        shared[key] = str((base_path.parent / value).resolve())

shared["file_id_map"] = str(file_id_map.resolve())
config["workflows"]["figures"]["shared_inputs"] = shared
config["workflows"]["figures"]["trait_program_gene_panel"] = {
    "outputs": {"output_dir": str(output_dir)},
    "parameters": {
        "lof_ids": [lof_id],
        "k": k,
        "program_n": program_n,
        "regulator_n": regulator_n,
        "max_genes_per_side": max_genes_per_side,
        "hit_abs_gamma_threshold": hit_abs_gamma_threshold,
        "loading_top_n": loading_top_n,
        "regulator_fdr_threshold": regulator_fdr_threshold,
        "random_iterations": random_iterations,
        "render_plot": render_plot,
    },
}

with open(temp_path, "w", encoding="utf-8") as handle:
    yaml.dump(config, handle, default_flow_style=False, sort_keys=False, allow_unicode=True)
PYEOF

paper-pipeline run --config "${TEMP_CONFIG}" figures-trait-program-gene-panel

printf 'time=%s\njob_id=%s\n' "$(date -Iseconds)" "${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${LOF_ID}.ok"
rm -f "${STATUS_DIR}/${LOF_ID}.failed" "${STATUS_DIR}/${LOF_ID}.running" "${FAILURE_DIR}/${LOF_ID}.failed"
SUCCESS=1
