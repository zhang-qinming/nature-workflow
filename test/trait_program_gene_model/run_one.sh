#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOF_ID="${LOF_ID:?LOF_ID not set}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BASE_CONFIG="${BASE_CONFIG:-${SCRIPT_DIR}/config.base.yaml}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-trait_program_gene_model}"
JOB_NAME="${JOB_NAME:-tpgm}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/trait_program_gene_model}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"
PLOT_ENV="${PLOT_ENV:-paper-pipeline-plot}"

TPGM_K="${TPGM_K:-60}"
TPGM_PROGRAM_N="${TPGM_PROGRAM_N:-auto}"
TPGM_REGULATOR_N="${TPGM_REGULATOR_N:-auto}"
TPGM_LOADING_TOP_N="${TPGM_LOADING_TOP_N:-200}"
TPGM_HIT_ABS_GAMMA_THRESHOLD="${TPGM_HIT_ABS_GAMMA_THRESHOLD:-0.1}"
TPGM_REGULATOR_FDR_THRESHOLD="${TPGM_REGULATOR_FDR_THRESHOLD:-0.05}"
TPGM_RANDOM_ITERATIONS="${TPGM_RANDOM_ITERATIONS:-10000}"
TPGM_PERMUTATION_ITERATIONS="${TPGM_PERMUTATION_ITERATIONS:-0}"
TPGM_MAX_GENES_PER_SIDE="${TPGM_MAX_GENES_PER_SIDE:-all}"
TPGM_RENDER_PLOT="${TPGM_RENDER_PLOT:-1}"
TPGM_PROGRAM_FDR_THRESHOLD="${TPGM_PROGRAM_FDR_THRESHOLD:-0.05}"
TPGM_REGULATOR_MAX_N="${TPGM_REGULATOR_MAX_N:-8}"

TPGM_POSTERIOR_DIR="${TPGM_POSTERIOR_DIR:-}"
TPGM_GENE_MAP="${TPGM_GENE_MAP:-}"
TPGM_SHET_PATH="${TPGM_SHET_PATH:-}"
TPGM_SPECTRA_PATH="${TPGM_SPECTRA_PATH:-}"
TPGM_REGULATION_DIR="${TPGM_REGULATION_DIR:-}"

SUCCESS=0
TEMP_META=""

mkdir -p "${STATUS_DIR}" "${FAILURE_DIR}" "${OUTPUT_DIR}/tables" "${OUTPUT_DIR}/plots" "${OUTPUT_DIR}/meta"

cleanup() {
    if [ -n "${TEMP_META}" ] && [ -f "${TEMP_META}" ]; then
        rm -f "${TEMP_META}"
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

TEMP_META="$(mktemp "/tmp/${TASK_NAME}_${LOF_ID}_${SLURM_JOB_ID:-$$}.XXXXXX.tsv")"

python3 - "${LOF_ID}" "${BASE_CONFIG}" "${FILE_ID_MAP}" "${TEMP_META}" "${TPGM_POSTERIOR_DIR}" "${TPGM_GENE_MAP}" "${TPGM_SHET_PATH}" "${TPGM_SPECTRA_PATH}" "${TPGM_REGULATION_DIR}" <<'PYEOF'
import sys
from pathlib import Path

import yaml

lof_id = sys.argv[1]
base_config = Path(sys.argv[2])
file_id_map = Path(sys.argv[3])
out_path = Path(sys.argv[4])
posterior_dir_override = sys.argv[5].strip()
gene_map_override = sys.argv[6].strip()
shet_path_override = sys.argv[7].strip()
spectra_path_override = sys.argv[8].strip()
regulation_dir_override = sys.argv[9].strip()

with open(base_config, encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

shared = config["workflows"]["figures"]["shared_inputs"]
if posterior_dir_override:
    shared["posterior_dir"] = posterior_dir_override
if gene_map_override:
    shared["gene_map"] = gene_map_override
if shet_path_override:
    shared["shet_path"] = shet_path_override
if spectra_path_override:
    shared["spectra_path"] = spectra_path_override
if regulation_dir_override:
    shared["regulation_dir"] = regulation_dir_override
posterior_dir = Path(shared["posterior_dir"])

with open(file_id_map, encoding="utf-8") as handle:
    header = handle.readline().rstrip("\n").split("\t")
    rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

try:
    id2_index = header.index("id2")
    path2_index = header.index("path2")
except ValueError as exc:
    raise SystemExit(f"file_id_map must contain id2 and path2 columns: {file_id_map}") from exc

match = next((row for row in rows if len(row) > max(id2_index, path2_index) and row[id2_index] == lof_id), None)
if match is None:
    raise SystemExit(f"LOF_ID {lof_id} was not found in {file_id_map} column id2")

trait_source = Path(match[path2_index])
stem = trait_source.name
for suffix in (
    ".summary_statistics.csv",
    ".per_gene_estimates.tsv",
    ".tsv.gz",
    ".txt.gz",
    ".csv.gz",
    ".tsv",
    ".csv",
):
    if stem.endswith(suffix):
        stem = stem[: -len(suffix)]
        break
posterior_name = f"{stem}.per_gene_estimates.tsv"
posterior_path = posterior_dir / posterior_name

fields = {
    "lof_id": lof_id,
    "posterior_path": str(posterior_path),
    "trait_source_path": str(trait_source),
    "gene_map": str(shared["gene_map"]),
    "shet_path": str(shared["shet_path"]),
    "spectra_path": str(shared["spectra_path"]),
    "regulation_dir": str(shared["regulation_dir"]),
}
with open(out_path, "w", encoding="utf-8") as handle:
    for key, value in fields.items():
        handle.write(f"{key}\t{value}\n")
PYEOF

lookup_value() {
    awk -F '\t' -v key="$1" '$1 == key { print $2 }' "${TEMP_META}"
}

POSTERIOR_PATH="$(lookup_value posterior_path)"
GENE_MAP="$(lookup_value gene_map)"
SHET_PATH="$(lookup_value shet_path)"
SPECTRA_PATH="$(lookup_value spectra_path)"
REGULATION_DIR="$(lookup_value regulation_dir)"

TABLE_PREFIX="${OUTPUT_DIR}/tables/${LOF_ID}"
PLOT_PREFIX="${OUTPUT_DIR}/plots/${LOF_ID}"

conda run --no-capture-output -n "${PLOT_ENV}" Rscript \
    "${SCRIPT_DIR}/trait_program_gene_model.R" \
    "${REGULATION_DIR}" \
    "${SPECTRA_PATH}" \
    "${GENE_MAP}" \
    "${POSTERIOR_PATH}" \
    "${LOF_ID}" \
    "${TPGM_K}" \
    "${TABLE_PREFIX}" \
    "${PLOT_PREFIX}" \
    "${SHET_PATH}" \
    "${TPGM_PROGRAM_N}" \
    "${TPGM_REGULATOR_N}" \
    "${TPGM_LOADING_TOP_N}" \
    "${TPGM_HIT_ABS_GAMMA_THRESHOLD}" \
    "${TPGM_REGULATOR_FDR_THRESHOLD}" \
    "${TPGM_RANDOM_ITERATIONS}" \
    "${TPGM_PERMUTATION_ITERATIONS}" \
    "${TPGM_MAX_GENES_PER_SIDE}" \
    "${TPGM_RENDER_PLOT}" \
    "${TPGM_PROGRAM_FDR_THRESHOLD}" \
    "${TPGM_REGULATOR_MAX_N}"

printf 'time=%s\njob_id=%s\n' "$(date -Iseconds)" "${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${LOF_ID}.ok"
rm -f "${STATUS_DIR}/${LOF_ID}.failed" "${STATUS_DIR}/${LOF_ID}.running" "${FAILURE_DIR}/${LOF_ID}.failed"
SUCCESS=1
