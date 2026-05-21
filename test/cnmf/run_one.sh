#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOF_ID="${LOF_ID:?LOF_ID not set}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BASE_CONFIG="${BASE_CONFIG:-${SCRIPT_DIR}/config.base.yaml}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cnmf}"
JOB_NAME="${JOB_NAME:-cnmf}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cnmf}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"

CNMF_K="${CNMF_K:-60}"
CNMF_PLOT_LABEL_PROGRAMS="${CNMF_PLOT_LABEL_PROGRAMS:-4,16,25,40}"
CNMF_CORREGULATION_PAIRS="${CNMF_CORREGULATION_PAIRS:-}"

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

python3 - "${LOF_ID}" "${PROJECT_ROOT}" "${BASE_CONFIG}" "${TEMP_CONFIG}" "${OUTPUT_ROOT}" "${OUTPUT_DIR}" "${CNMF_K}" "${CNMF_PLOT_LABEL_PROGRAMS}" "${CNMF_CORREGULATION_PAIRS}" <<'PYEOF'
import sys
from pathlib import Path

import yaml

lof_id = sys.argv[1]
proj_root = Path(sys.argv[2])
base_path = Path(sys.argv[3])
temp_path = Path(sys.argv[4])
output_root = Path(sys.argv[5])
output_dir = Path(sys.argv[6])
k = int(sys.argv[7])
plot_label_programs = [int(x.strip()) for x in sys.argv[8].split(",") if x.strip()]
corregulation_pairs_raw = sys.argv[9].strip()

with open(base_path, encoding="utf-8") as handle:
    config = yaml.safe_load(handle)

config["project_root"] = str(proj_root)
config["artifact_root"] = str(output_root)

shared = config.get("workflows", {}).get("figures", {}).get("shared_inputs", {})
for key in ("file_id_map", "program_association_dir", "cnmf_regulation_dir"):
    value = shared.get(key)
    if value and isinstance(value, str) and not value.startswith("/"):
        shared[key] = str((base_path.parent / value).resolve())

config["workflows"]["figures"]["shared_inputs"] = shared

trait_file = None
with open(shared["file_id_map"], encoding="utf-8") as handle:
    next(handle)
    for line in handle:
        parts = line.rstrip("\n\r").split("\t")
        if len(parts) >= 4 and parts[1] == lof_id:
            trait_file = Path(parts[3]).name
            break
if not trait_file:
    raise SystemExit(f"Unable to resolve trait file for {lof_id} from {shared['file_id_map']}")

pairs = []
if corregulation_pairs_raw:
    for raw in corregulation_pairs_raw.split(","):
        raw = raw.strip()
        if not raw:
            continue
        if ":" in raw:
            pair_id, pair_spec = raw.split(":", 1)
        else:
            pair_spec = raw
            pair_id = raw.replace(":", "__")
        if "|" in pair_spec:
            program_a, program_b = [x.strip() for x in pair_spec.split("|", 1)]
        elif ":" in pair_spec:
            program_a, program_b = [x.strip() for x in pair_spec.split(":", 1)]
        else:
            raise SystemExit(f"Invalid CNMF_CORREGULATION_PAIRS entry: {raw}")
        pairs.append({"output_id": pair_id.strip(), "program_a": program_a, "program_b": program_b})

config["workflows"]["figures"]["cnmf"] = {
    "outputs": {"output_dir": str(output_dir)},
    "parameters": {
        "k": k,
        "trait_targets": [{"trait_file": trait_file, "trait_id": lof_id}],
        "plot_label_programs": plot_label_programs,
        "corregulation_pairs": pairs,
    },
}

with open(temp_path, "w", encoding="utf-8") as handle:
    yaml.dump(config, handle, default_flow_style=False, sort_keys=False, allow_unicode=True)
PYEOF

paper-pipeline run --config "${TEMP_CONFIG}" figures-cnmf

printf 'time=%s\njob_id=%s\n' "$(date -Iseconds)" "${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${LOF_ID}.ok"
rm -f "${STATUS_DIR}/${LOF_ID}.failed" "${STATUS_DIR}/${LOF_ID}.running" "${FAILURE_DIR}/${LOF_ID}.failed"
SUCCESS=1
