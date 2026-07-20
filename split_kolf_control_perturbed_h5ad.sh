#!/usr/bin/env bash

set -euo pipefail

SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
PROJECT_ROOT="${PROJECT_ROOT:-/gpfs/chencao/qinminzhang/Nature/mine/code}"
INPUT_H5AD="${INPUT_H5AD:-${PROJECT_ROOT}/data/KOLF/KOLF_Pan_Genome_QC_Filtered.h5ad}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_ROOT}/data/KOLF/control_perturbed_target_genes}"
LOG_DIR="${LOG_DIR:-${PROJECT_ROOT}/scripts/kolf_split_logs}"
CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-paper-pipeline-perturbseq}"

PARTITION="${PARTITION:-fat}"
MEMORY="${MEMORY:-900G}"
CPUS="${CPUS:-16}"
TIME_LIMIT="${TIME_LIMIT:-7-00:00:00}"
MIN_FREE_GB="${MIN_FREE_GB:-100}"
OVERWRITE="${OVERWRITE:-0}"

usage() {
    cat <<'EOF'
Usage:
  bash split_kolf_control_perturbed_h5ad.sh          Submit itself with sbatch.
  bash split_kolf_control_perturbed_h5ad.sh --run    Run in the current allocation.

Environment overrides:
  PROJECT_ROOT INPUT_H5AD OUTPUT_DIR CONDA_SH CONDA_ENV
  PARTITION MEMORY CPUS TIME_LIMIT MIN_FREE_GB OVERWRITE

Outputs:
  KOLF.control.target_genes.h5ad
  KOLF.perturbed.target_genes.h5ad
EOF
}

case "${1:-}" in
    -h|--help)
        usage
        exit 0
        ;;
    ""|--run)
        ;;
    *)
        echo "Unknown argument: $1" >&2
        usage >&2
        exit 2
        ;;
esac

if [[ -z "${SLURM_JOB_ID:-}" && "${1:-}" != "--run" ]]; then
    command -v sbatch >/dev/null 2>&1 || {
        echo "sbatch is not available; use --run inside an allocated node." >&2
        exit 1
    }
    mkdir -p "${LOG_DIR}"
    job_id="$(sbatch \
        --parsable \
        --job-name=kolf_split_groups \
        --output="${LOG_DIR}/kolf_split_groups_%j.out" \
        --error="${LOG_DIR}/kolf_split_groups_%j.err" \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="${CPUS}" \
        --mem="${MEMORY}" \
        --partition="${PARTITION}" \
        --time="${TIME_LIMIT}" \
        --chdir="${PROJECT_ROOT}" \
        --export=ALL \
        "${SCRIPT_PATH}" --run | cut -d';' -f1)"
    echo "Submitted Slurm job: ${job_id}"
    echo "stdout: ${LOG_DIR}/kolf_split_groups_${job_id}.out"
    echo "stderr: ${LOG_DIR}/kolf_split_groups_${job_id}.err"
    exit 0
fi

[[ -s "${INPUT_H5AD}" ]] || {
    echo "Input H5AD is missing or empty: ${INPUT_H5AD}" >&2
    exit 1
}
[[ -r "${CONDA_SH}" ]] || {
    echo "Conda initialization file not found: ${CONDA_SH}" >&2
    exit 1
}
[[ "${MIN_FREE_GB}" =~ ^[0-9]+$ ]] || {
    echo "MIN_FREE_GB must be a non-negative integer" >&2
    exit 2
}
[[ "${OVERWRITE}" == "0" || "${OVERWRITE}" == "1" ]] || {
    echo "OVERWRITE must be 0 or 1" >&2
    exit 2
}

mkdir -p "${OUTPUT_DIR}"
exec 9>"${OUTPUT_DIR}/.split_kolf.lock"
if ! flock -n 9; then
    echo "Another KOLF split process is already using ${OUTPUT_DIR}" >&2
    exit 1
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${CPUS}}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${CPUS}}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${CPUS}}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${CPUS}}"
export HDF5_USE_FILE_LOCKING=FALSE

python -u - \
    "${INPUT_H5AD}" \
    "${OUTPUT_DIR}" \
    "${MIN_FREE_GB}" \
    "${OVERWRITE}" <<'PY'
from __future__ import annotations

import gc
import os
import re
import shutil
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


INPUT = Path(sys.argv[1]).resolve()
OUTPUT_DIR = Path(sys.argv[2]).resolve()
MIN_FREE_GB = int(sys.argv[3])
OVERWRITE = bool(int(sys.argv[4]))

OUTPUTS = {
    "control": OUTPUT_DIR / "KOLF.control.target_genes.h5ad",
    "perturbed": OUTPUT_DIR / "KOLF.perturbed.target_genes.h5ad",
}


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def normalize_perturbed(values: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    normalized = values.astype("string").str.strip().str.lower()
    perturbed = normalized.isin({"true", "1"}).to_numpy(dtype=bool)
    control = normalized.isin({"false", "0"}).to_numpy(dtype=bool)
    unknown = ~(perturbed | control)
    if unknown.any():
        examples = sorted(normalized[unknown].fillna("<NA>").unique().tolist())[:20]
        raise RuntimeError(
            "obs['perturbed'] contains values other than True/False: "
            + ", ".join(map(str, examples))
        )
    return control, perturbed


def strip_ensembl_version(value: object) -> str:
    return str(value).strip().split(".", 1)[0]


def split_identifiers(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [
        strip_ensembl_version(item)
        for item in re.split(r"[;,|]", str(value))
        if item.strip()
    ]


def select_target_genes(
    obs: pd.DataFrame, var: pd.DataFrame, perturbed_mask: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    required = {"gene_target", "gene_target_ensembl_id"}
    missing = sorted(required - set(obs.columns))
    if missing:
        raise RuntimeError("Missing required obs columns: " + ", ".join(missing))

    pairs = obs.loc[
        perturbed_mask, ["gene_target", "gene_target_ensembl_id"]
    ].copy()
    pairs["gene_target"] = pairs["gene_target"].astype("string").str.strip()
    pairs = pairs[pairs["gene_target"].notna() & pairs["gene_target"].ne("")]
    pairs = pairs[
        ~pairs["gene_target"].str.lower().isin(
            {"ntc", "non-targeting", "non_targeting", "control"}
        )
    ]
    if pairs.empty:
        raise RuntimeError("No target genes were found among perturbed cells")

    target_counts = pairs["gene_target"].value_counts().sort_index()
    var_names = pd.Index(var.index.astype(str))
    if not var_names.is_unique:
        raise RuntimeError("adata.var_names is not unique")
    name_to_position = {name: position for position, name in enumerate(var_names)}

    id_to_positions: dict[str, list[int]] = {}
    if "gene_ids" in var.columns:
        for position, value in enumerate(var["gene_ids"]):
            if pd.notna(value):
                gene_id = strip_ensembl_version(value)
                id_to_positions.setdefault(gene_id, []).append(position)

    ids_by_target: dict[str, set[str]] = {}
    for target, identifier in pairs.itertuples(index=False, name=None):
        ids_by_target.setdefault(str(target), set()).update(split_identifiers(identifier))

    position_to_targets: dict[int, list[str]] = {}
    unmatched: list[str] = []
    for target in target_counts.index.astype(str):
        position: int | None = name_to_position.get(target)
        if position is None:
            candidate_ids = set(ids_by_target.get(target, set()))
            if target.startswith("ENSG"):
                candidate_ids.add(strip_ensembl_version(target))
            candidates = {
                position
                for gene_id in candidate_ids
                for position in id_to_positions.get(gene_id, [])
            }
            if len(candidates) == 1:
                position = next(iter(candidates))
        if position is None:
            unmatched.append(target)
        else:
            position_to_targets.setdefault(position, []).append(target)

    match_rate = 1 - len(unmatched) / len(target_counts)
    if match_rate < 0.99:
        raise RuntimeError(
            f"Only {match_rate:.2%} of target genes matched the expression matrix"
        )

    positions = np.asarray(sorted(position_to_targets), dtype=np.int64)
    selected_var = var.iloc[positions].copy()
    selected_var["crispri_target"] = True
    selected_var["crispri_target_symbols"] = [
        ";".join(sorted(position_to_targets[int(position)]))
        for position in positions
    ]
    return positions, selected_var, unmatched


def validate_output(
    path: Path, expected_obs_names: pd.Index, expected_var_names: pd.Index
) -> None:
    with h5py.File(path, "r") as handle:
        expected_shape = (len(expected_obs_names), len(expected_var_names))
        for key in ("X", "layers/counts"):
            if key not in handle:
                raise RuntimeError(f"{path} does not contain {key}")
            matrix = handle[key]
            shape = tuple(int(value) for value in matrix.attrs.get("shape", ()))
            if shape != expected_shape:
                raise RuntimeError(
                    f"{path}:{key} has shape {shape}; expected {expected_shape}"
                )
            if matrix.attrs.get("encoding-type") not in {"csr_matrix", "csc_matrix"}:
                raise RuntimeError(f"{path}:{key} is not a sparse matrix")

        output_obs = ad.io.read_elem(handle["obs"])
        output_var = ad.io.read_elem(handle["var"])
        if not np.array_equal(
            output_obs.index.astype(str).to_numpy(),
            expected_obs_names.astype(str).to_numpy(),
        ):
            raise RuntimeError(f"{path} cell order does not match the source")
        if not np.array_equal(
            output_var.index.astype(str).to_numpy(),
            expected_var_names.astype(str).to_numpy(),
        ):
            raise RuntimeError(f"{path} gene order does not match the selection")


def output_is_valid(
    path: Path, expected_obs_names: pd.Index, expected_var_names: pd.Index
) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        validate_output(path, expected_obs_names, expected_var_names)
    except Exception as exc:
        if not OVERWRITE:
            raise RuntimeError(
                f"Existing output is invalid: {path}\n"
                "Set OVERWRITE=1 only after checking that this path may be replaced."
            ) from exc
        return False
    return True


free_gb = shutil.disk_usage(OUTPUT_DIR).free / 1024**3
log(f"Free disk space: {free_gb:.1f} GiB")
if free_gb < MIN_FREE_GB:
    raise RuntimeError(
        f"Only {free_gb:.1f} GiB is free; at least {MIN_FREE_GB} GiB is required"
    )

log(f"Loading the complete H5AD into memory: {INPUT}")
adata = ad.read_h5ad(INPUT)
if "perturbed" not in adata.obs.columns:
    raise RuntimeError("Input H5AD is missing obs['perturbed']")
if "counts" not in adata.layers:
    raise RuntimeError("Input H5AD is missing layers['counts']")
if not sparse.issparse(adata.X) or not sparse.issparse(adata.layers["counts"]):
    raise RuntimeError("X and layers['counts'] must both be sparse")
log(f"Loaded source: shape={adata.shape}, X={type(adata.X).__name__}")

control_mask, perturbed_mask = normalize_perturbed(adata.obs["perturbed"])
if int(control_mask.sum() + perturbed_mask.sum()) != adata.n_obs:
    raise RuntimeError("Control and perturbed masks do not cover every cell")
log(
    f"Cells: control={control_mask.sum():,}, "
    f"perturbed={perturbed_mask.sum():,}, total={adata.n_obs:,}"
)

var_positions, selected_var, unmatched = select_target_genes(
    adata.obs, adata.var, perturbed_mask
)
log(
    f"Selected {len(selected_var):,} expression features; "
    f"excluded {len(unmatched):,} unmatched target genes"
)
if unmatched:
    log("Unmatched target genes: " + ", ".join(unmatched))

masks = {"control": control_mask, "perturbed": perturbed_mask}
expected_var_names = pd.Index(selected_var.index.astype(str))

for group, mask in masks.items():
    final_path = OUTPUTS[group]
    expected_obs_names = pd.Index(adata.obs_names[mask].astype(str))
    if output_is_valid(final_path, expected_obs_names, expected_var_names) and not OVERWRITE:
        log(f"Valid output already exists; skipping {group}: {final_path}")
        continue

    partial_path = final_path.with_name(final_path.name + ".partial")
    if partial_path.exists():
        partial_path.unlink()

    log(f"Creating {group} subset")
    subset = adata[mask, var_positions].copy()
    subset.obs["split_group"] = group
    subset.var = selected_var.copy()
    subset.uns["kolf_split"] = {
        "source_h5ad": str(INPUT),
        "group": group,
        "group_definition": f"obs['perturbed'] == {'False' if group == 'control' else 'True'}",
        "target_gene_definition": "genes targeted among perturbed CRISPRi cells",
        "unmatched_target_genes": unmatched,
    }
    log(f"Writing {group}: shape={subset.shape} -> {partial_path}")
    subset.write_h5ad(partial_path)
    validate_output(partial_path, expected_obs_names, expected_var_names)
    os.replace(partial_path, final_path)
    log(f"Validated {group} output: {final_path}")
    del subset
    gc.collect()

log("Split completed successfully")
PY

echo "Outputs written under: ${OUTPUT_DIR}"
