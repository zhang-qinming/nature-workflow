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
MEMORY="${MEMORY:-300G}"
CPUS="${CPUS:-8}"
TIME_LIMIT="${TIME_LIMIT:-7-00:00:00}"
MIN_FREE_GB="${MIN_FREE_GB:-100}"
COMPRESSION_LEVEL="${COMPRESSION_LEVEL:-1}"
OVERWRITE="${OVERWRITE:-0}"

usage() {
    cat <<'EOF'
Usage:
  bash split_kolf_control_perturbed_h5ad.sh          Submit itself with sbatch.
  bash split_kolf_control_perturbed_h5ad.sh --run    Run in the current allocation.

Environment overrides:
  PROJECT_ROOT INPUT_H5AD OUTPUT_DIR CONDA_SH CONDA_ENV
  PARTITION MEMORY CPUS TIME_LIMIT MIN_FREE_GB
  COMPRESSION_LEVEL OVERWRITE

Outputs:
  KOLF.control.target_genes.h5ad
  KOLF.perturbed.target_genes.h5ad
  control_cells.tsv.gz
  perturbed_cells.tsv.gz
  target_genes.tsv
  gene_matching_report.tsv
  split_summary.tsv
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
    job_id="$({ sbatch \
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
        "${SCRIPT_PATH}" --run; } | cut -d';' -f1)"
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
[[ "${COMPRESSION_LEVEL}" =~ ^[1-9]$ ]] || {
    echo "COMPRESSION_LEVEL must be an integer from 1 to 9" >&2
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
    "${COMPRESSION_LEVEL}" \
    "${OVERWRITE}" <<'PY'
from __future__ import annotations

import gc
import os
import re
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


INPUT = Path(sys.argv[1]).resolve()
OUTPUT_DIR = Path(sys.argv[2]).resolve()
MIN_FREE_GB = int(sys.argv[3])
COMPRESSION_LEVEL = int(sys.argv[4])
OVERWRITE = bool(int(sys.argv[5]))

OUTPUTS = {
    "control": OUTPUT_DIR / "KOLF.control.target_genes.h5ad",
    "perturbed": OUTPUT_DIR / "KOLF.perturbed.target_genes.h5ad",
}


def log(message: str) -> None:
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), message, flush=True)


def atomic_to_csv(frame: pd.DataFrame, path: Path, **kwargs: object) -> None:
    suffix = ".tmp.gz" if path.suffix == ".gz" else ".tmp"
    temporary = path.with_name(path.name + suffix)
    frame.to_csv(temporary, **kwargs)
    os.replace(temporary, path)


def normalize_perturbed(values: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    normalized = values.astype("string").str.strip().str.lower()
    true_mask = normalized.isin({"true", "1"}).to_numpy(dtype=bool)
    false_mask = normalized.isin({"false", "0"}).to_numpy(dtype=bool)
    unknown = ~(true_mask | false_mask)
    if unknown.any():
        examples = sorted(normalized[unknown].fillna("<NA>").unique().tolist())[:20]
        raise RuntimeError(
            "obs['perturbed'] contains values other than True/False: "
            + ", ".join(map(str, examples))
        )
    return false_mask, true_mask


def strip_ensembl_version(value: object) -> str:
    text = str(value).strip()
    return text.split(".", 1)[0]


def split_identifiers(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [
        strip_ensembl_version(item)
        for item in re.split(r"[;,|]", str(value))
        if item.strip()
    ]


def build_gene_selection(
    obs: pd.DataFrame, var: pd.DataFrame, perturbed_mask: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
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

    name_to_position = {name: pos for pos, name in enumerate(var_names)}
    id_to_positions: dict[str, list[int]] = {}
    if "gene_ids" in var.columns:
        for pos, value in enumerate(var["gene_ids"]):
            if pd.isna(value):
                continue
            gene_id = strip_ensembl_version(value)
            id_to_positions.setdefault(gene_id, []).append(pos)

    ids_by_target: dict[str, set[str]] = {}
    for target, value in pairs.itertuples(index=False, name=None):
        ids_by_target.setdefault(str(target), set()).update(split_identifiers(value))

    report_rows: list[dict[str, object]] = []
    position_to_targets: dict[int, list[str]] = {}
    for target, cell_count in target_counts.items():
        target = str(target)
        position: int | None = None
        method = ""
        note = ""

        if target in name_to_position:
            position = name_to_position[target]
            method = "var_name"
        else:
            candidates: set[int] = set()
            candidate_ids = set(ids_by_target.get(target, set()))
            if target.startswith("ENSG"):
                candidate_ids.add(strip_ensembl_version(target))
            for gene_id in candidate_ids:
                candidates.update(id_to_positions.get(gene_id, []))
            if len(candidates) == 1:
                position = next(iter(candidates))
                method = "gene_target_ensembl_id"
            elif len(candidates) > 1:
                method = "ambiguous"
                note = "Multiple expression features match the target Ensembl ID"
            else:
                method = "unmatched"
                note = "No matching var_name or gene_ids value"

        if position is not None:
            position_to_targets.setdefault(position, []).append(target)
            matched_var_name = var_names[position]
            matched_gene_id = (
                str(var.iloc[position]["gene_ids"])
                if "gene_ids" in var.columns
                else ""
            )
        else:
            matched_var_name = ""
            matched_gene_id = ""

        report_rows.append(
            {
                "target_gene": target,
                "perturbed_cell_count": int(cell_count),
                "target_ensembl_ids": ";".join(sorted(ids_by_target.get(target, set()))),
                "status": "matched" if position is not None else method,
                "match_method": method,
                "var_name": matched_var_name,
                "var_gene_id": matched_gene_id,
                "note": note,
            }
        )

    report = pd.DataFrame(report_rows)
    match_rate = (report["status"] == "matched").mean()
    if match_rate < 0.99:
        raise RuntimeError(
            f"Only {match_rate:.2%} of target genes matched the expression matrix; "
            "inspect gene identifiers before continuing"
        )

    positions = np.array(sorted(position_to_targets), dtype=np.int64)
    selected_var = var.iloc[positions].copy()
    selected_var["crispri_target"] = True
    selected_var["crispri_target_symbols"] = [
        ";".join(sorted(position_to_targets[int(pos)])) for pos in positions
    ]
    return positions, selected_var, report


def write_empty_mapping(handle: h5py.File, key: str) -> None:
    ad.io.write_elem(handle, key, {})


@dataclass
class SparseCSCAppender:
    group: h5py.Group
    expected_rows: int
    expected_columns: int
    columns_written: int = 0
    nnz_written: int = 0
    buffered_nnz: int = 0

    def __post_init__(self) -> None:
        self._data_parts: list[np.ndarray] = []
        self._index_parts: list[np.ndarray] = []
        self._pointer_buffer: list[int] = []

    @classmethod
    def create(
        cls,
        parent: h5py.Group | h5py.File,
        key: str,
        shape: tuple[int, int],
        data_dtype: np.dtype,
    ) -> "SparseCSCAppender":
        group = parent.create_group(key)
        group.attrs["encoding-type"] = "csc_matrix"
        group.attrs["encoding-version"] = "0.1.0"
        group.attrs["shape"] = np.asarray(shape, dtype=np.int64)
        common = {
            "maxshape": (None,),
            "chunks": (1_000_000,),
            "compression": "gzip",
            "compression_opts": COMPRESSION_LEVEL,
        }
        group.create_dataset("data", shape=(0,), dtype=data_dtype, **common)
        group.create_dataset("indices", shape=(0,), dtype=np.int32, **common)
        group.create_dataset(
            "indptr",
            data=np.array([0], dtype=np.int64),
            maxshape=(None,),
            chunks=(262_144,),
            compression="gzip",
            compression_opts=COMPRESSION_LEVEL,
        )
        return cls(group=group, expected_rows=shape[0], expected_columns=shape[1])

    def append_column(self, data: np.ndarray, row_indices: np.ndarray) -> None:
        data = np.asarray(data)
        row_indices = np.asarray(row_indices, dtype=np.int32)
        if data.shape[0] != row_indices.shape[0]:
            raise RuntimeError("CSC data and row-index lengths differ")
        if row_indices.size and (
            int(row_indices[0]) < 0 or int(row_indices[-1]) >= self.expected_rows
        ):
            raise RuntimeError("A remapped CSC row index is out of range")
        if row_indices.size > 1 and np.any(row_indices[1:] < row_indices[:-1]):
            raise RuntimeError("CSC row indices are not sorted")

        count = int(data.shape[0])
        if count:
            self._data_parts.append(data)
            self._index_parts.append(row_indices)
        self.buffered_nnz += count
        self.columns_written += 1
        self._pointer_buffer.append(self.nnz_written + self.buffered_nnz)
        if self.buffered_nnz >= 5_000_000 or len(self._pointer_buffer) >= 256:
            self.flush()

    def flush(self) -> None:
        if not self._pointer_buffer:
            return
        if self.buffered_nnz:
            data = np.concatenate(self._data_parts)
            indices = np.concatenate(self._index_parts)
            for key, values in (("data", data), ("indices", indices)):
                dataset = self.group[key]
                dataset.resize((self.nnz_written + self.buffered_nnz,))
                dataset[self.nnz_written :] = values
        indptr = self.group["indptr"]
        old_length = indptr.shape[0]
        pointers = np.asarray(self._pointer_buffer, dtype=np.int64)
        indptr.resize((old_length + len(pointers),))
        indptr[old_length:] = pointers

        self.nnz_written += self.buffered_nnz
        self.buffered_nnz = 0
        self._data_parts.clear()
        self._index_parts.clear()
        self._pointer_buffer.clear()

    def finish(self) -> None:
        self.flush()
        if self.columns_written != self.expected_columns:
            raise RuntimeError(
                f"Wrote {self.columns_written} matrix columns; "
                f"expected {self.expected_columns}"
            )
        if self.group["indptr"].shape[0] != self.expected_columns + 1:
            raise RuntimeError("Invalid CSC indptr length")


class GroupOutput:
    def __init__(
        self,
        group_name: str,
        final_path: Path,
        obs: pd.DataFrame,
        var: pd.DataFrame,
        x_dtype: np.dtype,
        counts_dtype: np.dtype,
    ) -> None:
        self.group_name = group_name
        self.final_path = final_path
        self.partial_path = final_path.with_name(final_path.name + ".partial")
        if self.partial_path.exists():
            self.partial_path.unlink()
        self.handle = h5py.File(self.partial_path, "w")
        self.handle.attrs["encoding-type"] = "anndata"
        self.handle.attrs["encoding-version"] = "0.1.0"
        ad.io.write_elem(self.handle, "obs", obs)
        ad.io.write_elem(self.handle, "var", var)
        for key in ("obsm", "varm", "obsp", "varp", "uns"):
            write_empty_mapping(self.handle, key)
        write_empty_mapping(self.handle, "layers")

        shape = (len(obs), len(var))
        self.x = SparseCSCAppender.create(self.handle, "X", shape, x_dtype)
        self.counts = SparseCSCAppender.create(
            self.handle["layers"], "counts", shape, counts_dtype
        )

    def finish(self) -> None:
        self.x.finish()
        self.counts.finish()
        self.handle.flush()
        self.handle.close()
        os.replace(self.partial_path, self.final_path)

    def abort(self) -> None:
        try:
            self.handle.close()
        finally:
            if self.partial_path.exists():
                self.partial_path.unlink()


def validate_output(
    path: Path,
    expected_obs_names: pd.Index,
    expected_var_names: pd.Index,
) -> None:
    with h5py.File(path, "r") as handle:
        expected_shape = (len(expected_obs_names), len(expected_var_names))
        if "X" not in handle:
            raise RuntimeError(f"{path} does not contain X")
        if "layers" not in handle or "counts" not in handle["layers"]:
            raise RuntimeError(f"{path} does not contain layers['counts']")
        for key in ("X", "layers/counts"):
            matrix = handle[key]
            shape = tuple(int(value) for value in matrix.attrs.get("shape", ()))
            if shape != expected_shape:
                raise RuntimeError(
                    f"{path}:{key} has shape {shape}; expected {expected_shape}"
                )
            if matrix.attrs.get("encoding-type") != "csc_matrix":
                raise RuntimeError(f"{path}:{key} is not encoded as CSC")
            if matrix["indptr"].shape[0] != expected_shape[1] + 1:
                raise RuntimeError(f"{path}:{key} has an invalid indptr length")
            if matrix["data"].shape[0] != matrix["indices"].shape[0]:
                raise RuntimeError(f"{path}:{key} data/index lengths differ")
            if int(matrix["indptr"][-1]) != matrix["data"].shape[0]:
                raise RuntimeError(f"{path}:{key} has an invalid final pointer")

        result_obs = ad.io.read_elem(handle["obs"])
        result_var = ad.io.read_elem(handle["var"])
        if not np.array_equal(
            result_obs.index.astype(str).to_numpy(),
            expected_obs_names.astype(str).to_numpy(),
        ):
            raise RuntimeError(f"{path} cell order does not match the source H5AD")
        if not np.array_equal(
            result_var.index.astype(str).to_numpy(),
            expected_var_names.astype(str).to_numpy(),
        ):
            raise RuntimeError(f"{path} gene order does not match the selected genes")


def source_matrix_metadata(
    handle: h5py.File, key: str, expected_shape: tuple[int, int]
) -> tuple[h5py.Group, np.ndarray]:
    if key not in handle:
        raise RuntimeError(f"Input H5AD is missing {key}")
    matrix = handle[key]
    encoding = matrix.attrs.get("encoding-type")
    if encoding != "csc_matrix":
        raise RuntimeError(
            f"Input {key} uses {encoding!r}; this script requires csc_matrix"
        )
    shape = tuple(int(value) for value in matrix.attrs.get("shape", ()))
    if shape != expected_shape:
        raise RuntimeError(f"Input {key} has shape {shape}; expected {expected_shape}")
    pointers = np.asarray(matrix["indptr"][:], dtype=np.int64)
    if pointers.shape[0] != expected_shape[1] + 1:
        raise RuntimeError(f"Input {key} has an invalid CSC indptr length")
    if int(pointers[-1]) != matrix["data"].shape[0]:
        raise RuntimeError(f"Input {key} has an invalid final CSC pointer")
    return matrix, pointers


def append_source_matrix(
    source_matrix: h5py.Group,
    pointers: np.ndarray,
    var_positions: np.ndarray,
    masks: dict[str, np.ndarray],
    row_maps: dict[str, np.ndarray],
    writers: dict[str, GroupOutput],
    attribute: str,
    label: str,
) -> None:
    data_dataset = source_matrix["data"]
    index_dataset = source_matrix["indices"]
    total = len(var_positions)
    for output_column, source_column in enumerate(var_positions, start=1):
        start = int(pointers[source_column])
        end = int(pointers[source_column + 1])
        rows = np.asarray(index_dataset[start:end], dtype=np.int64)
        values = np.asarray(data_dataset[start:end])
        if rows.size > 1 and np.any(rows[1:] < rows[:-1]):
            raise RuntimeError(
                f"Input {label} column {source_column} has unsorted row indices"
            )
        for group, writer in writers.items():
            keep = masks[group][rows]
            output_rows = row_maps[group][rows[keep]]
            getattr(writer, attribute).append_column(values[keep], output_rows)

        if output_column == 1 or output_column % 500 == 0 or output_column == total:
            log(f"{label}: processed target genes {output_column:,}/{total:,}")


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
                "Set OVERWRITE=1 to replace it after checking the path."
            ) from exc
        return False
    return True


def selected_cell_table(obs: pd.DataFrame, mask: np.ndarray, group: str) -> pd.DataFrame:
    preferred = [
        "gene_target",
        "gene_target_ensembl_id",
        "gRNA",
        "perturbation",
        "perturbation_type",
        "batch",
        "channel",
        "celltype",
        "n_genes",
        "total_counts",
        "pct_counts_mt",
        "outlier",
    ]
    columns = [column for column in preferred if column in obs.columns]
    table = obs.loc[mask, columns].copy()
    table.insert(0, "split_group", group)
    table.index = table.index.astype(str)
    table.index.name = "cell_barcode"
    return table


log(f"Opening source metadata without loading expression matrices: {INPUT}")
with h5py.File(INPUT, "r") as source:
    if "obs" not in source or "var" not in source:
        raise RuntimeError("Input H5AD is missing obs or var")
    obs = ad.io.read_elem(source["obs"])
    var = ad.io.read_elem(source["var"])
    n_obs = len(obs)
    n_vars = len(var)

    if "perturbed" not in obs.columns:
        raise RuntimeError("Input H5AD is missing obs['perturbed']")
    x_source, x_pointers = source_matrix_metadata(source, "X", (n_obs, n_vars))
    counts_source, counts_pointers = source_matrix_metadata(
        source, "layers/counts", (n_obs, n_vars)
    )

    control_mask, perturbed_mask = normalize_perturbed(obs["perturbed"])
    if int(control_mask.sum() + perturbed_mask.sum()) != n_obs:
        raise RuntimeError("Control and perturbed masks do not cover every cell")

    log(
        f"Cells: control={control_mask.sum():,}, "
        f"perturbed={perturbed_mask.sum():,}, total={n_obs:,}"
    )
    var_positions, selected_var, match_report = build_gene_selection(
        obs, var, perturbed_mask
    )
    log(
        f"Target genes: matched={(match_report.status == 'matched').sum():,}, "
        f"unmatched={(match_report.status != 'matched').sum():,}, "
        f"expression features={len(var_positions):,}"
    )

    masks = {"control": control_mask, "perturbed": perturbed_mask}
    expected_obs_names = {
        group: pd.Index(obs.index[mask].astype(str)) for group, mask in masks.items()
    }
    expected_var_names = pd.Index(selected_var.index.astype(str))

    free_gb = shutil.disk_usage(OUTPUT_DIR).free / 1024**3
    log(f"Free disk space: {free_gb:.1f} GiB")
    if free_gb < MIN_FREE_GB:
        raise RuntimeError(
            f"Only {free_gb:.1f} GiB is free under {OUTPUT_DIR}; "
            f"at least {MIN_FREE_GB} GiB is required"
        )

    active_groups: list[str] = []
    for group, path in OUTPUTS.items():
        valid = output_is_valid(path, expected_obs_names[group], expected_var_names)
        if valid and not OVERWRITE:
            log(f"Valid output already exists; skipping {group}: {path}")
        else:
            active_groups.append(group)

    atomic_to_csv(
        match_report,
        OUTPUT_DIR / "gene_matching_report.tsv",
        sep="\t",
        index=False,
    )
    target_table = selected_var.copy()
    target_table.index = target_table.index.astype(str)
    target_table.index.name = "var_name"
    atomic_to_csv(target_table, OUTPUT_DIR / "target_genes.tsv", sep="\t")
    for group, mask in masks.items():
        atomic_to_csv(
            selected_cell_table(obs, mask, group),
            OUTPUT_DIR / f"{group}_cells.tsv.gz",
            sep="\t",
            compression="gzip",
        )

    if active_groups:
        x_dtype = np.dtype(x_source["data"].dtype)
        counts_dtype = np.dtype(counts_source["data"].dtype)
        writers: dict[str, GroupOutput] = {}
        try:
            for group in active_groups:
                group_obs = obs.loc[masks[group]].copy()
                group_obs["split_group"] = group
                writers[group] = GroupOutput(
                    group,
                    OUTPUTS[group],
                    group_obs,
                    selected_var,
                    x_dtype,
                    counts_dtype,
                )
                del group_obs
                gc.collect()

            row_maps: dict[str, np.ndarray] = {}
            for group in active_groups:
                row_map = np.full(n_obs, -1, dtype=np.int32)
                row_map[masks[group]] = np.arange(masks[group].sum(), dtype=np.int32)
                row_maps[group] = row_map

            append_source_matrix(
                x_source,
                x_pointers,
                var_positions,
                masks,
                row_maps,
                writers,
                "x",
                "X",
            )
            append_source_matrix(
                counts_source,
                counts_pointers,
                var_positions,
                masks,
                row_maps,
                writers,
                "counts",
                "layers/counts",
            )

            for writer in writers.values():
                writer.finish()
        except BaseException:
            for writer in writers.values():
                try:
                    writer.abort()
                except Exception:
                    pass
            raise

    for group, path in OUTPUTS.items():
        validate_output(path, expected_obs_names[group], expected_var_names)
        log(f"Validated {group} output: {path}")

    summary = pd.DataFrame(
        [
            {
                "group": group,
                "cells": int(mask.sum()),
                "target_expression_features": len(selected_var),
                "output_h5ad": str(OUTPUTS[group]),
                "output_bytes": OUTPUTS[group].stat().st_size,
            }
            for group, mask in masks.items()
        ]
    )
    atomic_to_csv(summary, OUTPUT_DIR / "split_summary.tsv", sep="\t", index=False)
    log("Split completed successfully")
PY

echo "Outputs written under: ${OUTPUT_DIR}"
