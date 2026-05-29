#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOF_ID="${LOF_ID:?LOF_ID not set}"

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
FILE_ID_MAP="${FILE_ID_MAP:-${PROJECT_ROOT}/configs/path.file_id_map.tsv}"

TASK_NAME="${TASK_NAME:-cross_trait_heatmap}"
JOB_NAME="${JOB_NAME:-ctheat}"

OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/cross_trait_heatmap}"
RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
STATUS_DIR="${STATUS_DIR:-${RUN_ROOT}/status}"
FAILURE_DIR="${FAILURE_DIR:-${RUN_ROOT}/failed}"

POSTERIOR_DIR="${POSTERIOR_DIR:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/run_all/outputs/genebayes/posterior}"
POSTERIOR_NAME_MAP="${POSTERIOR_NAME_MAP:-}"
GENE_MAP="${GENE_MAP:-/gpfs/chencao/qinminzhang/Nature/mine/code/data/gencode_v41_gname_gid_ALL_sorted_onlyID}"

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
CONTROL_ENV="${CONTROL_ENV:-paper-pipeline-control}"

SUCCESS=0

mkdir -p "${STATUS_DIR}" "${FAILURE_DIR}"

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
}

trap on_exit EXIT

if [ -f "${CONDA_SH}" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
    conda activate "${CONTROL_ENV}"
fi

python3 - "${LOF_ID}" "${FILE_ID_MAP}" "${POSTERIOR_DIR}" "${POSTERIOR_NAME_MAP}" "${GENE_MAP}" "${OUTPUT_DIR}" <<'PYEOF'
import csv
import glob
import sys
from pathlib import Path

lof_id = sys.argv[1]
file_id_map = Path(sys.argv[2])
posterior_dir = Path(sys.argv[3])
posterior_name_map = sys.argv[4].strip()
gene_map_path = Path(sys.argv[5])
output_dir = Path(sys.argv[6])


def derive_posterior_name(source_path: str) -> str:
    stem = Path(source_path).name
    for suffix in (
        ".summary_statistics.csv",
        ".csv.gz",
        ".tsv.gz",
        ".txt.gz",
        ".csv",
        ".tsv",
    ):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return f"{stem}.per_gene_estimates.tsv"


def load_name_overrides(path: str) -> dict[str, str]:
    if not path:
        return {}
    override_path = Path(path)
    if not override_path.exists():
        raise SystemExit(f"POSTERIOR_NAME_MAP does not exist: {override_path}")
    overrides: dict[str, str] = {}
    with override_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_index, row in enumerate(reader, start=1):
            if len(row) < 2:
                continue
            if row_index == 1 and row[0] in {"id", "id2", "lof_id"}:
                continue
            overrides[row[0].strip()] = row[1].strip()
    return overrides


def load_gene_lookup(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    lookup: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            ensg = row[0].split(".")[0].strip()
            gene = row[1].strip()
            if ensg and gene:
                lookup[ensg] = gene
    return lookup


map_row = None
with file_id_map.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"id1", "id2", "path1", "path2"}
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise SystemExit(f"{file_id_map} must contain columns: id1, id2, path1, path2")
    for row in reader:
        if row["id2"].strip() == lof_id:
            map_row = row
            break

if map_row is None:
    raise SystemExit(f"LoF ID not found in file_id_map id2 column: {lof_id}")

overrides = load_name_overrides(posterior_name_map)
posterior_name = overrides.get(lof_id) or derive_posterior_name(map_row["path2"].strip())
posterior_path = posterior_dir / posterior_name
if not posterior_path.exists():
    matches = sorted(glob.glob(str(posterior_dir / f"*{lof_id}*.per_gene_estimates.tsv")))
    if matches:
        posterior_path = Path(matches[0])
    else:
        raise SystemExit(f"Posterior file not found for {lof_id}: {posterior_path}")

gene_lookup = load_gene_lookup(gene_map_path)

effects_dir = output_dir / "tables" / "effects"
meta_dir = output_dir / "meta" / "traits"
effects_dir.mkdir(parents=True, exist_ok=True)
meta_dir.mkdir(parents=True, exist_ok=True)

effect_path = effects_dir / f"{lof_id}.tsv"
meta_path = meta_dir / f"{lof_id}.tsv"

with posterior_path.open(encoding="utf-8", newline="") as handle, effect_path.open("w", encoding="utf-8", newline="") as out:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"Posterior file has no header: {posterior_path}")
    missing = {"ensg", "post_mean"} - set(reader.fieldnames)
    if missing:
        raise SystemExit(f"Posterior file {posterior_path} is missing columns: {', '.join(sorted(missing))}")
    writer = csv.DictWriter(out, fieldnames=["ensg", "gene", "trait_id", "post_mean"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        ensg = (row.get("ensg") or "").split(".")[0].strip()
        if not ensg:
            continue
        writer.writerow(
            {
                "ensg": ensg,
                "gene": gene_lookup.get(ensg, ensg),
                "trait_id": lof_id,
                "post_mean": row.get("post_mean", ""),
            }
        )

with meta_path.open("w", encoding="utf-8", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=["source_id", "gwas_id", "source_path", "posterior_path", "effect_path"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerow(
        {
            "source_id": lof_id,
            "gwas_id": map_row["id1"].strip(),
            "source_path": map_row["path2"].strip(),
            "posterior_path": str(posterior_path),
            "effect_path": str(effect_path),
        }
    )
PYEOF

printf 'time=%s\njob_id=%s\n' "$(date -Iseconds)" "${SLURM_JOB_ID:-local}" > "${STATUS_DIR}/${LOF_ID}.ok"
rm -f "${STATUS_DIR}/${LOF_ID}.failed" "${STATUS_DIR}/${LOF_ID}.running" "${FAILURE_DIR}/${LOF_ID}.failed"
SUCCESS=1
