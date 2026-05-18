from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def filter_cells_for_cnmf(
    input_h5ad: Path,
    output_h5ad: Path,
    *,
    mito_column: str = "mitopercent",
    mito_threshold: float = 0.3,
    min_genes_per_cell: int = 100,
    min_counts_per_cell: int = 100,
    min_cells_per_gene: int = 10,
) -> None:
    import scanpy as sc

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(input_h5ad)
    if mito_column in adata.obs:
        adata = adata[adata.obs[mito_column] <= mito_threshold].copy()
    sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
    sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    adata.write(output_h5ad)


def merge_cnmf_tmp_runs(
    source_dirs: list[Path],
    source_names: list[str],
    target_dir: Path,
    *,
    target_name: str,
) -> None:
    if len(source_dirs) != len(source_names):
        raise ValueError("source_dirs and source_names must have the same length")

    target_tmp = target_dir / "cnmf_tmp"
    target_tmp.mkdir(parents=True, exist_ok=True)
    preserve_target_files = ("nmf_idvrun_params.yaml", "nmf_params.df.npz")

    for source_dir, source_name in zip(source_dirs, source_names, strict=True):
        source_tmp = source_dir / "cnmf_tmp"
        if not source_tmp.exists():
            raise FileNotFoundError(source_tmp)

        for path in source_tmp.iterdir():
            new_name = path.name.replace(source_name, target_name)
            if new_name == path.name and source_name != target_name:
                new_name = f"{target_name}.{path.name}"
            if any(new_name.endswith(suffix) for suffix in preserve_target_files):
                continue
            shutil.copy2(path, target_tmp / new_name)


def main_filter_cells(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Filter Perturbseq h5ad inputs before cNMF preparation."
    )
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--mito-column", default="mitopercent")
    parser.add_argument("--mito-threshold", type=float, default=0.3)
    parser.add_argument("--min-genes-per-cell", type=int, default=100)
    parser.add_argument("--min-counts-per-cell", type=int, default=100)
    parser.add_argument("--min-cells-per-gene", type=int, default=10)
    args = parser.parse_args(argv)

    filter_cells_for_cnmf(
        input_h5ad=Path(args.input_h5ad),
        output_h5ad=Path(args.output_h5ad),
        mito_column=args.mito_column,
        mito_threshold=args.mito_threshold,
        min_genes_per_cell=args.min_genes_per_cell,
        min_counts_per_cell=args.min_counts_per_cell,
        min_cells_per_gene=args.min_cells_per_gene,
    )


def main_merge_tmp(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Merge per-K cNMF temporary outputs into an aggregate run directory."
    )
    parser.add_argument("--source-dir", action="append", dest="source_dirs", required=True)
    parser.add_argument("--source-name", action="append", dest="source_names", required=True)
    parser.add_argument("--target-dir", required=True)
    parser.add_argument("--target-name", default="cNMF_all")
    args = parser.parse_args(argv)

    merge_cnmf_tmp_runs(
        source_dirs=[Path(value) for value in args.source_dirs],
        source_names=[str(value) for value in args.source_names],
        target_dir=Path(args.target_dir),
        target_name=str(args.target_name),
    )
