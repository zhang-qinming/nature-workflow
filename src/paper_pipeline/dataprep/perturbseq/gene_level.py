from __future__ import annotations

import argparse
from pathlib import Path


def extract_filtered_counts(
    input_h5ad: Path,
    output_csv: Path,
    *,
    min_genes_per_cell: int = 500,
    min_cells_per_gene: int = 500,
) -> None:
    import scanpy as sc

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(input_h5ad)
    sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    adata.to_df().to_csv(output_csv)


def main_extract_counts(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Extract filtered Perturbseq raw counts from an h5ad file."
    )
    parser.add_argument(
        "--input-h5ad",
        default="data/Perturbseq/K562_gwps_raw_singlecell_01.h5ad",
    )
    parser.add_argument(
        "--output-csv",
        default="data/Perturbseq/K562gwps_raw_count.csv",
    )
    parser.add_argument("--min-genes-per-cell", type=int, default=500)
    parser.add_argument("--min-cells-per-gene", type=int, default=500)
    args = parser.parse_args(argv)

    extract_filtered_counts(
        input_h5ad=Path(args.input_h5ad),
        output_csv=Path(args.output_csv),
        min_genes_per_cell=args.min_genes_per_cell,
        min_cells_per_gene=args.min_cells_per_gene,
    )
