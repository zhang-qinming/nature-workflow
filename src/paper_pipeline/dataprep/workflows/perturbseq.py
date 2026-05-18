from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from ..config import ConfigError, LoadedConfig
from ..perturbseq.resources import r_script_path
from ..runner import preview_command
from ..tasks import Task


def _ensure_directories(paths: list[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def _load_gene_level_context(config: LoadedConfig) -> dict[str, object]:
    perturbseq = config.workflow("perturbseq")
    workflow = perturbseq.get("gene_level", {})
    if not workflow:
        raise ConfigError("workflows.perturbseq.gene_level is required")

    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})

    required_input_keys = ["h5ad", "metadata"]
    for key in required_input_keys:
        if key not in inputs:
            raise ConfigError(f"workflows.perturbseq.gene_level.inputs.{key} is required")

    prepare_inputs_script = r_script_path("prepare_limma_inputs.R")
    limma_script = r_script_path("limmatrend.R")
    summarize_script = r_script_path("summarize_effect_sizes.R")
    perturbseq_runner = config.executable_command("perturbseq_runner", [])
    rscript_command = config.executable_command("rscript", "Rscript")

    h5ad_path = config.resolve_path(inputs["h5ad"])
    metadata_path = config.resolve_path(inputs["metadata"])
    raw_count_csv = config.resolve_path_or_artifact(
        outputs.get("raw_count_csv"),
        "perturbseq",
        "gene_level",
        "K562gwps_raw_count.csv",
    )
    limma_input_dir = config.resolve_path_or_artifact(
        outputs.get("limma_input_dir"),
        "perturbseq",
        "gene_level",
        "limma_input",
    )
    limma_dir = config.resolve_path_or_artifact(
        outputs.get("limma_dir"),
        "perturbseq",
        "gene_level",
        "gwps_limma",
    )
    summary_dir = config.resolve_path_or_artifact(
        outputs.get("summary_dir"),
        "perturbseq",
        "gene_level",
        "K562GW",
    )

    assert h5ad_path and metadata_path and raw_count_csv

    chunk_size = int(parameters.get("chunk_size", 50))
    min_genes = int(parameters.get("filter_min_genes", 500))
    min_cells = int(parameters.get("filter_min_cells", 500))

    total_chunks = int(parameters.get("total_chunks", 1))
    if metadata_path.exists():
        metadata = pd.read_csv(metadata_path, index_col=0)
        genes = sorted(set(metadata["gene"].dropna()) - {"non-targeting"})
        total_chunks = math.ceil(len(genes) / chunk_size)
    configured_chunks = parameters.get("chunk_indices") or list(range(1, total_chunks + 1))
    chunk_indices = [int(value) for value in configured_chunks]

    return {
        "workflow": workflow,
        "prepare_inputs_script": prepare_inputs_script,
        "limma_script": limma_script,
        "summarize_script": summarize_script,
        "perturbseq_runner": perturbseq_runner,
        "rscript_command": rscript_command,
        "h5ad_path": h5ad_path,
        "metadata_path": metadata_path,
        "raw_count_csv": raw_count_csv,
        "limma_input_dir": limma_input_dir,
        "limma_dir": limma_dir,
        "summary_dir": summary_dir,
        "chunk_size": chunk_size,
        "min_genes": min_genes,
        "min_cells": min_cells,
        "chunk_indices": chunk_indices,
        "project_root": config.project_root,
    }


def build_perturbseq_gene_level_prepare_tasks(config: LoadedConfig) -> list[Task]:
    context = _load_gene_level_context(config)
    raw_count_csv = context["raw_count_csv"]
    limma_input_dir = context["limma_input_dir"]
    limma_dir = context["limma_dir"]
    summary_dir = context["summary_dir"]
    assert isinstance(raw_count_csv, Path)
    assert isinstance(limma_input_dir, Path)
    assert isinstance(limma_dir, Path)
    assert isinstance(summary_dir, Path)

    perturbseq_runner = context["perturbseq_runner"]
    rscript_command = context["rscript_command"]
    h5ad_path = context["h5ad_path"]
    metadata_path = context["metadata_path"]
    prepare_inputs_script = context["prepare_inputs_script"]
    chunk_size = context["chunk_size"]
    min_genes = context["min_genes"]
    min_cells = context["min_cells"]
    project_root = context["project_root"]

    assert isinstance(perturbseq_runner, list)
    assert isinstance(rscript_command, list)
    assert isinstance(h5ad_path, Path)
    assert isinstance(metadata_path, Path)
    assert isinstance(prepare_inputs_script, Path)
    assert isinstance(chunk_size, int)
    assert isinstance(min_genes, int)
    assert isinstance(min_cells, int)
    assert isinstance(project_root, Path)

    return [
        Task(
            name="Create Perturbseq gene-level output directories",
            preview=(
                f"Create output directories: {raw_count_csv.parent}, {limma_input_dir}, "
                f"{limma_dir}, {summary_dir}"
            ),
            action=lambda: _ensure_directories(
                [raw_count_csv.parent, limma_input_dir, limma_dir, summary_dir]
            ),
        ),
        Task(
            name="Extract Perturbseq raw count matrix from h5ad",
            preview=preview_command(
                perturbseq_runner
                + [
                    "paper-pipeline-perturbseq-extract-counts",
                    "--input-h5ad",
                    str(h5ad_path),
                    "--output-csv",
                    str(raw_count_csv),
                    "--min-genes-per-cell",
                    str(min_genes),
                    "--min-cells-per-gene",
                    str(min_cells),
                ],
                cwd=project_root,
            ),
            command=perturbseq_runner
            + [
                "paper-pipeline-perturbseq-extract-counts",
                "--input-h5ad",
                str(h5ad_path),
                "--output-csv",
                str(raw_count_csv),
                "--min-genes-per-cell",
                str(min_genes),
                "--min-cells-per-gene",
                str(min_cells),
            ],
            cwd=project_root,
        ),
        Task(
            name="Build limma input tables from metadata and raw counts",
            preview=preview_command(
                rscript_command
                + [
                    str(prepare_inputs_script),
                    str(metadata_path),
                    str(raw_count_csv),
                    str(limma_input_dir),
                    str(chunk_size),
                ],
                cwd=project_root,
            ),
            command=rscript_command
            + [
                str(prepare_inputs_script),
                str(metadata_path),
                str(raw_count_csv),
                str(limma_input_dir),
                str(chunk_size),
            ],
            cwd=project_root,
        ),
    ]


def build_perturbseq_gene_level_limma_tasks(config: LoadedConfig) -> list[Task]:
    context = _load_gene_level_context(config)
    limma_dir = context["limma_dir"]
    limma_input_dir = context["limma_input_dir"]
    limma_script = context["limma_script"]
    rscript_command = context["rscript_command"]
    chunk_indices = context["chunk_indices"]
    project_root = context["project_root"]

    assert isinstance(limma_dir, Path)
    assert isinstance(limma_input_dir, Path)
    assert isinstance(limma_script, Path)
    assert isinstance(rscript_command, list)
    assert isinstance(chunk_indices, list)
    assert isinstance(project_root, Path)

    tasks: list[Task] = [
        Task(
            name="Create Perturbseq limma output directory",
            preview=f"Create output directory: {limma_dir}",
            action=lambda: _ensure_directories([limma_dir]),
        )
    ]

    for chunk_index in chunk_indices:
        limma_command = rscript_command + [
            str(limma_script),
            str(chunk_index),
            str(limma_input_dir),
            str(limma_dir),
        ]
        tasks.append(
            Task(
                name=f"Run limma differential model for chunk {chunk_index}",
                preview=preview_command(limma_command, cwd=project_root),
                command=limma_command,
                cwd=project_root,
            )
        )
    return tasks


def build_perturbseq_gene_level_summarize_tasks(config: LoadedConfig) -> list[Task]:
    context = _load_gene_level_context(config)
    limma_dir = context["limma_dir"]
    summary_dir = context["summary_dir"]
    summarize_script = context["summarize_script"]
    rscript_command = context["rscript_command"]
    project_root = context["project_root"]

    assert isinstance(limma_dir, Path)
    assert isinstance(summary_dir, Path)
    assert isinstance(summarize_script, Path)
    assert isinstance(rscript_command, list)
    assert isinstance(project_root, Path)

    return [
        Task(
            name="Create Perturbseq gene-level summary directory",
            preview=f"Create output directory: {summary_dir}",
            action=lambda: _ensure_directories([summary_dir]),
        ),
        Task(
            name="Summarize limma outputs into gene-level matrices",
            preview=preview_command(
                rscript_command
                + [
                    str(summarize_script),
                    str(limma_dir),
                    str(summary_dir),
                ],
                cwd=project_root,
            ),
            command=rscript_command
            + [
                str(summarize_script),
                str(limma_dir),
                str(summary_dir),
            ],
            cwd=project_root,
        ),
    ]


def build_perturbseq_gene_level_tasks(config: LoadedConfig) -> list[Task]:
    return [
        *build_perturbseq_gene_level_prepare_tasks(config),
        *build_perturbseq_gene_level_limma_tasks(config)[1:],
        *build_perturbseq_gene_level_summarize_tasks(config)[1:],
    ]
