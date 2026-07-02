from __future__ import annotations

from pathlib import Path

from ..config import ConfigError, LoadedConfig
from ..perturbseq.resources import r_script_path
from ..runner import preview_command
from ..tasks import Task


def _ensure_directories(paths: list[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def _threshold_label(value: float) -> str:
    return f"{value:g}".replace(".", "_")


def _trait_files(config: dict, lof_dir: Path) -> list[str]:
    explicit = config.get("trait_files")
    if explicit:
        return [str(value) for value in explicit]

    trait_glob = str(config.get("trait_glob", "*.per_gene_estimates.tsv"))
    matches = sorted(path.name for path in lof_dir.glob(trait_glob))
    if not matches:
        raise ConfigError(f"No trait files found in {lof_dir} matching {trait_glob}")
    return matches


def _command_with_runner(runner: list[str], *parts: str) -> list[str]:
    return [*runner, *parts]


def _parse_worker_indices(parameters: dict, total_workers: int, config_key: str) -> list[int]:
    raw = parameters.get("worker_indices")
    if raw is None:
        values = list(range(total_workers))
    else:
        values = [int(value) for value in raw]
    if not values:
        raise ConfigError(f"{config_key} must not be empty")

    invalid = [value for value in values if value < 0 or value >= total_workers]
    if invalid:
        raise ConfigError(
            f"{config_key} contains invalid worker indices {invalid}; valid range is 0..{total_workers - 1}"
        )

    deduped: list[int] = []
    seen: set[int] = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        deduped.append(value)
    return deduped


def _require_full_worker_set(
    worker_indices: list[int],
    total_workers: int,
    *,
    workflow_name: str,
    factorize_workflow_name: str,
    post_workflow_name: str,
) -> None:
    full_set = set(range(total_workers))
    current_set = set(worker_indices)
    if current_set == full_set and len(worker_indices) == total_workers:
        return
    raise ConfigError(
        f"`{workflow_name}` requires `worker_indices` to cover every worker 0..{total_workers - 1}. "
        f"Current value is {worker_indices}. For multi-node execution, run "
        f"`{factorize_workflow_name}` with subsets of workers, then run `{post_workflow_name}` "
        "after all factorization workers have completed."
    )


def _parse_program_indices(parameters: dict, total_programs: int, config_key: str) -> list[int]:
    raw = parameters.get("program_indices")
    if raw is None:
        values = list(range(1, total_programs + 1))
    else:
        values = [int(value) for value in raw]
    if not values:
        raise ConfigError(f"{config_key} must not be empty")

    invalid = [value for value in values if value < 1 or value > total_programs]
    if invalid:
        raise ConfigError(
            f"{config_key} contains invalid program indices {invalid}; valid range is 1..{total_programs}"
        )

    deduped: list[int] = []
    seen: set[int] = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        deduped.append(value)
    return deduped


def _require_full_program_set(
    program_indices: list[int],
    total_programs: int,
    *,
    workflow_name: str,
    regulation_workflow_name: str,
) -> None:
    full_set = set(range(1, total_programs + 1))
    current_set = set(program_indices)
    if current_set == full_set and len(program_indices) == total_programs:
        return
    raise ConfigError(
        f"`{workflow_name}` requires `program_indices` to cover every program 1..{total_programs}. "
        f"Current value is {program_indices}. For multi-node execution, run "
        f"`{regulation_workflow_name}` with subsets of programs after the corresponding "
        "postbase base workflow has completed."
    )


def _load_essential_context(config: LoadedConfig) -> dict:
    perturbseq = config.workflow("perturbseq")
    workflow = perturbseq.get("cnmf_essential", {})
    if not workflow:
        raise ConfigError("workflows.perturbseq.cnmf_essential is required")

    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})
    cell_types = [str(value) for value in workflow.get("cell_types", [])]
    if not cell_types:
        raise ConfigError("workflows.perturbseq.cnmf_essential.cell_types must not be empty")

    required_input_keys = ["h5ad_dir", "metadata_dir", "gene_map", "shet_path"]
    for key in required_input_keys:
        if key not in inputs:
            raise ConfigError(f"workflows.perturbseq.cnmf_essential.inputs.{key} is required")

    perturbseq_runner = config.executable_command("perturbseq_runner", [])
    rscript_runner = config.executable_command("rscript", "Rscript")

    h5ad_dir = config.resolve_path(inputs["h5ad_dir"])
    metadata_dir = config.resolve_path(inputs["metadata_dir"])
    gene_map = config.resolve_path(inputs["gene_map"])
    lof_dir = config.resolve_path_or_artifact(inputs.get("lof_dir"), "genebayes", "posterior")
    shet_path = config.resolve_path(inputs["shet_path"])

    filtered_h5ad_dir = config.resolve_path_or_artifact(
        outputs.get("filtered_h5ad_dir"),
        "perturbseq",
        "cnmf_essential",
        "filtered_data",
    )
    cnmf_dir = config.resolve_path_or_artifact(
        outputs.get("cnmf_dir"),
        "perturbseq",
        "cnmf_essential",
        "cNMF",
    )
    regulation_dir = config.resolve_path_or_artifact(
        outputs.get("regulation_dir"),
        "perturbseq",
        "cnmf_essential",
        "cNMF_regulation",
    )
    association_dir = config.resolve_path_or_artifact(
        outputs.get("association_dir"),
        "perturbseq",
        "cnmf_essential",
        "trait_association",
    )

    assert h5ad_dir and metadata_dir and gene_map and lof_dir and shet_path

    k = int(parameters.get("k", 60))
    total_workers = int(parameters.get("total_workers", 20))
    worker_indices = _parse_worker_indices(
        parameters,
        total_workers,
        "workflows.perturbseq.cnmf_essential.parameters.worker_indices",
    )
    n_iter = int(parameters.get("n_iter", 100))
    seed = int(parameters.get("seed", 14))
    numgenes = int(parameters.get("numgenes", 2000))
    density_threshold = float(parameters.get("consensus_density_threshold", 0.4))
    mito_column = str(parameters.get("mito_column", "mitopercent"))
    mito_threshold = float(parameters.get("mito_threshold", 0.3))
    min_genes = int(parameters.get("min_genes_per_cell", 100))
    min_counts = int(parameters.get("min_counts_per_cell", 100))
    min_cells = int(parameters.get("min_cells_per_gene", 10))
    random_iterations = int(parameters.get("random_iterations", 100000))
    top_n = int(parameters.get("top_n_program_genes", 100))
    trait_files = _trait_files(workflow, lof_dir)
    threshold = _threshold_label(density_threshold)
    program_indices = _parse_program_indices(
        parameters,
        k,
        "workflows.perturbseq.cnmf_essential.parameters.program_indices",
    )

    return {
        "workflow": workflow,
        "perturbseq_runner": perturbseq_runner,
        "rscript_runner": rscript_runner,
        "cell_types": cell_types,
        "h5ad_dir": h5ad_dir,
        "metadata_dir": metadata_dir,
        "gene_map": gene_map,
        "lof_dir": lof_dir,
        "shet_path": shet_path,
        "filtered_h5ad_dir": filtered_h5ad_dir,
        "cnmf_dir": cnmf_dir,
        "regulation_dir": regulation_dir,
        "association_dir": association_dir,
        "k": k,
        "total_workers": total_workers,
        "worker_indices": worker_indices,
        "n_iter": n_iter,
        "seed": seed,
        "numgenes": numgenes,
        "density_threshold": density_threshold,
        "threshold": threshold,
        "mito_column": mito_column,
        "mito_threshold": mito_threshold,
        "min_genes": min_genes,
        "min_counts": min_counts,
        "min_cells": min_cells,
        "random_iterations": random_iterations,
        "top_n": top_n,
        "trait_files": trait_files,
        "program_indices": program_indices,
        "reg_script": r_script_path("cnmf_regulatory_effects.R"),
        "burden_script": r_script_path("cnmf_burden_program_regulators.R"),
    }


def _build_essential_prepare_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = [
        Task(
            name="Create Perturbseq cNMF essential output directories",
            preview=(
                "Create output directories: "
                f"{ctx['filtered_h5ad_dir']}, {ctx['cnmf_dir']}, {ctx['regulation_dir']}, {ctx['association_dir']}"
            ),
            action=lambda: _ensure_directories(
                [ctx["filtered_h5ad_dir"], ctx["cnmf_dir"], ctx["regulation_dir"], ctx["association_dir"]]
            ),
        )
    ]

    for cell in ctx["cell_types"]:
        input_h5ad = ctx["h5ad_dir"] / f"{cell}.h5ad"
        filtered_h5ad = ctx["filtered_h5ad_dir"] / f"{cell}.h5ad"
        cell_cnmf_dir = ctx["cnmf_dir"] / cell

        filter_command = _command_with_runner(
            ctx["perturbseq_runner"],
            "paper-pipeline-perturbseq-filter-cnmf-cells",
            "--input-h5ad",
            str(input_h5ad),
            "--output-h5ad",
            str(filtered_h5ad),
            "--mito-column",
            ctx["mito_column"],
            "--mito-threshold",
            str(ctx["mito_threshold"]),
            "--min-genes-per-cell",
            str(ctx["min_genes"]),
            "--min-counts-per-cell",
            str(ctx["min_counts"]),
            "--min-cells-per-gene",
            str(ctx["min_cells"]),
        )
        tasks.append(
            Task(
                name=f"Filter cells for Perturbseq cNMF essential {cell}",
                preview=preview_command(filter_command, cwd=config.project_root),
                command=filter_command,
                cwd=config.project_root,
            )
        )

        prepare_command = _command_with_runner(
            ctx["perturbseq_runner"],
            "cnmf",
            "prepare",
            "--output-dir",
            str(cell_cnmf_dir),
            "--name",
            "test1",
            "-c",
            str(filtered_h5ad),
            "-k",
            str(ctx["k"]),
            "--n-iter",
            str(ctx["n_iter"]),
            "--seed",
            str(ctx["seed"]),
            "--numgenes",
            str(ctx["numgenes"]),
            "--total-workers",
            str(ctx["total_workers"]),
        )
        tasks.append(
            Task(
                name=f"Prepare cNMF essential run for {cell}",
                preview=preview_command(prepare_command, cwd=config.project_root),
                command=prepare_command,
                cwd=config.project_root,
            )
        )

    return tasks


def _build_essential_factorize_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = []
    for cell in ctx["cell_types"]:
        cell_cnmf_dir = ctx["cnmf_dir"] / cell
        for worker_index in ctx["worker_indices"]:
            factorize_command = _command_with_runner(
                ctx["perturbseq_runner"],
                "cnmf",
                "factorize",
                "--output-dir",
                str(cell_cnmf_dir),
                "--name",
                "test1",
                "--worker-index",
                str(worker_index),
                "--total-workers",
                str(ctx["total_workers"]),
            )
            tasks.append(
                Task(
                    name=f"Factorize cNMF essential {cell} worker {worker_index + 1}",
                    preview=preview_command(factorize_command, cwd=config.project_root),
                    command=factorize_command,
                    cwd=config.project_root,
                )
            )
    return tasks


def _build_essential_postbase_base_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = []
    for cell in ctx["cell_types"]:
        cell_cnmf_dir = ctx["cnmf_dir"] / cell
        cell_reg_dir = ctx["regulation_dir"] / cell
        cell_assoc_dir = ctx["association_dir"] / cell / "ProgramLevel"

        tasks.append(
            Task(
                name=f"Create cNMF essential downstream directories for {cell}",
                preview=f"Create output directories: {cell_reg_dir}, {cell_assoc_dir}",
                action=lambda reg_dir=cell_reg_dir, assoc_dir=cell_assoc_dir: _ensure_directories(
                    [reg_dir, assoc_dir]
                ),
            )
        )

        combine_command = _command_with_runner(
            ctx["perturbseq_runner"],
            "cnmf",
            "combine",
            "--output-dir",
            str(cell_cnmf_dir),
            "--name",
            "test1",
        )
        consensus_command = _command_with_runner(
            ctx["perturbseq_runner"],
            "cnmf",
            "consensus",
            "--output-dir",
            str(cell_cnmf_dir),
            "--name",
            "test1",
            "--components",
            str(ctx["k"]),
            "--local-density-threshold",
            str(ctx["density_threshold"]),
            "--show-clustering",
        )
        tasks.extend(
            [
                Task(
                    name=f"Combine cNMF essential workers for {cell}",
                    preview=preview_command(combine_command, cwd=config.project_root),
                    command=combine_command,
                    cwd=config.project_root,
                ),
                Task(
                    name=f"Consensus cNMF essential programs for {cell}",
                    preview=preview_command(consensus_command, cwd=config.project_root),
                    command=consensus_command,
                    cwd=config.project_root,
                ),
            ]
        )

    return tasks


def _build_essential_postbase_regulation_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = []
    for cell in ctx["cell_types"]:
        metadata_path = ctx["metadata_dir"] / f"{cell}_metadata.csv"
        cell_cnmf_dir = ctx["cnmf_dir"] / cell
        usages_path = cell_cnmf_dir / "test1" / f"test1.usages.k_{ctx['k']}.dt_{ctx['threshold']}.consensus.txt"
        cell_reg_dir = ctx["regulation_dir"] / cell

        tasks.append(
            Task(
                name=f"Create cNMF essential regulation directory for {cell}",
                preview=f"Create output directory: {cell_reg_dir}",
                action=lambda reg_dir=cell_reg_dir: _ensure_directories([reg_dir]),
            )
        )

        for program_index in ctx["program_indices"]:
            output_path = cell_reg_dir / f"K{ctx['k']}_program{program_index}_perturb_effects.txt"
            reg_command = [
                *ctx["rscript_runner"],
                str(ctx["reg_script"]),
                str(metadata_path),
                str(usages_path),
                str(program_index),
                str(output_path),
            ]
            tasks.append(
                Task(
                    name=f"Compute cNMF essential regulatory effects for {cell} program {program_index}",
                    preview=preview_command(reg_command, cwd=config.project_root),
                    command=reg_command,
                    cwd=config.project_root,
                )
            )

    return tasks


def _build_essential_postbase_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    return _build_essential_postbase_base_tasks(config, ctx) + _build_essential_postbase_regulation_tasks(
        config, ctx
    )


def _build_essential_association_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = []
    for cell in ctx["cell_types"]:
        cell_cnmf_dir = ctx["cnmf_dir"] / cell
        gep_path = cell_cnmf_dir / "test1" / f"test1.gene_spectra_score.k_{ctx['k']}.dt_{ctx['threshold']}.txt"
        cell_reg_dir = ctx["regulation_dir"] / cell
        cell_assoc_dir = ctx["association_dir"] / cell / "ProgramLevel"

        tasks.append(
            Task(
                name=f"Create cNMF essential association directory for {cell}",
                preview=f"Create output directory: {cell_assoc_dir}",
                action=lambda assoc_dir=cell_assoc_dir: _ensure_directories([assoc_dir]),
            )
        )

        for trait_file in ctx["trait_files"]:
            lof_path = ctx["lof_dir"] / trait_file
            out_reg = cell_assoc_dir / f"regulators_enrichment_K{ctx['k']}_{trait_file}"
            out_prog = cell_assoc_dir / f"programs_enrichment_K{ctx['k']}_{trait_file}"
            burden_command = [
                *ctx["rscript_runner"],
                str(ctx["burden_script"]),
                str(gep_path),
                str(ctx["gene_map"]),
                str(cell_reg_dir),
                str(lof_path),
                str(ctx["shet_path"]),
                str(ctx["k"]),
                str(out_reg),
                str(out_prog),
                str(ctx["random_iterations"]),
                str(ctx["top_n"]),
            ]
            tasks.append(
                Task(
                    name=f"Run cNMF essential program-trait association for {cell} and {trait_file}",
                    preview=preview_command(burden_command, cwd=config.project_root),
                    command=burden_command,
                    cwd=config.project_root,
                )
            )

    return tasks


def _build_essential_post_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    return _build_essential_postbase_tasks(config, ctx) + _build_essential_association_tasks(config, ctx)


def _load_essential_kselect_context(config: LoadedConfig) -> dict:
    perturbseq = config.workflow("perturbseq")
    workflow = perturbseq.get("cnmf_essential_kselect", {})
    if not workflow:
        raise ConfigError("workflows.perturbseq.cnmf_essential_kselect is required")

    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})

    if "h5ad" not in inputs:
        raise ConfigError("workflows.perturbseq.cnmf_essential_kselect.inputs.h5ad is required")

    perturbseq_runner = config.executable_command("perturbseq_runner", [])
    h5ad_path = config.resolve_path(inputs["h5ad"])
    cnmf_root = config.resolve_path_or_artifact(
        outputs.get("cnmf_root"),
        "perturbseq",
        "cnmf_essential_kselect",
        "cNMF",
    )

    assert h5ad_path

    ks = [int(value) for value in parameters.get("ks", [30, 40, 50, 60, 70, 80, 90, 120])]
    if not ks:
        raise ConfigError("workflows.perturbseq.cnmf_essential_kselect.parameters.ks must not be empty")
    aggregate_name = str(parameters.get("aggregate_name", "cNMF_all"))
    per_k_name_template = str(parameters.get("per_k_name_template", "cNMF_K{k}"))
    total_workers = int(parameters.get("total_workers", 20))
    worker_indices = _parse_worker_indices(
        parameters,
        total_workers,
        "workflows.perturbseq.cnmf_essential_kselect.parameters.worker_indices",
    )
    n_iter = int(parameters.get("n_iter", 100))
    seed = int(parameters.get("seed", 14))
    numgenes = int(parameters.get("numgenes", 2000))

    return {
        "workflow": workflow,
        "perturbseq_runner": perturbseq_runner,
        "h5ad_path": h5ad_path,
        "cnmf_root": cnmf_root,
        "ks": ks,
        "aggregate_name": aggregate_name,
        "per_k_name_template": per_k_name_template,
        "total_workers": total_workers,
        "worker_indices": worker_indices,
        "n_iter": n_iter,
        "seed": seed,
        "numgenes": numgenes,
    }


def _build_essential_kselect_prepare_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = [
        Task(
            name="Create Perturbseq cNMF essential K-selection output directory",
            preview=f"Create output directory: {ctx['cnmf_root']}",
            action=lambda: _ensure_directories([ctx["cnmf_root"]]),
        )
    ]

    for k in ctx["ks"]:
        run_name = ctx["per_k_name_template"].format(k=k)
        prepare_command = _command_with_runner(
            ctx["perturbseq_runner"],
            "cnmf",
            "prepare",
            "--output-dir",
            str(ctx["cnmf_root"]),
            "--name",
            run_name,
            "-c",
            str(ctx["h5ad_path"]),
            "-k",
            str(k),
            "--n-iter",
            str(ctx["n_iter"]),
            "--seed",
            str(ctx["seed"]),
            "--numgenes",
            str(ctx["numgenes"]),
            "--total-workers",
            str(ctx["total_workers"]),
        )
        tasks.append(
            Task(
                name=f"Prepare cNMF essential K-selection run {run_name}",
                preview=preview_command(prepare_command, cwd=config.project_root),
                command=prepare_command,
                cwd=config.project_root,
            )
        )

    aggregate_prepare_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "cnmf",
        "prepare",
        "--output-dir",
        str(ctx["cnmf_root"]),
        "--name",
        ctx["aggregate_name"],
        "-c",
        str(ctx["h5ad_path"]),
        "-k",
        *[str(value) for value in ctx["ks"]],
        "--n-iter",
        str(ctx["n_iter"]),
        "--seed",
        str(ctx["seed"]),
        "--numgenes",
        str(ctx["numgenes"]),
        "--total-workers",
        str(ctx["total_workers"]),
    )
    tasks.append(
        Task(
            name=f"Prepare aggregate cNMF essential K-selection run {ctx['aggregate_name']}",
            preview=preview_command(aggregate_prepare_command, cwd=config.project_root),
            command=aggregate_prepare_command,
            cwd=config.project_root,
        )
    )
    return tasks


def _build_essential_kselect_factorize_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = []
    for k in ctx["ks"]:
        run_name = ctx["per_k_name_template"].format(k=k)
        for worker_index in ctx["worker_indices"]:
            factorize_command = _command_with_runner(
                ctx["perturbseq_runner"],
                "cnmf",
                "factorize",
                "--output-dir",
                str(ctx["cnmf_root"]),
                "--name",
                run_name,
                "--worker-index",
                str(worker_index),
                "--total-workers",
                str(ctx["total_workers"]),
            )
            tasks.append(
                Task(
                    name=f"Factorize cNMF essential K-selection {run_name} worker {worker_index + 1}",
                    preview=preview_command(factorize_command, cwd=config.project_root),
                    command=factorize_command,
                    cwd=config.project_root,
                )
            )
    return tasks


def _build_essential_kselect_postbase_base_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    per_k_dirs: list[Path] = []
    per_k_names: list[str] = []
    for k in ctx["ks"]:
        run_name = ctx["per_k_name_template"].format(k=k)
        per_k_dirs.append(ctx["cnmf_root"] / run_name)
        per_k_names.append(run_name)

    merge_parts: list[str] = []
    for path, name in zip(per_k_dirs, per_k_names, strict=True):
        merge_parts.extend(["--source-dir", str(path), "--source-name", name])

    merge_tmp_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "paper-pipeline-perturbseq-merge-cnmf-tmp",
        *merge_parts,
        "--target-dir",
        str(ctx["cnmf_root"] / ctx["aggregate_name"]),
        "--target-name",
        ctx["aggregate_name"],
    )
    combine_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "cnmf",
        "combine",
        "--output-dir",
        str(ctx["cnmf_root"]),
        "--name",
        ctx["aggregate_name"],
    )
    k_selection_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "cnmf",
        "k_selection_plot",
        "--output-dir",
        str(ctx["cnmf_root"]),
        "--name",
        ctx["aggregate_name"],
    )

    return [
        Task(
            name=f"Merge per-K cNMF essential temporary outputs into {ctx['aggregate_name']}",
            preview=preview_command(merge_tmp_command, cwd=config.project_root),
            command=merge_tmp_command,
            cwd=config.project_root,
        ),
        Task(
            name=f"Combine aggregate cNMF essential workers for {ctx['aggregate_name']}",
            preview=preview_command(combine_command, cwd=config.project_root),
            command=combine_command,
            cwd=config.project_root,
        ),
        Task(
            name=f"Generate cNMF essential K-selection plot for {ctx['aggregate_name']}",
            preview=preview_command(k_selection_command, cwd=config.project_root),
            command=k_selection_command,
            cwd=config.project_root,
        ),
    ]


def build_perturbseq_cnmf_essential_kselect_prepare_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_kselect_context(config)
    return _build_essential_kselect_prepare_tasks(config, ctx)


def build_perturbseq_cnmf_essential_kselect_factorize_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_kselect_context(config)
    return _build_essential_kselect_factorize_tasks(config, ctx)


def build_perturbseq_cnmf_essential_kselect_postbase_base_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_kselect_context(config)
    return _build_essential_kselect_postbase_base_tasks(config, ctx)


def build_perturbseq_cnmf_essential_prepare_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    return _build_essential_prepare_tasks(config, ctx)


def build_perturbseq_cnmf_essential_factorize_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    return _build_essential_factorize_tasks(config, ctx)


def build_perturbseq_cnmf_essential_postbase_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    _require_full_program_set(
        ctx["program_indices"],
        ctx["k"],
        workflow_name="perturbseq-cnmf-essential-postbase",
        regulation_workflow_name="perturbseq-cnmf-essential-postbase-regulation",
    )
    return _build_essential_postbase_tasks(config, ctx)


def build_perturbseq_cnmf_essential_postbase_base_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    return _build_essential_postbase_base_tasks(config, ctx)


def build_perturbseq_cnmf_essential_postbase_regulation_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    return _build_essential_postbase_regulation_tasks(config, ctx)


def build_perturbseq_cnmf_essential_association_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    _require_full_program_set(
        ctx["program_indices"],
        ctx["k"],
        workflow_name="perturbseq-cnmf-essential-association",
        regulation_workflow_name="perturbseq-cnmf-essential-postbase-regulation",
    )
    return _build_essential_association_tasks(config, ctx)


def build_perturbseq_cnmf_essential_post_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    _require_full_program_set(
        ctx["program_indices"],
        ctx["k"],
        workflow_name="perturbseq-cnmf-essential-post",
        regulation_workflow_name="perturbseq-cnmf-essential-postbase-regulation",
    )
    return _build_essential_post_tasks(config, ctx)


def build_perturbseq_cnmf_essential_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_essential_context(config)
    _require_full_worker_set(
        ctx["worker_indices"],
        ctx["total_workers"],
        workflow_name="perturbseq-cnmf-essential",
        factorize_workflow_name="perturbseq-cnmf-essential-factorize",
        post_workflow_name="perturbseq-cnmf-essential-post",
    )
    _require_full_program_set(
        ctx["program_indices"],
        ctx["k"],
        workflow_name="perturbseq-cnmf-essential",
        regulation_workflow_name="perturbseq-cnmf-essential-postbase-regulation",
    )
    return (
        _build_essential_prepare_tasks(config, ctx)
        + _build_essential_factorize_tasks(config, ctx)
        + _build_essential_post_tasks(config, ctx)
    )


def _load_genomewide_context(config: LoadedConfig) -> dict:
    perturbseq = config.workflow("perturbseq")
    workflow = perturbseq.get("cnmf_genomewide", {})
    if not workflow:
        raise ConfigError("workflows.perturbseq.cnmf_genomewide is required")

    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})

    required_input_keys = ["h5ad", "metadata", "gene_map", "shet_path"]
    for key in required_input_keys:
        if key not in inputs:
            raise ConfigError(f"workflows.perturbseq.cnmf_genomewide.inputs.{key} is required")

    perturbseq_runner = config.executable_command("perturbseq_runner", [])
    rscript_runner = config.executable_command("rscript", "Rscript")

    h5ad_path = config.resolve_path(inputs["h5ad"])
    metadata_path = config.resolve_path(inputs["metadata"])
    gene_map = config.resolve_path(inputs["gene_map"])
    lof_dir = config.resolve_path_or_artifact(inputs.get("lof_dir"), "genebayes", "posterior")
    shet_path = config.resolve_path(inputs["shet_path"])

    cnmf_root = config.resolve_path_or_artifact(
        outputs.get("cnmf_root"),
        "perturbseq",
        "cnmf_genomewide",
        "cNMF",
    )
    regulation_dir = config.resolve_path_or_artifact(
        outputs.get("regulation_dir"),
        "perturbseq",
        "cnmf_genomewide",
        "cNMF_regulation",
        "K562GW",
    )
    association_dir = config.resolve_path_or_artifact(
        outputs.get("association_dir"),
        "perturbseq",
        "cnmf_genomewide",
        "trait_association",
        "K562GW",
        "ProgramLevel",
    )

    assert h5ad_path and metadata_path and gene_map and lof_dir and shet_path

    ks = [int(value) for value in parameters.get("ks", [30, 60, 90, 120])]
    if not ks:
        raise ConfigError("workflows.perturbseq.cnmf_genomewide.parameters.ks must not be empty")
    aggregate_name = str(parameters.get("aggregate_name", "cNMF_all"))
    per_k_name_template = str(parameters.get("per_k_name_template", "cNMF_K{k}"))
    total_workers = int(parameters.get("total_workers", 200))
    worker_indices = _parse_worker_indices(
        parameters,
        total_workers,
        "workflows.perturbseq.cnmf_genomewide.parameters.worker_indices",
    )
    n_iter = int(parameters.get("n_iter", 100))
    seed = int(parameters.get("seed", 14))
    numgenes = int(parameters.get("numgenes", 2000))
    consensus_ks = [int(value) for value in parameters.get("consensus_ks", ks)]
    consensus_density_threshold = float(parameters.get("consensus_density_threshold", 0.5))
    association_k = int(parameters.get("association_k", 60))
    random_iterations = int(parameters.get("random_iterations", 100000))
    top_n = int(parameters.get("top_n_program_genes", 100))
    trait_files = _trait_files(workflow, lof_dir)
    threshold = _threshold_label(consensus_density_threshold)
    program_indices = _parse_program_indices(
        parameters,
        association_k,
        "workflows.perturbseq.cnmf_genomewide.parameters.program_indices",
    )

    return {
        "workflow": workflow,
        "perturbseq_runner": perturbseq_runner,
        "rscript_runner": rscript_runner,
        "h5ad_path": h5ad_path,
        "metadata_path": metadata_path,
        "gene_map": gene_map,
        "lof_dir": lof_dir,
        "shet_path": shet_path,
        "cnmf_root": cnmf_root,
        "regulation_dir": regulation_dir,
        "association_dir": association_dir,
        "ks": ks,
        "aggregate_name": aggregate_name,
        "per_k_name_template": per_k_name_template,
        "total_workers": total_workers,
        "worker_indices": worker_indices,
        "n_iter": n_iter,
        "seed": seed,
        "numgenes": numgenes,
        "consensus_ks": consensus_ks,
        "consensus_density_threshold": consensus_density_threshold,
        "association_k": association_k,
        "random_iterations": random_iterations,
        "top_n": top_n,
        "trait_files": trait_files,
        "program_indices": program_indices,
        "threshold": threshold,
        "reg_script": r_script_path("cnmf_regulatory_effects.R"),
        "burden_script": r_script_path("cnmf_burden_program_regulators.R"),
    }


def _build_genomewide_prepare_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = [
        Task(
            name="Create Perturbseq cNMF genomewide output directories",
            preview=f"Create output directories: {ctx['cnmf_root']}, {ctx['regulation_dir']}, {ctx['association_dir']}",
            action=lambda: _ensure_directories([ctx["cnmf_root"], ctx["regulation_dir"], ctx["association_dir"]]),
        )
    ]

    for k in ctx["ks"]:
        run_name = ctx["per_k_name_template"].format(k=k)
        prepare_command = _command_with_runner(
            ctx["perturbseq_runner"],
            "cnmf",
            "prepare",
            "--output-dir",
            str(ctx["cnmf_root"]),
            "--name",
            run_name,
            "-c",
            str(ctx["h5ad_path"]),
            "-k",
            str(k),
            "--n-iter",
            str(ctx["n_iter"]),
            "--seed",
            str(ctx["seed"]),
            "--numgenes",
            str(ctx["numgenes"]),
            "--total-workers",
            str(ctx["total_workers"]),
        )
        tasks.append(
            Task(
                name=f"Prepare cNMF genomewide run {run_name}",
                preview=preview_command(prepare_command, cwd=config.project_root),
                command=prepare_command,
                cwd=config.project_root,
            )
        )

    aggregate_prepare_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "cnmf",
        "prepare",
        "--output-dir",
        str(ctx["cnmf_root"]),
        "--name",
        ctx["aggregate_name"],
        "-c",
        str(ctx["h5ad_path"]),
        "-k",
        *[str(value) for value in ctx["ks"]],
        "--n-iter",
        str(ctx["n_iter"]),
        "--seed",
        str(ctx["seed"]),
        "--numgenes",
        str(ctx["numgenes"]),
        "--total-workers",
        str(ctx["total_workers"]),
    )
    tasks.append(
        Task(
            name=f"Prepare aggregate cNMF genomewide run {ctx['aggregate_name']}",
            preview=preview_command(aggregate_prepare_command, cwd=config.project_root),
            command=aggregate_prepare_command,
            cwd=config.project_root,
        )
    )
    return tasks


def _build_genomewide_factorize_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = []
    for k in ctx["ks"]:
        run_name = ctx["per_k_name_template"].format(k=k)
        for worker_index in ctx["worker_indices"]:
            factorize_command = _command_with_runner(
                ctx["perturbseq_runner"],
                "cnmf",
                "factorize",
                "--output-dir",
                str(ctx["cnmf_root"]),
                "--name",
                run_name,
                "--worker-index",
                str(worker_index),
                "--total-workers",
                str(ctx["total_workers"]),
            )
            tasks.append(
                Task(
                    name=f"Factorize cNMF genomewide {run_name} worker {worker_index + 1}",
                    preview=preview_command(factorize_command, cwd=config.project_root),
                    command=factorize_command,
                    cwd=config.project_root,
                )
            )
    return tasks


def _build_genomewide_postbase_base_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = []

    per_k_dirs: list[Path] = []
    per_k_names: list[str] = []
    for k in ctx["ks"]:
        run_name = ctx["per_k_name_template"].format(k=k)
        per_k_dirs.append(ctx["cnmf_root"] / run_name)
        per_k_names.append(run_name)

    merge_parts: list[str] = []
    for path, name in zip(per_k_dirs, per_k_names, strict=True):
        merge_parts.extend(["--source-dir", str(path), "--source-name", name])
    merge_tmp_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "paper-pipeline-perturbseq-merge-cnmf-tmp",
        *merge_parts,
        "--target-dir",
        str(ctx["cnmf_root"] / ctx["aggregate_name"]),
        "--target-name",
        ctx["aggregate_name"],
    )
    tasks.append(
        Task(
            name=f"Merge per-K cNMF temporary outputs into {ctx['aggregate_name']}",
            preview=preview_command(merge_tmp_command, cwd=config.project_root),
            command=merge_tmp_command,
            cwd=config.project_root,
        )
    )

    combine_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "cnmf",
        "combine",
        "--output-dir",
        str(ctx["cnmf_root"]),
        "--name",
        ctx["aggregate_name"],
    )
    k_selection_command = _command_with_runner(
        ctx["perturbseq_runner"],
        "cnmf",
        "k_selection_plot",
        "--output-dir",
        str(ctx["cnmf_root"]),
        "--name",
        ctx["aggregate_name"],
    )
    tasks.extend(
        [
            Task(
                name=f"Combine aggregate cNMF genomewide workers for {ctx['aggregate_name']}",
                preview=preview_command(combine_command, cwd=config.project_root),
                command=combine_command,
                cwd=config.project_root,
            ),
            Task(
                name=f"Generate cNMF K selection plot for {ctx['aggregate_name']}",
                preview=preview_command(k_selection_command, cwd=config.project_root),
                command=k_selection_command,
                cwd=config.project_root,
            ),
        ]
    )

    for k in ctx["consensus_ks"]:
        consensus_command = _command_with_runner(
            ctx["perturbseq_runner"],
            "cnmf",
            "consensus",
            "--output-dir",
            str(ctx["cnmf_root"]),
            "--name",
            ctx["aggregate_name"],
            "--components",
            str(k),
            "--local-density-threshold",
            str(ctx["consensus_density_threshold"]),
            "--show-clustering",
        )
        tasks.append(
            Task(
                name=f"Consensus aggregate cNMF genomewide programs for K={k}",
                preview=preview_command(consensus_command, cwd=config.project_root),
                command=consensus_command,
                cwd=config.project_root,
            )
        )

    return tasks


def _build_genomewide_postbase_regulation_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = [
        Task(
            name="Create cNMF genomewide regulation directory",
            preview=f"Create output directory: {ctx['regulation_dir']}",
            action=lambda: _ensure_directories([ctx["regulation_dir"]]),
        )
    ]

    usages_path = (
        ctx["cnmf_root"]
        / ctx["aggregate_name"]
        / f"{ctx['aggregate_name']}.usages.k_{ctx['association_k']}.dt_{ctx['threshold']}.consensus.txt"
    )

    for program_index in ctx["program_indices"]:
        output_path = ctx["regulation_dir"] / f"K{ctx['association_k']}_program{program_index}_perturb_effects.txt"
        reg_command = [
            *ctx["rscript_runner"],
            str(ctx["reg_script"]),
            str(ctx["metadata_path"]),
            str(usages_path),
            str(program_index),
            str(output_path),
        ]
        tasks.append(
            Task(
                name=f"Compute cNMF genomewide regulatory effects for program {program_index}",
                preview=preview_command(reg_command, cwd=config.project_root),
                command=reg_command,
                cwd=config.project_root,
            )
        )

    return tasks


def _build_genomewide_postbase_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    return _build_genomewide_postbase_base_tasks(config, ctx) + _build_genomewide_postbase_regulation_tasks(
        config, ctx
    )


def _build_genomewide_association_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    tasks: list[Task] = [
        Task(
            name="Create cNMF genomewide association directory",
            preview=f"Create output directory: {ctx['association_dir']}",
            action=lambda: _ensure_directories([ctx["association_dir"]]),
        )
    ]

    gep_path = (
        ctx["cnmf_root"]
        / ctx["aggregate_name"]
        / f"{ctx['aggregate_name']}.gene_spectra_score.k_{ctx['association_k']}.dt_{ctx['threshold']}.txt"
    )

    for trait_file in ctx["trait_files"]:
        lof_path = ctx["lof_dir"] / trait_file
        out_reg = ctx["association_dir"] / f"regulators_enrichment_K{ctx['association_k']}_{trait_file}"
        out_prog = ctx["association_dir"] / f"programs_enrichment_K{ctx['association_k']}_{trait_file}"
        burden_command = [
            *ctx["rscript_runner"],
            str(ctx["burden_script"]),
            str(gep_path),
            str(ctx["gene_map"]),
            str(ctx["regulation_dir"]),
            str(lof_path),
            str(ctx["shet_path"]),
            str(ctx["association_k"]),
            str(out_reg),
            str(out_prog),
            str(ctx["random_iterations"]),
            str(ctx["top_n"]),
        ]
        tasks.append(
            Task(
                name=f"Run cNMF genomewide program-trait association for {trait_file}",
                preview=preview_command(burden_command, cwd=config.project_root),
                command=burden_command,
                cwd=config.project_root,
            )
        )

    return tasks


def _build_genomewide_post_tasks(config: LoadedConfig, ctx: dict) -> list[Task]:
    return _build_genomewide_postbase_tasks(config, ctx) + _build_genomewide_association_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_prepare_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    return _build_genomewide_prepare_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_factorize_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    return _build_genomewide_factorize_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_postbase_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    _require_full_program_set(
        ctx["program_indices"],
        ctx["association_k"],
        workflow_name="perturbseq-cnmf-genomewide-postbase",
        regulation_workflow_name="perturbseq-cnmf-genomewide-postbase-regulation",
    )
    return _build_genomewide_postbase_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_postbase_base_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    return _build_genomewide_postbase_base_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_postbase_regulation_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    return _build_genomewide_postbase_regulation_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_association_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    _require_full_program_set(
        ctx["program_indices"],
        ctx["association_k"],
        workflow_name="perturbseq-cnmf-genomewide-association",
        regulation_workflow_name="perturbseq-cnmf-genomewide-postbase-regulation",
    )
    return _build_genomewide_association_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_post_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    _require_full_program_set(
        ctx["program_indices"],
        ctx["association_k"],
        workflow_name="perturbseq-cnmf-genomewide-post",
        regulation_workflow_name="perturbseq-cnmf-genomewide-postbase-regulation",
    )
    return _build_genomewide_post_tasks(config, ctx)


def build_perturbseq_cnmf_genomewide_tasks(config: LoadedConfig) -> list[Task]:
    ctx = _load_genomewide_context(config)
    _require_full_worker_set(
        ctx["worker_indices"],
        ctx["total_workers"],
        workflow_name="perturbseq-cnmf-genomewide",
        factorize_workflow_name="perturbseq-cnmf-genomewide-factorize",
        post_workflow_name="perturbseq-cnmf-genomewide-post",
    )
    _require_full_program_set(
        ctx["program_indices"],
        ctx["association_k"],
        workflow_name="perturbseq-cnmf-genomewide",
        regulation_workflow_name="perturbseq-cnmf-genomewide-postbase-regulation",
    )
    return (
        _build_genomewide_prepare_tasks(config, ctx)
        + _build_genomewide_factorize_tasks(config, ctx)
        + _build_genomewide_post_tasks(config, ctx)
    )
