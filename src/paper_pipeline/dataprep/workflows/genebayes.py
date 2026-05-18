from __future__ import annotations

from pathlib import Path

from ..config import ConfigError, LoadedConfig
from ..runner import preview_command
from ..tasks import Task


def _trait_labels(config: dict) -> list[str]:
    traits = config.get("traits", [])
    if not traits:
        raise ConfigError("workflows.genebayes.traits must not be empty")

    trait_prefix = str(config.get("trait_prefix", ""))
    labels: list[str] = []
    for trait in traits:
        label = str(trait)
        if trait_prefix and not label.startswith(trait_prefix):
            label = f"{trait_prefix}{label}"
        labels.append(label)
    return labels


def _ensure_directories(paths: list[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def _ensure_files_exist(paths: list[Path], *, label: str) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        joined = ", ".join(missing)
        raise ConfigError(f"{label} are required but missing: {joined}")


def _resolve_legacy_runner(config: LoadedConfig, workflow: dict) -> list[str]:
    runner_name = str(workflow.get("runner", "genebayes_runner"))
    runner = config.executable_command(runner_name, None if runner_name != "genebayes_runner" else [])
    if runner_name != "genebayes_runner" and not runner:
        raise ConfigError(f"executables.{runner_name} is required")
    return runner


def _resolve_phase_runner(
    config: LoadedConfig,
    workflow: dict,
    *,
    phase: str,
    legacy_runner: list[str],
) -> list[str]:
    workflow_key = f"{phase}_runner"
    executable_key = f"genebayes_{phase}_runner"

    if workflow_key in workflow:
        runner_name = str(workflow[workflow_key])
        runner = config.executable_command(runner_name, None if runner_name != executable_key else [])
        if runner_name != executable_key and not runner:
            raise ConfigError(f"executables.{runner_name} is required")
        return runner

    runner = config.executable_command(executable_key, None)
    if runner:
        return runner
    return legacy_runner


def _build_prediction_tasks(
    *,
    config: LoadedConfig,
    trait_labels: list[str],
    traits_dir: Path,
    response_template: str,
    runner: list[str],
    s_het: Path,
    embeddings: Path,
    ntpm: Path,
    geneformer: Path,
    prediction_dir: Path,
    prefer_gpu: bool,
) -> list[Task]:
    tasks: list[Task] = []
    for trait_label in trait_labels:
        trait_file = traits_dir / response_template.format(trait=trait_label)
        signed_command = runner + [
            "paper-pipeline-genebayes-signed",
            trait_label,
            "--s-het",
            str(s_het),
            "--embeddings",
            str(embeddings),
            "--celltype-ntpm",
            str(ntpm),
            "--geneformer-mixture",
            str(geneformer),
            "--trait-file",
            str(trait_file),
            "--output-dir",
            str(prediction_dir),
        ]
        magnitude_command = runner + [
            "paper-pipeline-genebayes-magnitude",
            trait_label,
            "--s-het",
            str(s_het),
            "--embeddings",
            str(embeddings),
            "--celltype-ntpm",
            str(ntpm),
            "--geneformer-mixture",
            str(geneformer),
            "--trait-file",
            str(trait_file),
            "--output-dir",
            str(prediction_dir),
        ]
        if not prefer_gpu:
            signed_command.append("--no-prefer-gpu")
            magnitude_command.append("--no-prefer-gpu")

        tasks.append(
            Task(
                name=f"GeneBayes signed feature prediction for {trait_label}",
                preview=preview_command(signed_command, cwd=config.project_root),
                command=signed_command,
                cwd=config.project_root,
            )
        )
        tasks.append(
            Task(
                name=f"GeneBayes magnitude feature prediction for {trait_label}",
                preview=preview_command(magnitude_command, cwd=config.project_root),
                command=magnitude_command,
                cwd=config.project_root,
            )
        )
    return tasks


def _build_merge_task(
    *,
    config: LoadedConfig,
    trait_labels: list[str],
    runner: list[str],
    prediction_dir: Path,
    feature_dir: Path,
) -> Task:
    merge_command = runner + [
        "paper-pipeline-genebayes-merge",
        "--prediction-dir",
        str(prediction_dir),
        "--feature-dir",
        str(feature_dir),
        "--traits",
        *trait_labels,
    ]
    return Task(
        name="Merge GeneBayes signed and magnitude feature predictions",
        preview=preview_command(merge_command, cwd=config.project_root),
        command=merge_command,
        cwd=config.project_root,
    )


def _build_training_tasks(
    *,
    config: LoadedConfig,
    trait_labels: list[str],
    traits_dir: Path,
    response_template: str,
    runner: list[str],
    feature_dir: Path,
    posterior_dir: Path,
    train_genes: Path,
    val_genes: Path,
    parameters: dict,
) -> list[Task]:
    numeric_parameter_flags = {
        "integration_lb": "--integration_lb",
        "integration_ub": "--integration_ub",
        "batch_size": "--batch_size",
        "n_integration_pts": "--n_integration_pts",
        "lr": "--lr",
        "total_iterations": "--total_iterations",
        "early_stopping_iter": "--early_stopping_iter",
        "n_trees_per_iteration": "--n_trees_per_iteration",
        "reg_alpha": "--reg_alpha",
        "reg_lambda": "--reg_lambda",
        "subsample": "--subsample",
        "min_child_weight": "--min_child_weight",
        "max_depth": "--max_depth",
    }

    tasks: list[Task] = []
    for trait_label in trait_labels:
        trait_file = traits_dir / response_template.format(trait=trait_label)
        out_prefix = posterior_dir / trait_label
        train_command = runner + [
            "paper-pipeline-genebayes-train",
            "--response",
            str(trait_file),
            "--features",
            str(feature_dir / f"features.{trait_label}.tsv"),
            "--train_genes",
            str(train_genes),
            "--val_genes",
            str(val_genes),
            "--out",
            str(out_prefix),
        ]
        for key, flag in numeric_parameter_flags.items():
            value = parameters.get(key)
            if value is not None:
                train_command.extend([flag, str(value)])

        tasks.append(
            Task(
                name=f"Train GeneBayes posterior model for {trait_label}",
                preview=preview_command(train_command, cwd=config.project_root),
                command=train_command,
                cwd=config.project_root,
            )
        )
    return tasks


def build_genebayes_tasks(config: LoadedConfig) -> list[Task]:
    workflow = config.workflow("genebayes")
    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})
    trait_labels = _trait_labels(workflow)

    required_input_keys = [
        "s_het",
        "reduced_embeddings",
        "celltype_ntpm",
        "geneformer_mixture",
        "traits_dir",
        "train_gene_list",
        "val_gene_list",
    ]
    for key in required_input_keys:
        if key not in inputs:
            raise ConfigError(f"workflows.genebayes.inputs.{key} is required")

    s_het = config.resolve_path(inputs["s_het"])
    embeddings = config.resolve_path(inputs["reduced_embeddings"])
    ntpm = config.resolve_path(inputs["celltype_ntpm"])
    geneformer = config.resolve_path(inputs["geneformer_mixture"])
    traits_dir = config.resolve_path(inputs["traits_dir"])
    train_genes = config.resolve_path(inputs["train_gene_list"])
    val_genes = config.resolve_path(inputs["val_gene_list"])
    prediction_dir = config.resolve_path_or_artifact(outputs.get("prediction_dir"), "genebayes", "predictions")
    feature_dir = config.resolve_path_or_artifact(outputs.get("feature_dir"), "genebayes", "features")
    posterior_dir = config.resolve_path_or_artifact(outputs.get("posterior_dir"), "genebayes", "posterior")

    if not (s_het and embeddings and ntpm and geneformer):
        raise ConfigError("workflows.genebayes.inputs: s_het, reduced_embeddings, celltype_ntpm, geneformer_mixture are all required")
    if not (traits_dir and train_genes and val_genes):
        raise ConfigError("workflows.genebayes.inputs: traits_dir, train_gene_list, val_gene_list are all required")

    response_template = str(parameters.get("response_file_template", "{trait}.summary_statistics.csv"))
    prefer_gpu = bool(parameters.get("prefer_gpu", True))
    genebayes_runner = _resolve_legacy_runner(config, workflow)

    tasks: list[Task] = [
        Task(
            name="Create GeneBayes output directories",
            preview=f"Create output directories: {prediction_dir}, {feature_dir}, {posterior_dir}",
            action=lambda: _ensure_directories([prediction_dir, feature_dir, posterior_dir]),
        )
    ]

    tasks.extend(
        _build_prediction_tasks(
            config=config,
            trait_labels=trait_labels,
            traits_dir=traits_dir,
            response_template=response_template,
            runner=genebayes_runner,
            s_het=s_het,
            embeddings=embeddings,
            ntpm=ntpm,
            geneformer=geneformer,
            prediction_dir=prediction_dir,
            prefer_gpu=prefer_gpu,
        )
    )
    tasks.append(
        _build_merge_task(
            config=config,
            trait_labels=trait_labels,
            runner=genebayes_runner,
            prediction_dir=prediction_dir,
            feature_dir=feature_dir,
        )
    )
    tasks.extend(
        _build_training_tasks(
            config=config,
            trait_labels=trait_labels,
            traits_dir=traits_dir,
            response_template=response_template,
            runner=genebayes_runner,
            feature_dir=feature_dir,
            posterior_dir=posterior_dir,
            train_genes=train_genes,
            val_genes=val_genes,
            parameters=parameters,
        )
    )
    return tasks


def build_genebayes_cpu_tasks(config: LoadedConfig) -> list[Task]:
    workflow = config.workflow("genebayes")
    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})
    trait_labels = _trait_labels(workflow)

    required_input_keys = [
        "s_het",
        "reduced_embeddings",
        "celltype_ntpm",
        "geneformer_mixture",
        "traits_dir",
    ]
    for key in required_input_keys:
        if key not in inputs:
            raise ConfigError(f"workflows.genebayes.inputs.{key} is required")

    s_het = config.resolve_path(inputs["s_het"])
    embeddings = config.resolve_path(inputs["reduced_embeddings"])
    ntpm = config.resolve_path(inputs["celltype_ntpm"])
    geneformer = config.resolve_path(inputs["geneformer_mixture"])
    traits_dir = config.resolve_path(inputs["traits_dir"])
    prediction_dir = config.resolve_path_or_artifact(outputs.get("prediction_dir"), "genebayes", "predictions")
    feature_dir = config.resolve_path_or_artifact(outputs.get("feature_dir"), "genebayes", "features")
    posterior_dir = config.resolve_path_or_artifact(outputs.get("posterior_dir"), "genebayes", "posterior")

    if not (s_het and embeddings and ntpm and geneformer and traits_dir):
        raise ConfigError("workflows.genebayes.inputs: s_het, reduced_embeddings, celltype_ntpm, geneformer_mixture, traits_dir are all required")

    response_template = str(parameters.get("response_file_template", "{trait}.summary_statistics.csv"))
    legacy_runner = _resolve_legacy_runner(config, workflow)
    cpu_runner = _resolve_phase_runner(config, workflow, phase="cpu", legacy_runner=legacy_runner)

    tasks: list[Task] = [
        Task(
            name="Create GeneBayes output directories",
            preview=f"Create output directories: {prediction_dir}, {feature_dir}, {posterior_dir}",
            action=lambda: _ensure_directories([prediction_dir, feature_dir, posterior_dir]),
        )
    ]
    tasks.extend(
        _build_prediction_tasks(
            config=config,
            trait_labels=trait_labels,
            traits_dir=traits_dir,
            response_template=response_template,
            runner=cpu_runner,
            s_het=s_het,
            embeddings=embeddings,
            ntpm=ntpm,
            geneformer=geneformer,
            prediction_dir=prediction_dir,
            prefer_gpu=False,
        )
    )
    tasks.append(
        _build_merge_task(
            config=config,
            trait_labels=trait_labels,
            runner=cpu_runner,
            prediction_dir=prediction_dir,
            feature_dir=feature_dir,
        )
    )
    return tasks


def build_genebayes_gpu_tasks(config: LoadedConfig) -> list[Task]:
    workflow = config.workflow("genebayes")
    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})
    trait_labels = _trait_labels(workflow)

    required_input_keys = [
        "traits_dir",
        "train_gene_list",
        "val_gene_list",
    ]
    for key in required_input_keys:
        if key not in inputs:
            raise ConfigError(f"workflows.genebayes.inputs.{key} is required")

    traits_dir = config.resolve_path(inputs["traits_dir"])
    train_genes = config.resolve_path(inputs["train_gene_list"])
    val_genes = config.resolve_path(inputs["val_gene_list"])
    feature_dir = config.resolve_path_or_artifact(outputs.get("feature_dir"), "genebayes", "features")
    posterior_dir = config.resolve_path_or_artifact(outputs.get("posterior_dir"), "genebayes", "posterior")
    prediction_dir = config.resolve_path_or_artifact(outputs.get("prediction_dir"), "genebayes", "predictions")

    if not (traits_dir and train_genes and val_genes):
        raise ConfigError("workflows.genebayes.inputs: traits_dir, train_gene_list, val_gene_list are all required")

    response_template = str(parameters.get("response_file_template", "{trait}.summary_statistics.csv"))
    legacy_runner = _resolve_legacy_runner(config, workflow)
    gpu_runner = _resolve_phase_runner(config, workflow, phase="gpu", legacy_runner=legacy_runner)

    tasks: list[Task] = [
        Task(
            name="Create GeneBayes output directories",
            preview=f"Create output directories: {prediction_dir}, {feature_dir}, {posterior_dir}",
            action=lambda: _ensure_directories([prediction_dir, feature_dir, posterior_dir]),
        ),
        Task(
            name="Verify merged GeneBayes feature tables exist",
            preview=(
                "Verify merged feature tables exist: "
                + ", ".join(str(feature_dir / f"features.{trait_label}.tsv") for trait_label in trait_labels)
            ),
            action=lambda: _ensure_files_exist(
                [feature_dir / f"features.{trait_label}.tsv" for trait_label in trait_labels],
                label="Merged GeneBayes feature tables",
            ),
        ),
    ]
    tasks.extend(
        _build_training_tasks(
            config=config,
            trait_labels=trait_labels,
            traits_dir=traits_dir,
            response_template=response_template,
            runner=gpu_runner,
            feature_dir=feature_dir,
            posterior_dir=posterior_dir,
            train_genes=train_genes,
            val_genes=val_genes,
            parameters=parameters,
        )
    )
    return tasks
