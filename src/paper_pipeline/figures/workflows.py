from __future__ import annotations

import csv
from collections.abc import Callable
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

from ..dataprep.config import ConfigError, LoadedConfig
from ..dataprep.runner import preview_command
from ..dataprep.tasks import Task


SCRIPT_DIR = Path(__file__).resolve().parent / "scripts"
DEFAULT_HIGHLIGHT_GENESETS: tuple[str, ...] = (
    "HALLMARK_HEME_METABOLISM",
    "Hematopoiesisgenes",
    "mitotic_cell_cycle",
    "positive_macromolecule_synthesis",
)


@dataclass(frozen=True)
class FileIdMappingEntry:
    id1: str
    id2: str
    path1: str | Path
    path2: str | Path


@dataclass(frozen=True)
class CnmfTraitTarget:
    trait_file: str
    output_id: str


@dataclass(frozen=True)
class FileTarget:
    source_id: str
    source_path: str | Path


@dataclass(frozen=True)
class CorregulationTarget:
    program_a: str
    program_b: str
    output_id: str


@dataclass(frozen=True)
class ResolvedFiguresCnmf:
    config: LoadedConfig
    program_association_dir: Path
    cnmf_regulation_dir: Path
    output_dir: Path
    k: int
    trait_targets: tuple[CnmfTraitTarget, ...]
    plot_label_programs: tuple[int, ...]
    corregulation_targets: tuple[CorregulationTarget, ...]


@dataclass(frozen=True)
class ResolvedFiguresProgramHeatmap:
    config: LoadedConfig
    program_association_dir: Path
    output_dir: Path
    k: int
    trait_targets: tuple[CnmfTraitTarget, ...]
    output_id: str
    metrics: tuple[str, ...]


@dataclass(frozen=True)
class ResolvedFiguresBurdenVolcano:
    config: LoadedConfig
    geneset_dir: Path
    gene_map: Path
    spectra_path: Path
    output_dir: Path
    targets: tuple[FileTarget, ...]
    highlight_genesets: tuple[str, ...]
    data_genesets: tuple[str, ...]
    label_fdr_threshold: float
    line_fdr_threshold: float
    k: int
    top_n_program_genes: int


@dataclass(frozen=True)
class ResolvedFiguresGwasManhattan:
    config: LoadedConfig
    gene_annotation: Path
    geneset_dir: Path
    gene_map: Path
    spectra_path: Path
    output_dir: Path
    targets: tuple[FileTarget, ...]
    highlight_genesets: tuple[str, ...]
    data_genesets: tuple[str, ...]
    flank_bp: int
    label_p_threshold: float
    genomewide_threshold: float
    k: int
    top_n_program_genes: int


@dataclass(frozen=True)
class ResolvedFiguresPosteriorVolcano:
    config: LoadedConfig
    geneset_dir: Path
    gene_map: Path
    spectra_path: Path
    output_dir: Path
    targets: tuple[GeneLevelScatterTarget, ...]
    highlight_genesets: tuple[str, ...]
    data_genesets: tuple[str, ...]
    label_fdr_threshold: float
    line_fdr_threshold: float
    k: int
    top_n_program_genes: int


@dataclass(frozen=True)
class GeneLevelScatterTarget:
    source_id: str
    trait_stem: str
    posterior_path: Path


@dataclass(frozen=True)
class ResolvedFiguresGeneLevelScatter:
    config: LoadedConfig
    limma_path: Path
    shet_path: Path
    gene_map: Path
    output_dir: Path
    targets: tuple[GeneLevelScatterTarget, ...]
    highlight_genes: tuple[str, ...]
    top_n_labels: int
    y_limit: float


@dataclass(frozen=True)
class CrossTraitTarget:
    output_id: str
    x_id: str
    y_id: str
    x_label: str
    y_label: str
    posterior_path_x: Path
    posterior_path_y: Path


@dataclass(frozen=True)
class ResolvedFiguresCrossTrait:
    config: LoadedConfig
    gene_map: Path
    output_dir: Path
    targets: tuple[CrossTraitTarget, ...]
    highlight_genes: tuple[str, ...]
    top_n_labels: int


@dataclass(frozen=True)
class ResolvedFiguresGeneLevelQq:
    config: LoadedConfig
    limma_path: Path
    shet_path: Path
    gene_map: Path
    output_dir: Path
    targets: tuple[GeneLevelScatterTarget, ...]
    y_limit: float
    render_plot: bool


@dataclass(frozen=True)
class ResolvedFiguresProgramRankings:
    config: LoadedConfig
    program_association_dir: Path
    output_dir: Path
    k: int
    trait_targets: tuple[CnmfTraitTarget, ...]
    top_n: int


@dataclass(frozen=True)
class ResolvedFiguresCnmfProgramTopGenes:
    config: LoadedConfig
    spectra_path: Path
    gene_map: Path
    output_dir: Path
    k: int
    programs: tuple[int, ...]
    top_n: int
    output_id: str


@dataclass(frozen=True)
class ResolvedFiguresCnmfProgramEnrichment:
    config: LoadedConfig
    spectra_path: Path
    gene_map: Path
    geneset_dir: Path
    output_dir: Path
    k: int
    programs: tuple[int, ...]
    top_n: int
    genesets: tuple[str, ...]
    output_id: str


@dataclass(frozen=True)
class ResolvedFiguresTraitProgramGenePanel:
    config: LoadedConfig
    program_association_dir: Path
    regulation_dir: Path
    spectra_path: Path
    gene_map: Path
    output_dir: Path
    targets: tuple[GeneLevelScatterTarget, ...]
    k: int
    max_programs: int
    max_genes_per_side: int
    hit_abs_gamma_threshold: float
    loading_top_n: int
    regulator_fdr_threshold: float
    min_abs_score: float
    render_plot: bool


@dataclass(frozen=True)
class GwasLocusTarget:
    output_id: str
    source_id: str
    source_path: str | Path
    chrom: str
    start: int
    end: int
    flank_bp: int
    locus_label: str


@dataclass(frozen=True)
class ResolvedFiguresGwasLocusZoom:
    config: LoadedConfig
    gene_annotation: Path
    output_dir: Path
    targets: tuple[GwasLocusTarget, ...]
    genomewide_threshold: float
    label_top_n: int


@dataclass(frozen=True)
class ResolvedFiguresCrossTraitHeatmap:
    config: LoadedConfig
    gene_map: Path
    output_dir: Path
    targets: tuple[CnmfTraitTarget, ...]
    output_id: str
    method: str


def _string_list(value: Any, label: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ConfigError(f"{label} must be a list of strings")
    return tuple(item for item in value)


def _int_list(value: Any, label: str) -> tuple[int, ...]:
    if value is None:
        return ()
    if not isinstance(value, list):
        raise ConfigError(f"{label} must be a list of integers")
    try:
        return tuple(int(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{label} must be a list of integers") from exc


def _mapping_list(value: Any, label: str) -> tuple[dict[str, Any], ...]:
    if value is None:
        return ()
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise ConfigError(f"{label} must be a list of mappings")
    return tuple(value)


def _workflow_root(config: LoadedConfig) -> dict[str, Any]:
    figures = config.workflow("figures")
    if not isinstance(figures, dict):
        raise ConfigError("workflows.figures must be a mapping")
    return figures


def _workflow_mapping(config: LoadedConfig, key: str) -> dict[str, Any]:
    figures = _workflow_root(config)
    workflow = figures.get(key)
    if workflow is None:
        raise ConfigError(f"workflows.figures.{key} is not configured")
    if not isinstance(workflow, dict):
        raise ConfigError(f"workflows.figures.{key} must be a mapping")

    # Merge shared_inputs / shared_parameters from the figures root into this
    # sub-workflow.  Sub-workflow values take precedence over shared values.
    shared_inputs = figures.get("shared_inputs")
    if isinstance(shared_inputs, dict):
        merged_inputs = {**shared_inputs, **_validate_mapping(workflow.get("inputs"), f"workflows.figures.{key}.inputs")}
        workflow = {**workflow, "inputs": merged_inputs}

    shared_parameters = figures.get("shared_parameters")
    if isinstance(shared_parameters, dict):
        merged_parameters = {**shared_parameters, **_validate_mapping(workflow.get("parameters"), f"workflows.figures.{key}.parameters")}
        workflow = {**workflow, "parameters": merged_parameters}

    return workflow


def _validate_mapping(value: Any, label: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(f"{label} must be a mapping")
    return value


def _normalize_output_id(value: str, label: str) -> str:
    normalized = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value.strip())
    if not normalized:
        raise ConfigError(f"{label} must not be empty")
    return normalized


def _resolve_external_path(base_dir: Path, raw_path: str) -> str | Path:
    if raw_path.startswith("/") and not raw_path.startswith("//"):
        return raw_path
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def _load_file_id_map(config: LoadedConfig, mapping_path_value: str | Path | None) -> tuple[FileIdMappingEntry, ...]:
    mapping_path = config.resolve_path(mapping_path_value)
    if mapping_path is None:
        return ()
    if not mapping_path.exists():
        raise ConfigError(f"File ID mapping does not exist: {mapping_path}")

    rows: list[FileIdMappingEntry] = []
    with mapping_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"id1", "id2", "path1", "path2"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            missing = ", ".join(sorted(required - set(reader.fieldnames or [])))
            raise ConfigError(f"{mapping_path} must contain TSV columns: id1, id2, path1, path2; missing: {missing}")
        for index, row in enumerate(reader, start=2):
            id1 = (row.get("id1") or "").strip()
            id2 = (row.get("id2") or "").strip()
            path1 = (row.get("path1") or "").strip()
            path2 = (row.get("path2") or "").strip()
            if not (id1 and id2 and path1 and path2):
                raise ConfigError(f"{mapping_path}:{index} must provide non-empty id1, id2, path1, and path2")
            rows.append(
                FileIdMappingEntry(
                    id1=id1,
                    id2=id2,
                    path1=_resolve_external_path(mapping_path.parent, path1),
                    path2=_resolve_external_path(mapping_path.parent, path2),
                )
            )
    return tuple(rows)


def _load_two_column_name_map(
    config: LoadedConfig,
    mapping_path_value: str | Path | None,
    *,
    label: str,
) -> dict[str, str]:
    mapping_path = config.resolve_path(mapping_path_value)
    if mapping_path is None:
        return {}
    if not mapping_path.exists():
        raise ConfigError(f"{label} does not exist: {mapping_path}")

    rows: dict[str, str] = {}
    with mapping_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for index, row in enumerate(reader, start=1):
            if not row or all(not str(cell).strip() for cell in row):
                continue
            if index == 1 and row[0].strip().lower() in {"id", "id2", "lof_id"}:
                continue
            if len(row) < 2:
                raise ConfigError(f"{label}:{index} must contain at least 2 tab-separated columns")
            source_id = row[0].strip()
            output_name = row[1].strip()
            if not source_id or not output_name:
                raise ConfigError(f"{label}:{index} must provide non-empty source id and filename")
            rows[source_id] = output_name
    return rows


def _targets_from_mapping(
    mapping_entries: tuple[FileIdMappingEntry, ...],
    ids: tuple[str, ...],
    *,
    id_field: str,
    path_field: str,
    label: str,
) -> tuple[FileTarget, ...]:
    if not ids:
        return ()
    if not mapping_entries:
        raise ConfigError(f"{label} requires workflows.figures.<workflow>.inputs.file_id_map")

    targets: list[FileTarget] = []
    for source_id in ids:
        match = next((entry for entry in mapping_entries if getattr(entry, id_field) == source_id), None)
        if match is None:
            raise ConfigError(f"{label} value '{source_id}' was not found in file_id_map column {id_field}")
        targets.append(
            FileTarget(
                source_id=_normalize_output_id(source_id, label),
                source_path=getattr(match, path_field),
            )
        )
    return tuple(targets)


def _targets_from_paths(paths: tuple[str, ...], label: str, *, resolve: Callable[[str], Path | None]) -> tuple[FileTarget, ...]:
    targets: list[FileTarget] = []
    for raw_path in paths:
        source_path = resolve(raw_path)
        if source_path is None:
            raise ConfigError(f"{label} path '{raw_path}' could not be resolved")
        output_id = _file_label(raw_path)
        targets.append(FileTarget(source_id=output_id, source_path=source_path))
    return tuple(targets)


def _mapping_entries_by_id(
    mapping_entries: tuple[FileIdMappingEntry, ...],
    *,
    id_field: str,
) -> dict[str, FileIdMappingEntry]:
    return {str(getattr(entry, id_field)): entry for entry in mapping_entries}


def _posterior_filename_from_trait_source(source_path: str | Path) -> str:
    return f"{_file_label(str(source_path))}.per_gene_estimates.tsv"


def _correlation_filename_from_trait_source(source_path: str | Path) -> str:
    return f"{_posterior_filename_from_trait_source(source_path)}_geneRegulation_correlation.txt"


def _posterior_filename_for_entry(entry: FileIdMappingEntry, overrides: dict[str, str]) -> str:
    override = overrides.get(entry.id2)
    if override:
        return override
    return _posterior_filename_from_trait_source(entry.path2)


def _trait_targets_from_parameters(
    parameters: dict[str, Any],
    label: str,
) -> tuple[CnmfTraitTarget, ...]:
    target_mappings = _mapping_list(parameters.get("trait_targets"), f"{label}.trait_targets")
    targets: list[CnmfTraitTarget] = []
    if target_mappings:
        for index, target in enumerate(target_mappings, start=1):
            trait_file = target.get("trait_file")
            trait_id = target.get("trait_id")
            if not isinstance(trait_file, str) or not trait_file:
                raise ConfigError(f"{label}.trait_targets[{index}] must provide non-empty trait_file")
            if trait_id is None:
                trait_id = _trait_stem(trait_file)
            if not isinstance(trait_id, str) or not trait_id:
                raise ConfigError(f"{label}.trait_targets[{index}] trait_id must be a non-empty string")
            targets.append(
                CnmfTraitTarget(
                    trait_file=trait_file,
                    output_id=_normalize_output_id(trait_id, f"{label}.trait_targets[{index}].trait_id"),
                )
            )
        return tuple(targets)

    trait_files = _string_list(parameters.get("trait_files"), f"{label}.trait_files")
    if not trait_files:
        trait_files = ("Backman_2021_86.per_gene_estimates.tsv",)
    return tuple(
        CnmfTraitTarget(
            trait_file=trait_file,
            output_id=_normalize_output_id(_trait_stem(trait_file), f"{label}.trait_files"),
        )
        for trait_file in trait_files
    )


def _resolve_figures_cnmf(config: LoadedConfig) -> ResolvedFiguresCnmf:
    workflow = _workflow_mapping(config, "cnmf")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.cnmf.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.cnmf.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.cnmf.parameters")

    k = int(parameters.get("k", 60))
    plot_label_programs = _int_list(
        parameters.get("plot_label_programs"),
        "workflows.figures.cnmf.parameters.plot_label_programs",
    ) or (4, 16, 25, 40)
    corregulation_pairs_value = parameters.get("corregulation_pairs")
    if corregulation_pairs_value is None:
        corregulation_pairs_raw = ({"program_a": "P25", "program_b": "P16"},)
    else:
        corregulation_pairs_raw = _mapping_list(
            corregulation_pairs_value,
            "workflows.figures.cnmf.parameters.corregulation_pairs",
        )

    pairs: list[CorregulationTarget] = []
    for index, pair in enumerate(corregulation_pairs_raw, start=1):
        program_a = pair.get("program_a")
        program_b = pair.get("program_b")
        output_id = pair.get("output_id")
        if not isinstance(program_a, str) or not isinstance(program_b, str):
            raise ConfigError(
                f"workflows.figures.cnmf.parameters.corregulation_pairs[{index}] must provide string program_a/program_b"
            )
        pair_id = output_id if isinstance(output_id, str) and output_id else f"{program_a}__{program_b}"
        pairs.append(
            CorregulationTarget(
                program_a=program_a,
                program_b=program_b,
                output_id=_normalize_output_id(pair_id, f"workflows.figures.cnmf.parameters.corregulation_pairs[{index}].output_id"),
            )
        )

    return ResolvedFiguresCnmf(
        config=config,
        program_association_dir=config.resolve_path_or_artifact(
            inputs.get("program_association_dir"),
            "perturbseq",
            "cnmf_genomewide",
            "trait_association",
            "K562GW",
            "ProgramLevel",
        ),
        cnmf_regulation_dir=config.resolve_path_or_artifact(
            inputs.get("cnmf_regulation_dir"),
            "perturbseq",
            "cnmf_genomewide",
            "cNMF_regulation",
            "K562GW",
        ),
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "cnmf"),
        k=k,
        trait_targets=_trait_targets_from_parameters(parameters, "workflows.figures.cnmf.parameters"),
        plot_label_programs=plot_label_programs,
        corregulation_targets=tuple(pairs),
    )


def _resolve_figures_burden_volcano(config: LoadedConfig) -> ResolvedFiguresBurdenVolcano:
    workflow = _workflow_mapping(config, "burden_volcano")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.burden_volcano.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.burden_volcano.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.burden_volcano.parameters")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    lof_ids = _string_list(parameters.get("lof_ids"), "workflows.figures.burden_volcano.parameters.lof_ids")
    burden_files = _string_list(
        parameters.get("burden_files"),
        "workflows.figures.burden_volcano.parameters.burden_files",
    )

    if lof_ids:
        targets = _targets_from_mapping(
            mapping_entries,
            lof_ids,
            id_field="id2",
            path_field="path2",
            label="workflows.figures.burden_volcano.parameters.lof_ids",
        )
    else:
        if not burden_files:
            burden_files = ("Backman_2021_86_M1_001.summary_statistics.csv",)
        raw_burden_dir = config.resolve_path(inputs.get("raw_burden_dir")) or config.project_root / "data" / "LoF" / "raw_burden"
        targets = tuple(
            FileTarget(
                source_id=_file_label(path),
                source_path=raw_burden_dir / path,
            )
            for path in burden_files
        )

    highlight_genesets = _string_list(
        parameters.get("highlight_genesets"),
        "workflows.figures.burden_volcano.parameters.highlight_genesets",
    ) or DEFAULT_HIGHLIGHT_GENESETS
    data_genesets = _string_list(
        parameters.get("data_genesets"),
        "workflows.figures.burden_volcano.parameters.data_genesets",
    ) or highlight_genesets
    k = int(parameters.get("k", 60))

    return ResolvedFiguresBurdenVolcano(
        config=config,
        geneset_dir=config.resolve_path(inputs.get("geneset_dir")) or config.project_root / "data" / "geneset",
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        spectra_path=config.resolve_path(inputs.get("spectra_path"))
        or config.resolve_path_or_artifact(
            inputs.get("spectra_path"),
            "perturbseq",
            "cnmf_genomewide",
            "cNMF",
            "cNMF_all",
            f"cNMF_all.gene_spectra_score.k_{k}.dt_0_5.txt",
        ),
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "burden_volcano"),
        targets=targets,
        highlight_genesets=highlight_genesets,
        data_genesets=data_genesets,
        label_fdr_threshold=float(parameters.get("label_fdr_threshold", 0.01)),
        line_fdr_threshold=float(parameters.get("line_fdr_threshold", 0.1)),
        k=k,
        top_n_program_genes=int(parameters.get("top_n_program_genes", 100)),
    )


def _resolve_figures_posterior_volcano(config: LoadedConfig) -> ResolvedFiguresPosteriorVolcano:
    workflow = _workflow_mapping(config, "posterior_volcano")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.posterior_volcano.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.posterior_volcano.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.posterior_volcano.parameters")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    lof_ids = _string_list(parameters.get("lof_ids"), "workflows.figures.posterior_volcano.parameters.lof_ids")
    if not lof_ids:
        raise ConfigError("workflows.figures.posterior_volcano.parameters.lof_ids must not be empty")
    if not mapping_entries:
        raise ConfigError("workflows.figures.posterior_volcano.inputs.file_id_map is required")

    id_lookup = _mapping_entries_by_id(mapping_entries, id_field="id2")
    posterior_dir = config.resolve_path_or_artifact(inputs.get("posterior_dir"), "genebayes", "posterior")
    posterior_name_map = _load_two_column_name_map(
        config,
        inputs.get("posterior_name_map"),
        label="workflows.figures.posterior_volcano.inputs.posterior_name_map",
    )

    targets: list[GeneLevelScatterTarget] = []
    for source_id in lof_ids:
        entry = id_lookup.get(source_id)
        if entry is None:
            raise ConfigError(
                f"workflows.figures.posterior_volcano.parameters.lof_ids value '{source_id}' "
                "was not found in file_id_map column id2"
            )
        trait_stem = _file_label(str(entry.path2))
        targets.append(
            GeneLevelScatterTarget(
                source_id=_normalize_output_id(source_id, "workflows.figures.posterior_volcano.parameters.lof_ids"),
                trait_stem=trait_stem,
                posterior_path=posterior_dir / _posterior_filename_for_entry(entry, posterior_name_map),
            )
        )

    highlight_genesets = _string_list(
        parameters.get("highlight_genesets"),
        "workflows.figures.posterior_volcano.parameters.highlight_genesets",
    ) or DEFAULT_HIGHLIGHT_GENESETS
    data_genesets = _string_list(
        parameters.get("data_genesets"),
        "workflows.figures.posterior_volcano.parameters.data_genesets",
    ) or highlight_genesets
    k = int(parameters.get("k", 60))

    return ResolvedFiguresPosteriorVolcano(
        config=config,
        geneset_dir=config.resolve_path(inputs.get("geneset_dir")) or config.project_root / "data" / "geneset",
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        spectra_path=config.resolve_path(inputs.get("spectra_path"))
        or config.resolve_path_or_artifact(
            inputs.get("spectra_path"),
            "perturbseq",
            "cnmf_genomewide",
            "cNMF",
            "cNMF_all",
            f"cNMF_all.gene_spectra_score.k_{k}.dt_0_5.txt",
        ),
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "posterior_volcano"),
        targets=tuple(targets),
        highlight_genesets=highlight_genesets,
        data_genesets=data_genesets,
        label_fdr_threshold=float(parameters.get("label_fdr_threshold", 0.01)),
        line_fdr_threshold=float(parameters.get("line_fdr_threshold", 0.1)),
        k=k,
        top_n_program_genes=int(parameters.get("top_n_program_genes", 100)),
    )


def _resolve_figures_program_heatmap(config: LoadedConfig) -> ResolvedFiguresProgramHeatmap:
    workflow = _workflow_mapping(config, "program_heatmap")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.program_heatmap.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.program_heatmap.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.program_heatmap.parameters")

    metric_values = _string_list(parameters.get("metrics"), "workflows.figures.program_heatmap.parameters.metrics") or (
        "program_score",
        "regulator_score",
    )
    allowed_metrics = {"program_score", "regulator_score"}
    invalid = [metric for metric in metric_values if metric not in allowed_metrics]
    if invalid:
        raise ConfigError(
            "workflows.figures.program_heatmap.parameters.metrics contains unsupported values: "
            + ", ".join(invalid)
        )

    trait_targets = _trait_targets_from_parameters(parameters, "workflows.figures.program_heatmap.parameters")
    output_id = parameters.get("output_id")
    if output_id is None:
        output_id = "__".join(target.output_id for target in trait_targets)
    if not isinstance(output_id, str) or not output_id:
        raise ConfigError("workflows.figures.program_heatmap.parameters.output_id must be a non-empty string")

    return ResolvedFiguresProgramHeatmap(
        config=config,
        program_association_dir=config.resolve_path_or_artifact(
            inputs.get("program_association_dir"),
            "perturbseq",
            "cnmf_genomewide",
            "trait_association",
            "K562GW",
            "ProgramLevel",
        ),
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "program_heatmap"),
        k=int(parameters.get("k", 60)),
        trait_targets=trait_targets,
        output_id=_normalize_output_id(output_id, "workflows.figures.program_heatmap.parameters.output_id"),
        metrics=metric_values,
    )


def _resolve_figures_gwas_manhattan(config: LoadedConfig) -> ResolvedFiguresGwasManhattan:
    workflow = _workflow_mapping(config, "gwas_manhattan")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.gwas_manhattan.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.gwas_manhattan.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.gwas_manhattan.parameters")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    gwas_ids = _string_list(parameters.get("gwas_ids"), "workflows.figures.gwas_manhattan.parameters.gwas_ids")
    gwas_files = _string_list(
        parameters.get("gwas_files"),
        "workflows.figures.gwas_manhattan.parameters.gwas_files",
    )
    gwas_dir = config.resolve_path(inputs.get("gwas_dir")) or config.project_root / "data" / "GWAS"

    if gwas_ids:
        targets = _targets_from_mapping(
            mapping_entries,
            gwas_ids,
            id_field="id1",
            path_field="path1",
            label="workflows.figures.gwas_manhattan.parameters.gwas_ids",
        )
    else:
        if not gwas_files:
            gwas_files = ("30050_irnt.tsv.gz",)
        targets = tuple(
            FileTarget(
                source_id=_file_label(path),
                source_path=gwas_dir / path,
            )
            for path in gwas_files
        )

    highlight_genesets = _string_list(
        parameters.get("highlight_genesets"),
        "workflows.figures.gwas_manhattan.parameters.highlight_genesets",
    ) or DEFAULT_HIGHLIGHT_GENESETS
    data_genesets = _string_list(
        parameters.get("data_genesets"),
        "workflows.figures.gwas_manhattan.parameters.data_genesets",
    ) or highlight_genesets
    k = int(parameters.get("k", 60))

    return ResolvedFiguresGwasManhattan(
        config=config,
        gene_annotation=config.resolve_path(inputs.get("gene_annotation")) or gwas_dir / "genes.protein_coding.v39.gtf",
        geneset_dir=config.resolve_path(inputs.get("geneset_dir")) or config.project_root / "data" / "geneset",
        gene_map=config.resolve_path(inputs.get("gene_map")) or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        spectra_path=config.resolve_path(inputs.get("spectra_path"))
        or config.resolve_path_or_artifact(
            inputs.get("spectra_path"),
            "perturbseq",
            "cnmf_genomewide",
            "cNMF",
            "cNMF_all",
            f"cNMF_all.gene_spectra_score.k_{k}.dt_0_5.txt",
        ),
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "gwas_manhattan"),
        targets=targets,
        highlight_genesets=highlight_genesets,
        data_genesets=data_genesets,
        flank_bp=int(parameters.get("flank_bp", 50000)),
        label_p_threshold=float(parameters.get("label_p_threshold", 1e-30)),
        genomewide_threshold=float(parameters.get("genomewide_threshold", 5e-8)),
        k=k,
        top_n_program_genes=int(parameters.get("top_n_program_genes", 100)),
    )


def _resolve_figures_gene_level_scatter(config: LoadedConfig) -> ResolvedFiguresGeneLevelScatter:
    workflow = _workflow_mapping(config, "gene_level_scatter")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.gene_level_scatter.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.gene_level_scatter.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.gene_level_scatter.parameters")

    limma_path = config.resolve_path(inputs.get("limma_path"))
    if limma_path is None:
        raise ConfigError("workflows.figures.gene_level_scatter.inputs.limma_path is required")
    shet_path = config.resolve_path(inputs.get("shet_path"))
    if shet_path is None:
        raise ConfigError("workflows.figures.gene_level_scatter.inputs.shet_path is required")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    lof_ids = _string_list(parameters.get("lof_ids"), "workflows.figures.gene_level_scatter.parameters.lof_ids")
    if not lof_ids:
        raise ConfigError("workflows.figures.gene_level_scatter.parameters.lof_ids must not be empty")
    if not mapping_entries:
        raise ConfigError("workflows.figures.gene_level_scatter.inputs.file_id_map is required")

    id_lookup = _mapping_entries_by_id(mapping_entries, id_field="id2")
    posterior_dir = config.resolve_path_or_artifact(inputs.get("posterior_dir"), "genebayes", "posterior")
    posterior_name_map = _load_two_column_name_map(
        config,
        inputs.get("posterior_name_map"),
        label="workflows.figures.gene_level_scatter.inputs.posterior_name_map",
    )
    targets: list[GeneLevelScatterTarget] = []
    for source_id in lof_ids:
        entry = id_lookup.get(source_id)
        if entry is None:
            raise ConfigError(
                f"workflows.figures.gene_level_scatter.parameters.lof_ids value '{source_id}' "
                "was not found in file_id_map column id2"
            )
        trait_stem = _file_label(str(entry.path2))
        targets.append(
            GeneLevelScatterTarget(
                source_id=_normalize_output_id(source_id, "workflows.figures.gene_level_scatter.parameters.lof_ids"),
                trait_stem=trait_stem,
                posterior_path=posterior_dir / _posterior_filename_for_entry(entry, posterior_name_map),
            )
        )

    return ResolvedFiguresGeneLevelScatter(
        config=config,
        limma_path=limma_path,
        shet_path=shet_path,
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "gene_level_scatter"),
        targets=tuple(targets),
        highlight_genes=_string_list(
            parameters.get("highlight_genes"),
            "workflows.figures.gene_level_scatter.parameters.highlight_genes",
        ),
        top_n_labels=int(parameters.get("top_n_labels", 8)),
        y_limit=float(parameters.get("y_limit", 8.0)),
    )


def _resolve_figures_cross_trait(config: LoadedConfig) -> ResolvedFiguresCrossTrait:
    workflow = _workflow_mapping(config, "cross_trait")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.cross_trait.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.cross_trait.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.cross_trait.parameters")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    if not mapping_entries:
        raise ConfigError("workflows.figures.cross_trait.inputs.file_id_map is required")

    id_lookup = _mapping_entries_by_id(mapping_entries, id_field="id2")
    posterior_dir = config.resolve_path_or_artifact(inputs.get("posterior_dir"), "genebayes", "posterior")
    posterior_name_map = _load_two_column_name_map(
        config,
        inputs.get("posterior_name_map"),
        label="workflows.figures.cross_trait.inputs.posterior_name_map",
    )
    pair_mappings = _mapping_list(parameters.get("lof_pairs"), "workflows.figures.cross_trait.parameters.lof_pairs")
    lof_ids = _string_list(parameters.get("lof_ids"), "workflows.figures.cross_trait.parameters.lof_ids")

    pair_specs: list[dict[str, str]] = []
    if pair_mappings:
        for index, pair in enumerate(pair_mappings, start=1):
            id_x = pair.get("id_x")
            id_y = pair.get("id_y")
            if not isinstance(id_x, str) or not id_x or not isinstance(id_y, str) or not id_y:
                raise ConfigError(
                    f"workflows.figures.cross_trait.parameters.lof_pairs[{index}] must provide non-empty id_x and id_y"
                )
            pair_specs.append(
                {
                    "id_x": id_x,
                    "id_y": id_y,
                    "label_x": str(pair.get("label_x") or id_x),
                    "label_y": str(pair.get("label_y") or id_y),
                    "output_id": str(pair.get("output_id") or f"{id_x}__{id_y}"),
                }
            )
    else:
        if len(lof_ids) < 2:
            raise ConfigError(
                "workflows.figures.cross_trait.parameters.lof_pairs is required unless lof_ids provides at least two IDs"
            )
        for id_x, id_y in combinations(lof_ids, 2):
            pair_specs.append(
                {
                    "id_x": id_x,
                    "id_y": id_y,
                    "label_x": id_x,
                    "label_y": id_y,
                    "output_id": f"{id_x}__{id_y}",
                }
            )

    targets: list[CrossTraitTarget] = []
    for pair in pair_specs:
        entry_x = id_lookup.get(pair["id_x"])
        entry_y = id_lookup.get(pair["id_y"])
        if entry_x is None:
            raise ConfigError(
                f"workflows.figures.cross_trait parameter id_x '{pair['id_x']}' was not found in file_id_map column id2"
            )
        if entry_y is None:
            raise ConfigError(
                f"workflows.figures.cross_trait parameter id_y '{pair['id_y']}' was not found in file_id_map column id2"
            )
        targets.append(
            CrossTraitTarget(
                output_id=_normalize_output_id(pair["output_id"], "workflows.figures.cross_trait.parameters.output_id"),
                x_id=_normalize_output_id(pair["id_x"], "workflows.figures.cross_trait.parameters.id_x"),
                y_id=_normalize_output_id(pair["id_y"], "workflows.figures.cross_trait.parameters.id_y"),
                x_label=pair["label_x"],
                y_label=pair["label_y"],
                posterior_path_x=posterior_dir / _posterior_filename_for_entry(entry_x, posterior_name_map),
                posterior_path_y=posterior_dir / _posterior_filename_for_entry(entry_y, posterior_name_map),
            )
        )

    return ResolvedFiguresCrossTrait(
        config=config,
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "cross_trait"),
        targets=tuple(targets),
        highlight_genes=_string_list(
            parameters.get("highlight_genes"),
            "workflows.figures.cross_trait.parameters.highlight_genes",
        ),
        top_n_labels=int(parameters.get("top_n_labels", 12)),
    )


def _resolve_figures_gene_level_qq(config: LoadedConfig) -> ResolvedFiguresGeneLevelQq:
    workflow = _workflow_mapping(config, "gene_level_qq")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.gene_level_qq.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.gene_level_qq.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.gene_level_qq.parameters")

    limma_path = config.resolve_path(inputs.get("limma_path"))
    if limma_path is None:
        raise ConfigError("workflows.figures.gene_level_qq.inputs.limma_path is required")
    shet_path = config.resolve_path(inputs.get("shet_path"))
    if shet_path is None:
        raise ConfigError("workflows.figures.gene_level_qq.inputs.shet_path is required")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    lof_ids = _string_list(parameters.get("lof_ids"), "workflows.figures.gene_level_qq.parameters.lof_ids")
    if not lof_ids:
        raise ConfigError("workflows.figures.gene_level_qq.parameters.lof_ids must not be empty")
    if not mapping_entries:
        raise ConfigError("workflows.figures.gene_level_qq.inputs.file_id_map is required")

    id_lookup = _mapping_entries_by_id(mapping_entries, id_field="id2")
    posterior_dir = config.resolve_path_or_artifact(inputs.get("posterior_dir"), "genebayes", "posterior")
    posterior_name_map = _load_two_column_name_map(
        config,
        inputs.get("posterior_name_map"),
        label="workflows.figures.gene_level_qq.inputs.posterior_name_map",
    )
    targets: list[GeneLevelScatterTarget] = []
    for source_id in lof_ids:
        entry = id_lookup.get(source_id)
        if entry is None:
            raise ConfigError(
                f"workflows.figures.gene_level_qq.parameters.lof_ids value '{source_id}' "
                "was not found in file_id_map column id2"
            )
        trait_stem = _file_label(str(entry.path2))
        targets.append(
            GeneLevelScatterTarget(
                source_id=_normalize_output_id(source_id, "workflows.figures.gene_level_qq.parameters.lof_ids"),
                trait_stem=trait_stem,
                posterior_path=posterior_dir / _posterior_filename_for_entry(entry, posterior_name_map),
            )
        )

    return ResolvedFiguresGeneLevelQq(
        config=config,
        limma_path=limma_path,
        shet_path=shet_path,
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "gene_level_qq"),
        targets=tuple(targets),
        y_limit=float(parameters.get("y_limit", 10.0)),
        render_plot=bool(parameters.get("render_plot", True)),
    )


def _resolve_figures_program_rankings(config: LoadedConfig) -> ResolvedFiguresProgramRankings:
    workflow = _workflow_mapping(config, "program_rankings")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.program_rankings.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.program_rankings.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.program_rankings.parameters")

    return ResolvedFiguresProgramRankings(
        config=config,
        program_association_dir=config.resolve_path_or_artifact(
            inputs.get("program_association_dir"),
            "perturbseq",
            "cnmf_genomewide",
            "trait_association",
            "K562GW",
            "ProgramLevel",
        ),
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "program_rankings"),
        k=int(parameters.get("k", 60)),
        trait_targets=_trait_targets_from_parameters(parameters, "workflows.figures.program_rankings.parameters"),
        top_n=int(parameters.get("top_n", 12)),
    )


def _resolve_figures_gwas_locus_zoom(config: LoadedConfig) -> ResolvedFiguresGwasLocusZoom:
    workflow = _workflow_mapping(config, "gwas_locus_zoom")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.gwas_locus_zoom.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.gwas_locus_zoom.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.gwas_locus_zoom.parameters")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    id_lookup = _mapping_entries_by_id(mapping_entries, id_field="id1")
    target_mappings = _mapping_list(
        parameters.get("locus_targets"),
        "workflows.figures.gwas_locus_zoom.parameters.locus_targets",
    )
    if not target_mappings:
        raise ConfigError("workflows.figures.gwas_locus_zoom.parameters.locus_targets must not be empty")

    targets: list[GwasLocusTarget] = []
    for index, target in enumerate(target_mappings, start=1):
        gwas_id = target.get("gwas_id")
        gwas_file = target.get("gwas_file")
        chrom = target.get("chrom")
        start = target.get("start")
        end = target.get("end")
        if not isinstance(chrom, str) or not chrom:
            raise ConfigError(
                f"workflows.figures.gwas_locus_zoom.parameters.locus_targets[{index}] must provide non-empty chrom"
            )
        if start is None or end is None:
            raise ConfigError(
                f"workflows.figures.gwas_locus_zoom.parameters.locus_targets[{index}] must provide start and end"
            )
        try:
            start_i = int(start)
            end_i = int(end)
        except (TypeError, ValueError) as exc:
            raise ConfigError(
                f"workflows.figures.gwas_locus_zoom.parameters.locus_targets[{index}] start/end must be integers"
            ) from exc

        if isinstance(gwas_id, str) and gwas_id:
            if not mapping_entries:
                raise ConfigError("workflows.figures.gwas_locus_zoom.inputs.file_id_map is required when using gwas_id")
            entry = id_lookup.get(gwas_id)
            if entry is None:
                raise ConfigError(
                    f"workflows.figures.gwas_locus_zoom.parameters.locus_targets[{index}].gwas_id '{gwas_id}' "
                    "was not found in file_id_map column id1"
                )
            source_id = _normalize_output_id(gwas_id, "workflows.figures.gwas_locus_zoom.parameters.gwas_id")
            source_path: str | Path = entry.path1
        elif isinstance(gwas_file, str) and gwas_file:
            gwas_dir = config.resolve_path(inputs.get("gwas_dir")) or config.project_root / "data" / "GWAS"
            source_path = gwas_dir / gwas_file
            source_id = _file_label(gwas_file)
        else:
            raise ConfigError(
                f"workflows.figures.gwas_locus_zoom.parameters.locus_targets[{index}] must provide gwas_id or gwas_file"
            )

        flank_bp = int(target.get("flank_bp", parameters.get("flank_bp", 250000)))
        locus_label = str(target.get("locus_label") or f"chr{chrom}:{start_i}-{end_i}")
        output_id = str(target.get("output_id") or f"{source_id}__chr{chrom}_{start_i}_{end_i}")
        targets.append(
            GwasLocusTarget(
                output_id=_normalize_output_id(
                    output_id,
                    f"workflows.figures.gwas_locus_zoom.parameters.locus_targets[{index}].output_id",
                ),
                source_id=source_id,
                source_path=source_path,
                chrom=chrom,
                start=start_i,
                end=end_i,
                flank_bp=flank_bp,
                locus_label=locus_label,
            )
        )

    gwas_dir = config.resolve_path(inputs.get("gwas_dir")) or config.project_root / "data" / "GWAS"
    return ResolvedFiguresGwasLocusZoom(
        config=config,
        gene_annotation=config.resolve_path(inputs.get("gene_annotation")) or gwas_dir / "genes.protein_coding.v39.gtf",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "gwas_locus_zoom"),
        targets=tuple(targets),
        genomewide_threshold=float(parameters.get("genomewide_threshold", 5e-8)),
        label_top_n=int(parameters.get("label_top_n", 6)),
    )


def _resolve_figures_cross_trait_heatmap(config: LoadedConfig) -> ResolvedFiguresCrossTraitHeatmap:
    workflow = _workflow_mapping(config, "cross_trait_heatmap")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.cross_trait_heatmap.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.cross_trait_heatmap.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.cross_trait_heatmap.parameters")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    if not mapping_entries:
        raise ConfigError("workflows.figures.cross_trait_heatmap.inputs.file_id_map is required")
    lof_ids = _string_list(parameters.get("lof_ids"), "workflows.figures.cross_trait_heatmap.parameters.lof_ids")
    if len(lof_ids) < 2:
        raise ConfigError("workflows.figures.cross_trait_heatmap.parameters.lof_ids must contain at least two IDs")

    id_lookup = _mapping_entries_by_id(mapping_entries, id_field="id2")
    posterior_dir = config.resolve_path_or_artifact(inputs.get("posterior_dir"), "genebayes", "posterior")
    posterior_name_map = _load_two_column_name_map(
        config,
        inputs.get("posterior_name_map"),
        label="workflows.figures.cross_trait_heatmap.inputs.posterior_name_map",
    )
    targets: list[CnmfTraitTarget] = []
    for source_id in lof_ids:
        entry = id_lookup.get(source_id)
        if entry is None:
            raise ConfigError(
                f"workflows.figures.cross_trait_heatmap.parameters.lof_ids value '{source_id}' "
                "was not found in file_id_map column id2"
            )
        targets.append(
            CnmfTraitTarget(
                trait_file=str(posterior_dir / _posterior_filename_for_entry(entry, posterior_name_map)),
                output_id=_normalize_output_id(source_id, "workflows.figures.cross_trait_heatmap.parameters.lof_ids"),
            )
        )

    output_id = parameters.get("output_id")
    if output_id is None:
        output_id = "__".join(target.output_id for target in targets)
    method = str(parameters.get("method", "pearson")).lower()
    if method not in {"pearson", "spearman"}:
        raise ConfigError("workflows.figures.cross_trait_heatmap.parameters.method must be 'pearson' or 'spearman'")

    return ResolvedFiguresCrossTraitHeatmap(
        config=config,
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "cross_trait_heatmap"),
        targets=tuple(targets),
        output_id=_normalize_output_id(output_id, "workflows.figures.cross_trait_heatmap.parameters.output_id"),
        method=method,
    )


def _resolve_figures_cnmf_program_top_genes(config: LoadedConfig) -> ResolvedFiguresCnmfProgramTopGenes:
    workflow = _workflow_mapping(config, "cnmf_program_top_genes")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.cnmf_program_top_genes.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.cnmf_program_top_genes.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.cnmf_program_top_genes.parameters")

    k = int(parameters.get("k", 60))
    return ResolvedFiguresCnmfProgramTopGenes(
        config=config,
        spectra_path=config.resolve_path_or_artifact(
            inputs.get("spectra_path"),
            "perturbseq",
            "cnmf_genomewide",
            "cNMF",
            "cNMF_all",
            f"cNMF_all.gene_spectra_score.k_{k}.dt_0_5.txt",
        ),
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "cnmf_program_top_genes"),
        k=k,
        programs=_int_list(parameters.get("programs"), "workflows.figures.cnmf_program_top_genes.parameters.programs"),
        top_n=int(parameters.get("top_n", 15)),
        output_id=_normalize_output_id(parameters.get("output_id", "top_genes"), "workflows.figures.cnmf_program_top_genes.parameters.output_id"),
    )


def _resolve_figures_cnmf_program_enrichment(config: LoadedConfig) -> ResolvedFiguresCnmfProgramEnrichment:
    workflow = _workflow_mapping(config, "cnmf_program_enrichment")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.cnmf_program_enrichment.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.cnmf_program_enrichment.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.cnmf_program_enrichment.parameters")

    k = int(parameters.get("k", 60))
    genesets = _string_list(parameters.get("genesets"), "workflows.figures.cnmf_program_enrichment.parameters.genesets")
    if not genesets:
        raise ConfigError("workflows.figures.cnmf_program_enrichment.parameters.genesets must not be empty")

    return ResolvedFiguresCnmfProgramEnrichment(
        config=config,
        spectra_path=config.resolve_path_or_artifact(
            inputs.get("spectra_path"),
            "perturbseq",
            "cnmf_genomewide",
            "cNMF",
            "cNMF_all",
            f"cNMF_all.gene_spectra_score.k_{k}.dt_0_5.txt",
        ),
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        geneset_dir=config.resolve_path(inputs.get("geneset_dir")) or config.project_root / "data" / "geneset",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "cnmf_program_enrichment"),
        k=k,
        programs=_int_list(parameters.get("programs"), "workflows.figures.cnmf_program_enrichment.parameters.programs"),
        top_n=int(parameters.get("top_n", 200)),
        genesets=genesets,
        output_id=_normalize_output_id(parameters.get("output_id", "enrichment"), "workflows.figures.cnmf_program_enrichment.parameters.output_id"),
    )


def _resolve_figures_trait_program_gene_panel(config: LoadedConfig) -> ResolvedFiguresTraitProgramGenePanel:
    workflow = _workflow_mapping(config, "trait_program_gene_panel")
    inputs = _validate_mapping(workflow.get("inputs"), "workflows.figures.trait_program_gene_panel.inputs")
    outputs = _validate_mapping(workflow.get("outputs"), "workflows.figures.trait_program_gene_panel.outputs")
    parameters = _validate_mapping(workflow.get("parameters"), "workflows.figures.trait_program_gene_panel.parameters")

    mapping_entries = _load_file_id_map(config, inputs.get("file_id_map"))
    if not mapping_entries:
        raise ConfigError("workflows.figures.trait_program_gene_panel.inputs.file_id_map is required")

    posterior_dir = config.resolve_path_or_artifact(inputs.get("posterior_dir"), "genebayes", "posterior")
    posterior_name_map = _load_two_column_name_map(
        config,
        inputs.get("posterior_name_map"),
        label="workflows.figures.trait_program_gene_panel.inputs.posterior_name_map",
    )
    lof_ids = _string_list(parameters.get("lof_ids"), "workflows.figures.trait_program_gene_panel.parameters.lof_ids")
    if not lof_ids:
        lof_ids = tuple(entry.id2 for entry in mapping_entries)

    targets: list[GeneLevelScatterTarget] = []
    for source_id in lof_ids:
        entry = next((item for item in mapping_entries if item.id2 == source_id), None)
        if entry is None:
            raise ConfigError(
                f"workflows.figures.trait_program_gene_panel.parameters.lof_ids value '{source_id}' "
                "was not found in file_id_map column id2"
            )
        trait_stem = _file_label(str(entry.path2))
        targets.append(
            GeneLevelScatterTarget(
                source_id=_normalize_output_id(source_id, "workflows.figures.trait_program_gene_panel.parameters.lof_ids"),
                trait_stem=trait_stem,
                posterior_path=posterior_dir / _posterior_filename_for_entry(entry, posterior_name_map),
            )
        )

    k = int(parameters.get("k", 60))
    return ResolvedFiguresTraitProgramGenePanel(
        config=config,
        program_association_dir=config.resolve_path(inputs.get("program_association_dir"))
        or config.project_root / "outputs" / "perturbseq" / "cnmf_genomewide" / "trait_association" / "K562GW" / "ProgramLevel",
        regulation_dir=config.resolve_path(inputs.get("regulation_dir"))
        or config.project_root / "outputs" / "perturbseq" / "cnmf_genomewide" / "cNMF_regulation" / "K562GW",
        spectra_path=config.resolve_path(inputs.get("spectra_path"))
        or config.project_root / "outputs" / "perturbseq" / "cnmf_genomewide" / "cNMF" / "cNMF_all" / f"cNMF_all.gene_spectra_score.k_{k}.dt_0_5.txt",
        gene_map=config.resolve_path(inputs.get("gene_map"))
        or config.project_root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID",
        output_dir=config.resolve_path_or_artifact(outputs.get("output_dir"), "trait_program_gene_panel"),
        targets=tuple(targets),
        k=k,
        max_programs=int(parameters.get("max_programs", 8)),
        max_genes_per_side=int(parameters.get("max_genes_per_side", 8)),
        hit_abs_gamma_threshold=float(parameters.get("hit_abs_gamma_threshold", 0.1)),
        loading_top_n=int(parameters.get("loading_top_n", 200)),
        regulator_fdr_threshold=float(parameters.get("regulator_fdr_threshold", 0.05)),
        min_abs_score=float(parameters.get("min_abs_score", 1.3)),
        render_plot=bool(parameters.get("render_plot", True)),
    )


def _script_path(name: str) -> Path:
    path = SCRIPT_DIR / name
    if not path.exists():
        raise ConfigError(f"Figure script is missing: {path}")
    return path


def _rscript_command(config: LoadedConfig, script_name: str, *args: str | int | float | Path) -> list[str]:
    command = config.executable_command("plot_rscript", ["Rscript"])
    command.append(str(_script_path(script_name)))
    command.extend(str(arg) for arg in args)
    return command


def _trait_stem(trait_file: str) -> str:
    stem = trait_file.removesuffix(".per_gene_estimates.tsv")
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in stem)


def _file_label(path_like: str) -> str:
    stem = Path(path_like).name
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
    return _normalize_output_id(stem, "output_id")


def _ensure_directory_task(path: Path, label: str) -> Task:
    return Task(
        name=f"Create {label} output directory",
        preview=f"Create output directory: {path}",
        action=lambda: path.mkdir(parents=True, exist_ok=True),
    )


def _write_manifest_task(path: Path, rows: list[dict[str, str]], label: str) -> Task:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    def _action() -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    return Task(
        name=f"Write {label} manifest",
        preview=f"Write manifest: {path}",
        action=_action,
    )


def build_figures_cnmf_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_cnmf(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir / "program_regulator", "figures-cnmf tables"),
        _ensure_directory_task(plots_dir / "program_regulator", "figures-cnmf plots"),
        _ensure_directory_task(tables_dir / "corregulation", "figures-cnmf tables"),
        _ensure_directory_task(plots_dir / "corregulation", "figures-cnmf plots"),
        _ensure_directory_task(meta_dir, "figures-cnmf meta"),
    ]

    label_arg = ",".join(str(program) for program in resolved.plot_label_programs)
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.trait_targets:
        table_path = tables_dir / "program_regulator" / f"{target.output_id}.tsv"
        plot_prefix = plots_dir / "program_regulator" / target.output_id
        command = _rscript_command(
            config,
            "cnmf_program_regulator_scatter.R",
            resolved.program_association_dir,
            target.trait_file,
            resolved.k,
            label_arg,
            table_path,
            plot_prefix,
        )
        tasks.append(
            Task(
                name=f"Render figures-cnmf program/regulator scatter for {target.output_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.output_id,
                "figure_kind": "cnmf_program_regulator_scatter",
                "source_id": target.output_id,
                "source_ref": target.trait_file,
                "table_path": str(table_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
            }
        )

    for target in resolved.corregulation_targets:
        table_path = tables_dir / "corregulation" / f"{target.output_id}.tsv"
        plot_prefix = plots_dir / "corregulation" / target.output_id
        command = _rscript_command(
            config,
            "cnmf_corregulation_scatter.R",
            resolved.cnmf_regulation_dir,
            resolved.k,
            target.program_a,
            target.program_b,
            table_path,
            plot_prefix,
        )
        tasks.append(
            Task(
                name=f"Render figures-cnmf corregulation scatter for {target.output_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.output_id,
                "figure_kind": "cnmf_corregulation_scatter",
                "source_id": target.output_id,
                "source_ref": f"{target.program_a} vs {target.program_b}",
                "table_path": str(table_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-cnmf"))
    return tasks


def build_figures_burden_volcano_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_burden_volcano(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-burden-volcano tables"),
        _ensure_directory_task(plots_dir, "figures-burden-volcano plots"),
        _ensure_directory_task(meta_dir, "figures-burden-volcano meta"),
    ]
    geneset_arg = ",".join(resolved.highlight_genesets)
    data_geneset_arg = ",".join(resolved.data_genesets)
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.targets:
        genes_path = tables_dir / f"{target.source_id}_genes.tsv"
        hits_path = tables_dir / f"{target.source_id}_hits.tsv"
        plot_prefix = plots_dir / target.source_id
        command = _rscript_command(
            config,
            "burden_volcano.R",
            target.source_path,
            resolved.gene_map,
            resolved.geneset_dir,
            geneset_arg,
            genes_path,
            plot_prefix,
            resolved.label_fdr_threshold,
            resolved.line_fdr_threshold,
            data_geneset_arg,
            resolved.spectra_path,
            resolved.k,
            resolved.top_n_program_genes,
        )
        tasks.append(
            Task(
                name=f"Render figures-burden-volcano for {target.source_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.source_id,
                "figure_kind": "burden_volcano",
                "source_id": target.source_id,
                "source_path": str(target.source_path),
                "table_path": str(genes_path),
                "genes_path": str(genes_path),
                "hits_path": str(hits_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
                "plot_png": str(plot_prefix.with_suffix(".png")),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-burden-volcano"))
    return tasks


def build_figures_posterior_volcano_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_posterior_volcano(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-posterior-volcano tables"),
        _ensure_directory_task(plots_dir, "figures-posterior-volcano plots"),
        _ensure_directory_task(meta_dir, "figures-posterior-volcano meta"),
    ]
    geneset_arg = ",".join(resolved.highlight_genesets)
    data_geneset_arg = ",".join(resolved.data_genesets)
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.targets:
        genes_path = tables_dir / f"{target.source_id}_genes.tsv"
        hits_path = tables_dir / f"{target.source_id}_hits.tsv"
        plot_prefix = plots_dir / target.source_id
        command = _rscript_command(
            config,
            "posterior_volcano.R",
            target.posterior_path,
            resolved.gene_map,
            resolved.geneset_dir,
            geneset_arg,
            genes_path,
            plot_prefix,
            resolved.label_fdr_threshold,
            resolved.line_fdr_threshold,
            data_geneset_arg,
            resolved.spectra_path,
            resolved.k,
            resolved.top_n_program_genes,
        )
        tasks.append(
            Task(
                name=f"Render figures-posterior-volcano for {target.source_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.source_id,
                "figure_kind": "posterior_volcano",
                "source_id": target.source_id,
                "posterior_path": str(target.posterior_path),
                "table_path": str(genes_path),
                "genes_path": str(genes_path),
                "hits_path": str(hits_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
                "plot_png": str(plot_prefix.with_suffix(".png")),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-posterior-volcano"))
    return tasks


def build_figures_program_heatmap_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_program_heatmap(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-program-heatmap tables"),
        _ensure_directory_task(plots_dir, "figures-program-heatmap plots"),
        _ensure_directory_task(meta_dir, "figures-program-heatmap meta"),
    ]

    trait_rows = [
        {"trait_id": target.output_id, "trait_file": target.trait_file}
        for target in resolved.trait_targets
    ]
    trait_target_path = meta_dir / "trait_targets.tsv"
    tasks.append(_write_manifest_task(trait_target_path, trait_rows, "figures-program-heatmap trait target"))

    table_prefix = tables_dir / resolved.output_id
    plot_prefix = plots_dir / resolved.output_id
    metrics_arg = ",".join(resolved.metrics)
    command = _rscript_command(
        config,
        "program_heatmap.R",
        resolved.program_association_dir,
        resolved.k,
        trait_target_path,
        table_prefix,
        plot_prefix,
        metrics_arg,
    )
    tasks.append(
        Task(
            name=f"Render figures-program-heatmap for {resolved.output_id}",
            preview=preview_command(command),
            command=command,
        )
    )

    manifest_rows = [
        {
            "figure_id": resolved.output_id,
            "figure_kind": "program_heatmap",
            "trait_target_path": str(trait_target_path),
            "table_prefix": str(table_prefix),
            "plot_prefix": str(plot_prefix),
            "metrics": ",".join(resolved.metrics),
        }
    ]
    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-program-heatmap"))
    return tasks


def build_figures_gwas_manhattan_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_gwas_manhattan(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-gwas-manhattan tables"),
        _ensure_directory_task(plots_dir, "figures-gwas-manhattan plots"),
        _ensure_directory_task(meta_dir, "figures-gwas-manhattan meta"),
    ]
    geneset_arg = ",".join(resolved.highlight_genesets)
    data_geneset_arg = ",".join(resolved.data_genesets)
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.targets:
        variants_path = tables_dir / f"{target.source_id}_variants.tsv"
        hits_path = tables_dir / f"{target.source_id}_hits.tsv"
        plot_prefix = plots_dir / target.source_id
        command = _rscript_command(
            config,
            "gwas_manhattan.R",
            target.source_path,
            resolved.gene_annotation,
            resolved.geneset_dir,
            geneset_arg,
            variants_path,
            plot_prefix,
            resolved.flank_bp,
            resolved.label_p_threshold,
            resolved.genomewide_threshold,
            data_geneset_arg,
            resolved.spectra_path,
            resolved.gene_map,
            resolved.k,
            resolved.top_n_program_genes,
        )
        tasks.append(
            Task(
                name=f"Render figures-gwas-manhattan for {target.source_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.source_id,
                "figure_kind": "gwas_manhattan",
                "source_id": target.source_id,
                "source_path": str(target.source_path),
                "table_path": str(variants_path),
                "variants_path": str(variants_path),
                "hits_path": str(hits_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
                "plot_png": str(plot_prefix.with_suffix(".png")),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-gwas-manhattan"))
    return tasks


def build_figures_gene_level_scatter_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_gene_level_scatter(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-gene-level-scatter tables"),
        _ensure_directory_task(plots_dir, "figures-gene-level-scatter plots"),
        _ensure_directory_task(meta_dir, "figures-gene-level-scatter meta"),
    ]
    highlight_arg = ",".join(resolved.highlight_genes)
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.targets:
        table_path = tables_dir / f"{target.source_id}.tsv"
        plot_prefix = plots_dir / target.source_id
        command = _rscript_command(
            config,
            "gene_level_scatter.R",
            target.posterior_path,
            resolved.limma_path,
            resolved.shet_path,
            target.source_id,
            table_path,
            plot_prefix,
            resolved.top_n_labels,
            highlight_arg,
            resolved.y_limit,
            resolved.gene_map,
        )
        tasks.append(
            Task(
                name=f"Render figures-gene-level-scatter for {target.source_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.source_id,
                "figure_kind": "gene_level_scatter",
                "trait_stem": target.trait_stem,
                "posterior_path": str(target.posterior_path),
                "limma_path": str(resolved.limma_path),
                "shet_path": str(resolved.shet_path),
                "table_path": str(table_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
                "plot_png": str(plot_prefix.with_suffix(".png")),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-gene-level-scatter"))
    return tasks


def build_figures_cross_trait_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_cross_trait(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-cross-trait tables"),
        _ensure_directory_task(plots_dir, "figures-cross-trait plots"),
        _ensure_directory_task(meta_dir, "figures-cross-trait meta"),
    ]
    highlight_arg = ",".join(resolved.highlight_genes)
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.targets:
        table_path = tables_dir / f"{target.output_id}.tsv"
        plot_prefix = plots_dir / target.output_id
        command = _rscript_command(
            config,
            "cross_trait_scatter.R",
            target.posterior_path_x,
            target.posterior_path_y,
            target.x_label,
            target.y_label,
            resolved.gene_map,
            table_path,
            plot_prefix,
            resolved.top_n_labels,
            highlight_arg,
        )
        tasks.append(
            Task(
                name=f"Render figures-cross-trait for {target.output_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.output_id,
                "figure_kind": "cross_trait",
                "x_id": target.x_id,
                "y_id": target.y_id,
                "x_label": target.x_label,
                "y_label": target.y_label,
                "posterior_path_x": str(target.posterior_path_x),
                "posterior_path_y": str(target.posterior_path_y),
                "table_path": str(table_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
                "plot_png": str(plot_prefix.with_suffix(".png")),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-cross-trait"))
    return tasks


def build_figures_gene_level_qq_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_gene_level_qq(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-gene-level-qq tables"),
        _ensure_directory_task(plots_dir, "figures-gene-level-qq plots"),
        _ensure_directory_task(meta_dir, "figures-gene-level-qq meta"),
    ]
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.targets:
        table_path = tables_dir / f"{target.source_id}.tsv"
        plot_prefix = plots_dir / target.source_id
        command = _rscript_command(
            resolved.config,
            "gene_level_qq.R",
            target.posterior_path,
            resolved.limma_path,
            resolved.shet_path,
            resolved.gene_map,
            target.source_id,
            table_path,
            plot_prefix,
            resolved.y_limit,
            int(resolved.render_plot),
        )
        tasks.append(
            Task(
                name=f"Render figures-gene-level-qq for {target.source_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.source_id,
                "figure_kind": "gene_level_qq",
                "trait_stem": target.trait_stem,
                "posterior_path": str(target.posterior_path),
                "limma_path": str(resolved.limma_path),
                "shet_path": str(resolved.shet_path),
                "table_path": str(table_path),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")) if resolved.render_plot else "",
                "plot_png": str(plot_prefix.with_suffix(".png")) if resolved.render_plot else "",
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-gene-level-qq"))
    return tasks


def build_figures_program_rankings_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_program_rankings(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-program-rankings tables"),
        _ensure_directory_task(plots_dir, "figures-program-rankings plots"),
        _ensure_directory_task(meta_dir, "figures-program-rankings meta"),
    ]
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.trait_targets:
        table_prefix = tables_dir / target.output_id
        plot_prefix = plots_dir / target.output_id
        command = _rscript_command(
            resolved.config,
            "program_rankings.R",
            resolved.program_association_dir,
            target.trait_file,
            target.output_id,
            resolved.k,
            resolved.top_n,
            table_prefix,
            plot_prefix,
        )
        tasks.append(
            Task(
                name=f"Render figures-program-rankings for {target.output_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.output_id,
                "figure_kind": "program_rankings",
                "source_ref": target.trait_file,
                "table_prefix": str(table_prefix),
                "plot_program_pdf": str(plots_dir / f"{target.output_id}_program_score.pdf"),
                "plot_regulator_pdf": str(plots_dir / f"{target.output_id}_regulator_score.pdf"),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-program-rankings"))
    return tasks


def build_figures_gwas_locus_zoom_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_gwas_locus_zoom(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-gwas-locus-zoom tables"),
        _ensure_directory_task(plots_dir, "figures-gwas-locus-zoom plots"),
        _ensure_directory_task(meta_dir, "figures-gwas-locus-zoom meta"),
    ]
    manifest_rows: list[dict[str, str]] = []

    for target in resolved.targets:
        table_prefix = tables_dir / target.output_id
        plot_prefix = plots_dir / target.output_id
        command = _rscript_command(
            resolved.config,
            "gwas_locus_zoom.R",
            target.source_path,
            resolved.gene_annotation,
            target.source_id,
            target.locus_label,
            target.chrom,
            target.start,
            target.end,
            target.flank_bp,
            table_prefix,
            plot_prefix,
            resolved.genomewide_threshold,
            resolved.label_top_n,
        )
        tasks.append(
            Task(
                name=f"Render figures-gwas-locus-zoom for {target.output_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.output_id,
                "figure_kind": "gwas_locus_zoom",
                "source_id": target.source_id,
                "source_path": str(target.source_path),
                "chrom": target.chrom,
                "start": str(target.start),
                "end": str(target.end),
                "table_prefix": str(table_prefix),
                "summary_path": str(Path(f"{table_prefix}_summary.tsv")),
                "variants_path": str(Path(f"{table_prefix}_variants.tsv")),
                "genes_path": str(Path(f"{table_prefix}_genes.tsv")),
                "top_hits_path": str(Path(f"{table_prefix}_top_hits.tsv")),
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
                "plot_png": str(plot_prefix.with_suffix(".png")),
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-gwas-locus-zoom"))
    return tasks


def build_figures_cross_trait_heatmap_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_cross_trait_heatmap(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-cross-trait-heatmap tables"),
        _ensure_directory_task(plots_dir, "figures-cross-trait-heatmap plots"),
        _ensure_directory_task(meta_dir, "figures-cross-trait-heatmap meta"),
    ]

    target_rows = [
        {
            "trait_id": target.output_id,
            "posterior_path": target.trait_file,
        }
        for target in resolved.targets
    ]
    target_path = meta_dir / "trait_targets.tsv"
    tasks.append(_write_manifest_task(target_path, target_rows, "figures-cross-trait-heatmap trait target"))

    table_prefix = tables_dir / resolved.output_id
    plot_prefix = plots_dir / resolved.output_id
    command = _rscript_command(
        resolved.config,
        "cross_trait_heatmap.R",
        target_path,
        resolved.gene_map,
        table_prefix,
        plot_prefix,
        resolved.method,
    )
    tasks.append(
        Task(
            name=f"Render figures-cross-trait-heatmap for {resolved.output_id}",
            preview=preview_command(command),
            command=command,
        )
    )

    manifest_rows = [
        {
            "figure_id": resolved.output_id,
            "figure_kind": "cross_trait_heatmap",
            "trait_target_path": str(target_path),
            "table_prefix": str(table_prefix),
            "plot_pdf": str(plot_prefix.with_suffix(".pdf")),
            "plot_png": str(plot_prefix.with_suffix(".png")),
            "method": resolved.method,
        }
    ]
    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-cross-trait-heatmap"))
    return tasks


def build_figures_cnmf_program_top_genes_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_cnmf_program_top_genes(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-cnmf-program-top-genes tables"),
        _ensure_directory_task(plots_dir, "figures-cnmf-program-top-genes plots"),
        _ensure_directory_task(meta_dir, "figures-cnmf-program-top-genes meta"),
    ]

    programs_arg = ",".join(str(p) for p in resolved.programs)
    table_prefix = tables_dir / resolved.output_id
    plot_prefix = plots_dir / resolved.output_id

    command = _rscript_command(
        resolved.config,
        "cnmf_program_top_genes.R",
        resolved.spectra_path,
        resolved.gene_map,
        resolved.k,
        programs_arg,
        resolved.top_n,
        table_prefix,
        plot_prefix,
    )
    tasks.append(
        Task(
            name=f"Render cNMF program top genes for {resolved.output_id}",
            preview=preview_command(command),
            command=command,
        )
    )

    manifest_rows = [
        {
            "figure_id": resolved.output_id,
            "figure_kind": "cnmf_program_top_genes",
            "table_path": f"{table_prefix}_top_genes.tsv",
            "plot_pdf": f"{plot_prefix}_top_genes.pdf",
            "plot_png": f"{plot_prefix}_top_genes.png",
        }
    ]
    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-cnmf-program-top-genes"))
    return tasks


def build_figures_cnmf_program_enrichment_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_cnmf_program_enrichment(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-cnmf-program-enrichment tables"),
        _ensure_directory_task(plots_dir, "figures-cnmf-program-enrichment plots"),
        _ensure_directory_task(meta_dir, "figures-cnmf-program-enrichment meta"),
    ]

    programs_arg = ",".join(str(p) for p in resolved.programs)
    genesets_arg = ",".join(resolved.genesets)
    table_prefix = tables_dir / resolved.output_id
    plot_prefix = plots_dir / resolved.output_id

    command = _rscript_command(
        resolved.config,
        "cnmf_program_enrichment.R",
        resolved.spectra_path,
        resolved.gene_map,
        resolved.k,
        programs_arg,
        resolved.top_n,
        resolved.geneset_dir,
        genesets_arg,
        table_prefix,
        plot_prefix,
    )
    tasks.append(
        Task(
            name=f"Render cNMF program enrichment for {resolved.output_id}",
            preview=preview_command(command),
            command=command,
        )
    )

    manifest_rows = [
        {
            "figure_id": resolved.output_id,
            "figure_kind": "cnmf_program_enrichment",
            "table_path": f"{table_prefix}_enrichment.tsv",
            "plot_pdf": f"{plot_prefix}_enrichment.pdf",
            "plot_png": f"{plot_prefix}_enrichment.png",
        }
    ]
    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-cnmf-program-enrichment"))
    return tasks


def build_figures_trait_program_gene_panel_tasks(config: LoadedConfig) -> list[Task]:
    resolved = _resolve_figures_trait_program_gene_panel(config)
    tables_dir = resolved.output_dir / "tables"
    plots_dir = resolved.output_dir / "plots"
    meta_dir = resolved.output_dir / "meta"
    tasks: list[Task] = [
        _ensure_directory_task(tables_dir, "figures-trait-program-gene-panel tables"),
        _ensure_directory_task(plots_dir, "figures-trait-program-gene-panel plots"),
        _ensure_directory_task(meta_dir, "figures-trait-program-gene-panel meta"),
    ]

    manifest_rows: list[dict[str, str]] = []
    for target in resolved.targets:
        table_prefix = tables_dir / target.source_id
        plot_prefix = plots_dir / target.source_id
        command = _rscript_command(
            resolved.config,
            "trait_program_gene_panel.R",
            resolved.program_association_dir,
            resolved.regulation_dir,
            resolved.spectra_path,
            resolved.gene_map,
            target.posterior_path,
            target.source_id,
            resolved.k,
            table_prefix,
            plot_prefix,
            resolved.max_programs,
            resolved.max_genes_per_side,
            resolved.hit_abs_gamma_threshold,
            resolved.loading_top_n,
            resolved.regulator_fdr_threshold,
            resolved.min_abs_score,
            int(resolved.render_plot),
        )
        tasks.append(
            Task(
                name=f"Render trait-program-gene panel for {target.source_id}",
                preview=preview_command(command),
                command=command,
            )
        )
        manifest_rows.append(
            {
                "figure_id": target.source_id,
                "figure_kind": "trait_program_gene_panel",
                "source_id": target.source_id,
                "trait_stem": target.trait_stem,
                "table_long_path": f"{table_prefix}_long.tsv",
                "table_program_path": f"{table_prefix}_programs.tsv",
                "plot_pdf": str(plot_prefix.with_suffix(".pdf")) if resolved.render_plot else "",
                "plot_png": str(plot_prefix.with_suffix(".png")) if resolved.render_plot else "",
            }
        )

    tasks.append(_write_manifest_task(meta_dir / "manifest.tsv", manifest_rows, "figures-trait-program-gene-panel"))
    return tasks


def _consolidate_manifests_task(config: LoadedConfig, figure_kinds: list[str]) -> Task:
    """Produce a single ``manifest_all.tsv`` that merges every sub-workflow manifest."""
    figures_root = config.resolve_path_or_artifact(None)

    def _action() -> None:
        all_rows: list[dict[str, str]] = []
        all_fieldnames: list[str] = ["figure_kind"]
        for kind in figure_kinds:
            manifest_path = figures_root / kind / "meta" / "manifest.tsv"
            if not manifest_path.exists():
                continue
            with manifest_path.open("r", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    row.setdefault("figure_kind", kind)
                    for col in row:
                        if col not in all_fieldnames:
                            all_fieldnames.append(col)
                    all_rows.append(row)

        output_path = figures_root / "manifest_all.tsv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=all_fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in all_rows:
                writer.writerow(row)

    return Task(
        name="Consolidate all figures manifests",
        preview=f"Merge manifests from {len(figure_kinds)} figure sub-workflows into manifest_all.tsv",
        action=_action,
    )


def build_figures_tasks(config: LoadedConfig) -> list[Task]:
    figures = _workflow_root(config)
    tasks: list[Task] = []
    active_kinds: list[str] = []
    if figures.get("cnmf"):
        tasks.extend(build_figures_cnmf_tasks(config))
        active_kinds.append("cnmf")
    if figures.get("program_heatmap"):
        tasks.extend(build_figures_program_heatmap_tasks(config))
        active_kinds.append("program_heatmap")
    if figures.get("burden_volcano"):
        tasks.extend(build_figures_burden_volcano_tasks(config))
        active_kinds.append("burden_volcano")
    if figures.get("posterior_volcano"):
        tasks.extend(build_figures_posterior_volcano_tasks(config))
        active_kinds.append("posterior_volcano")
    if figures.get("gwas_manhattan"):
        tasks.extend(build_figures_gwas_manhattan_tasks(config))
        active_kinds.append("gwas_manhattan")
    if figures.get("gene_level_scatter"):
        tasks.extend(build_figures_gene_level_scatter_tasks(config))
        active_kinds.append("gene_level_scatter")
    if figures.get("cross_trait"):
        tasks.extend(build_figures_cross_trait_tasks(config))
        active_kinds.append("cross_trait")
    if figures.get("gene_level_qq"):
        tasks.extend(build_figures_gene_level_qq_tasks(config))
        active_kinds.append("gene_level_qq")
    if figures.get("program_rankings"):
        tasks.extend(build_figures_program_rankings_tasks(config))
        active_kinds.append("program_rankings")
    if figures.get("gwas_locus_zoom"):
        tasks.extend(build_figures_gwas_locus_zoom_tasks(config))
        active_kinds.append("gwas_locus_zoom")
    if figures.get("cross_trait_heatmap"):
        tasks.extend(build_figures_cross_trait_heatmap_tasks(config))
        active_kinds.append("cross_trait_heatmap")
    if figures.get("cnmf_program_top_genes"):
        tasks.extend(build_figures_cnmf_program_top_genes_tasks(config))
        active_kinds.append("cnmf_program_top_genes")
    if figures.get("cnmf_program_enrichment"):
        tasks.extend(build_figures_cnmf_program_enrichment_tasks(config))
        active_kinds.append("cnmf_program_enrichment")
    if figures.get("trait_program_gene_panel"):
        tasks.extend(build_figures_trait_program_gene_panel_tasks(config))
        active_kinds.append("trait_program_gene_panel")
    if not tasks:
        raise ConfigError("workflows.figures must contain at least one configured subworkflow")

    # Append a final task that consolidates all sub-workflow manifests.
    tasks.append(_consolidate_manifests_task(config, active_kinds))
    return tasks


__all__ = [
    "build_figures_burden_volcano_tasks",
    "build_figures_posterior_volcano_tasks",
    "build_figures_cnmf_tasks",
    "build_figures_cross_trait_tasks",
    "build_figures_cross_trait_heatmap_tasks",
    "build_figures_gene_level_qq_tasks",
    "build_figures_gene_level_scatter_tasks",
    "build_figures_gwas_locus_zoom_tasks",
    "build_figures_gwas_manhattan_tasks",
    "build_figures_program_rankings_tasks",
    "build_figures_program_heatmap_tasks",
    "build_figures_cnmf_program_top_genes_tasks",
    "build_figures_cnmf_program_enrichment_tasks",
    "build_figures_trait_program_gene_panel_tasks",
    "build_figures_tasks",
]

