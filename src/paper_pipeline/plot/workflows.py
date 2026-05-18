from __future__ import annotations

import csv
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from ..dataprep.config import ConfigError, LoadedConfig
from ..dataprep.tasks import Task
from .legacy import (
    ensure_directories,
    format_command_preview,
    link_alias,
    link_output_path,
    link_path,
    populate_directory,
    replace_tree,
    run_command,
    run_commands,
    write_text,
)
from .resources import script_path as plot_script_path


ScriptAction = Callable[[], None]


@dataclass(frozen=True)
class PlotStageSpec:
    key: str
    name: str
    summary: str
    scripts: tuple[str, ...]
    upstream_products: tuple[str, ...]
    expected_outputs: tuple[str, ...]
    r_packages: tuple[str, ...] = ()
    external_commands: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class PlotWorkflowSpec:
    key: str
    cli_name: str
    title: str
    summary: str
    input_descriptions: tuple[tuple[str, str], ...]
    upstream_workflows: tuple[str, ...]
    stages: tuple[PlotStageSpec, ...]
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class ResolvedPlotWorkflow:
    config: LoadedConfig
    spec: PlotWorkflowSpec
    inputs: dict[str, Path]
    input_rows: tuple[tuple[str, str, Path], ...]
    output_dir: Path
    parameters: dict[str, Any]

    def input_path(self, key: str) -> Path:
        return self.inputs[key]

    def parameter(self, key: str, default: Any = None) -> Any:
        return self.parameters.get(key, default)


@dataclass(frozen=True)
class PlannedStage:
    spec: PlotStageSpec
    preview_lines: tuple[str, ...]
    action: ScriptAction


GENE_BURDEN_WORKFLOW = PlotWorkflowSpec(
    key="gene_burden",
    cli_name="plot-gene-burden",
    title="plot-gene-burden",
    summary=(
        "Figure2-centric plotting workflow that executes legacy GeneBayes/GWAS figure "
        "scripts in a workspace-compatible layout."
    ),
    input_descriptions=(
        ("genebayes_posterior_dir", "GeneBayes posterior per-gene estimates."),
        ("raw_burden_dir", "Raw burden summary statistics for volcano and enrichment scripts."),
        ("gwas_dir", "GWAS summary statistics, closest-gene references, and gene annotations."),
        ("geneset_dir", "Gene-set collections highlighted in Figure2."),
        ("gene_map", "Gene symbol to Ensembl mapping file."),
        ("trait_correspondence", "Backman/Niele trait correspondence table."),
        ("ld_reference_dir", "1000G LD reference used by GWAS clumping."),
    ),
    upstream_workflows=("genebayes",),
    stages=(
        PlotStageSpec(
            key="gwas_closest_genes",
            name="GWAS clumping and closest-gene preparation",
            summary="Run Figure2 shell scripts that clump GWAS loci and annotate closest genes.",
            scripts=(
                "Figure2/3_GWAS_clumping.sh",
                "Figure2/summarize_clumped_variants.R",
                "Figure2/4_annotate_closest_genes.sh",
            ),
            upstream_products=("External GWAS summary statistics and LD reference panels.",),
            expected_outputs=(
                "clump/*.assoc",
                "clump/*.clumped",
                "clump/*.top.100kbpmerged",
                "data/GWAS/closest_genes/*_genes_closest.txt",
            ),
            external_commands=("bash", "plink", "bedtools", "Rscript"),
        ),
        PlotStageSpec(
            key="gene_burden_figures",
            name="Figure2 plotting",
            summary="Render the main Figure2 Manhattan, volcano, enrichment, and cross-trait panels.",
            scripts=(
                "Figure2/1_manhattanplot.R",
                "Figure2/2_volcanoplot.R",
                "Figure2/5_pathway_enrichment_GWAS.R",
                "Figure2/6_pathway_enrichment_LoF.R",
                "Figure2/7_pathway_enrichment_plot.R",
                "Figure2/8_CrossTrait_plot.R",
            ),
            upstream_products=(
                "GeneBayes posterior outputs from `genebayes`.",
                "Closest-gene files from the clumping stage.",
            ),
            expected_outputs=("Fig2A.png", "Fig2A.pdf", "Fig2B.pdf", "Fig2C.pdf", "Fig2D.pdf"),
            r_packages=(
                "clusterProfiler",
                "data.table",
                "dplyr",
                "ggbreak",
                "ggrastr",
                "ggplot2",
                "ggrepel",
                "ggsci",
                "msigdbr",
                "org.Hs.eg.db",
                "reshape2",
                "scattermore",
                "topr",
            ),
        ),
    ),
)


PERTURBSEQ_GENE_LEVEL_WORKFLOW = PlotWorkflowSpec(
    key="perturbseq_gene_level",
    cli_name="plot-perturbseq-gene-level",
    title="plot-perturbseq-gene-level",
    summary=(
        "Figure3-centric plotting workflow for limma gene-level summaries combined with "
        "GeneBayes posterior and GWAS context."
    ),
    input_descriptions=(
        ("gene_level_summary_dir", "Summarized limma outputs from `perturbseq-gene-level-summarize`."),
        ("genebayes_posterior_dir", "GeneBayes posterior per-gene estimates."),
        ("raw_burden_dir", "Raw burden summary statistics."),
        ("gwas_dir", "GWAS closest-gene files and annotations."),
        ("gene_map", "Gene symbol to Ensembl mapping file."),
        ("shet_path", "s_het covariate table."),
        ("burden_reg_correlation_dir", "Legacy QQ-plot burden/regulation correlation inputs."),
        ("gene_level_trait_metadata", "Legacy Backman/Niele trait metadata for Figure3 QQ plots."),
    ),
    upstream_workflows=("perturbseq-gene-level", "genebayes"),
    stages=(
        PlotStageSpec(
            key="gene_level_trait_association",
            name="Gene-level burden/regulation association",
            summary="Run Figure3 association scripts that correlate limma effects with GeneBayes traits.",
            scripts=("Figure3/1_gamma_beta_cor.R",),
            upstream_products=("Perturbseq limma summaries.", "GeneBayes posterior outputs."),
            expected_outputs=(
                "data/Perturbseq/trait_association/K562GW/GeneLevel/*_geneRegulation_correlation.txt",
            ),
            r_packages=("boot", "data.table"),
        ),
        PlotStageSpec(
            key="gene_level_enrichment",
            name="GWAS and LoF enrichment",
            summary="Run Figure3 enrichment script around target gene regulators.",
            scripts=("Figure3/2_GWAS_enrich.R",),
            upstream_products=("Gene-level limma summaries.", "GWAS closest genes.", "Raw burden statistics."),
            expected_outputs=("MCH_GWAS_HBA1regulator_enrichment.txt", "MCH_LoF_HBA1regulator_enrichment.txt", "Fig3B.pdf"),
            r_packages=("boot", "data.table", "ggplot2"),
        ),
        PlotStageSpec(
            key="gene_level_rendering",
            name="Figure3 rendering",
            summary="Render the Figure3 scatter and QQ plots.",
            scripts=("Figure3/3_plot_gamma_vs_cor.R", "Figure3/4_plot_QQ.R"),
            upstream_products=("Gene-level association tables.", "Legacy burden/regulation QQ inputs."),
            expected_outputs=("Fig3c.pdf", "Fig3D.pdf"),
            r_packages=("ggplot2", "ggrepel", "ggsci", "RColorBrewer"),
        ),
    ),
)


PERTURBSEQ_CNMF_WORKFLOW = PlotWorkflowSpec(
    key="perturbseq_cnmf",
    cli_name="plot-perturbseq-cnmf",
    title="plot-perturbseq-cnmf",
    summary=(
        "Figure4-centric plotting workflow for cNMF program-level outputs, regulator "
        "association summaries, and trans-eQTL follow-up."
    ),
    input_descriptions=(
        ("cnmf_root", "cNMF consensus outputs produced by `perturbseq-cnmf-genomewide`."),
        ("cnmf_regulation_dir", "Program-level perturbation effect summaries."),
        ("program_association_dir", "Existing program/regulator trait association outputs."),
        ("genebayes_posterior_dir", "GeneBayes posterior per-gene estimates."),
        ("gwas_dir", "GWAS summary statistics."),
        ("gene_map", "Gene symbol to Ensembl mapping file."),
        ("shet_path", "s_het covariate table."),
        ("metadata", "Perturbseq metadata table."),
        ("trans_eqtl_path", "External trans-eQTL summary table."),
    ),
    upstream_workflows=("perturbseq-cnmf-genomewide", "genebayes"),
    stages=(
        PlotStageSpec(
            key="program_trait_association",
            name="Program/regulator burden association",
            summary="Run the Figure4 program/regulator burden association script.",
            scripts=("Figure4/1_burden_program_regulators_test.R",),
            upstream_products=("cNMF spectra and regulatory effects.", "GeneBayes posterior outputs."),
            expected_outputs=(
                "data/Perturbseq/trait_association/K562GW/ProgramLevel/regulators_enrichment_K*_*.tsv",
                "data/Perturbseq/trait_association/K562GW/ProgramLevel/programs_enrichment_K*_*.tsv",
            ),
        ),
        PlotStageSpec(
            key="program_level_rendering",
            name="Figure4 program-level rendering",
            summary="Render Figure4 program/regulator and regression plots.",
            scripts=(
                "Figure4/2_plot_burden_program_regulators.R",
                "Figure4/3_corregulation_plot.R",
                "Figure4/4_multipleRegression.R",
            ),
            upstream_products=("Program/regulator enrichment tables.", "cNMF regulation outputs."),
            expected_outputs=("*_Fig4C.pdf", "*_Fig5.pdf", "Fig5F.pdf"),
            r_packages=("ggallin", "ggplot2", "ggrepel"),
        ),
        PlotStageSpec(
            key="trans_eqtl_follow_up",
            name="trans-eQTL follow-up",
            summary="Prepare trans-eQTL gene sets, align alleles, and render the overlay plot.",
            scripts=(
                "Figure4/5_1_transeQTL_prepare.R",
                "Figure4/5_2_transeQTL_GWASalleles.R",
                "Figure4/5_3_plot_transeQTL_vs_Perturbmodel.R",
            ),
            upstream_products=("cNMF spectra outputs.", "Program/regulator enrichment tables.", "trans-eQTL and GWAS summary statistics."),
            expected_outputs=(
                "program_trans/P*_top100genes_trans.txt",
                "program_trans/result/P*_transeQTL_*_directionalTest.txt",
                "*_transeQTLvsLoFburden_Fig4j.pdf",
            ),
            r_packages=("data.table", "ggplot2", "ggrepel", "plyr"),
        ),
    ),
)


SUPPLEMENTARY_WORKFLOW = PlotWorkflowSpec(
    key="supplementary",
    cli_name="plot-supplementary",
    title="plot-supplementary",
    summary=(
        "Supplementary plotting workflow that executes legacy Figure1, Figure5, EDFig4/5, "
        "and multi-cell supplementary scripts."
    ),
    input_descriptions=(
        ("ldsc_dir", "Legacy LDSC working directory with bed files and plotting metadata."),
        ("ldsc_ld_reference_dir", "LD reference directory used by LDSC annotation and rg scripts."),
        ("ldsc_baseline_dir", "LDSC baseline directory containing `baseline.*` resources and chromosome `.snp` files."),
        ("ldsc_weights_dir", "LDSC weights reference prefix."),
        ("ldsc_frq_dir", "LDSC allele-frequency reference prefix."),
        ("ldsc_sumstats_dir", "Directory containing LDSC-formatted GWAS sumstats."),
        ("chip_dir", "TF ChIP input directory with metadata and existing downloads."),
        ("gene_level_summary_dir", "Perturbseq limma summary outputs."),
        ("cnmf_root", "K562GW cNMF consensus outputs."),
        ("cnmf_regulation_dir", "K562GW cNMF regulation outputs."),
        ("program_association_dir", "K562GW program-level trait association outputs."),
        ("multicell_cnmf_root", "Multi-cell cNMF root used by supplementary cross-cell scripts."),
        ("multicell_cnmf_regulation_dir", "Multi-cell cNMF regulation root."),
        ("multicell_trait_association_dir", "Multi-cell trait-association root."),
        ("genebayes_posterior_dir", "GeneBayes posterior per-gene estimates."),
        ("gwas_dir", "GWAS summary statistics and closest-gene files."),
        ("growth_screening_csv", "Legacy growth-screening comparison CSV."),
        ("gene_map", "Gene symbol to Ensembl mapping file."),
        ("shet_path", "s_het covariate table."),
        ("metadata", "Perturbseq metadata CSV."),
        ("geneset_dir", "Gene-set directory used by supplementary figure scripts."),
        ("liftOver_chain", "Chain file for hg38->hg19 liftOver in the ChIP collection step."),
    ),
    upstream_workflows=("genebayes", "perturbseq-gene-level", "perturbseq-cnmf-genomewide"),
    stages=(
        PlotStageSpec(
            key="ldsc_supplementary",
            name="Figure1 LDSC supplementary workflow",
            summary="Run the legacy Figure1 LDSC shell chain and the downstream plotting script.",
            scripts=(
                "Figure1/1_annotation_LDSC.sh",
                "Figure1/2_stratify_LDSC.sh",
                "Figure1/3_summarize_stratify_results.sh",
                "Figure1/4_geneticCorrelation.sh",
                "Figure1/5_makePlot.R",
            ),
            upstream_products=("LDSC reference resources and GWAS sumstats.",),
            expected_outputs=("data/LDSC/figure/Fig1B.pdf", "data/LDSC/figure/Fig1E.pdf", "data/LDSC/figure/FigS1A.pdf"),
            external_commands=("bash", "bedtools", "ldsc.py"),
            r_packages=("ComplexHeatmap", "circlize", "ggplot2", "ggrepel", "RColorBrewer", "reshape2"),
        ),
        PlotStageSpec(
            key="chip_supplementary",
            name="EDFig4 ChIP supplementary workflow",
            summary="Run the supplementary ChIP collection, chip-score, enrichment, and heatmap scripts.",
            scripts=(
                "supplementary_figures/EDFig4_ChIP/1_collect_ChIP_K562.sh",
                "supplementary_figures/EDFig4_ChIP/2_calculate_TF_chipScore.R",
                "supplementary_figures/EDFig4_ChIP/3_ChIPscore_program_enrichment.R",
                "supplementary_figures/EDFig4_ChIP/4_select_useful_experiments.R",
                "supplementary_figures/EDFig4_ChIP/5_program_enrichment_test.R",
                "supplementary_figures/EDFig4_ChIP/6_plot_program_TF_enrichment.R",
            ),
            upstream_products=("Perturbseq gene-level summaries.", "cNMF spectra and regulation outputs.", "TF ChIP metadata/files."),
            expected_outputs=("data/TF_ChIP/chip_score/*.txt", "data/TF_ChIP/cNMF_enrichment/*.txt", "ChIP_score_enrich_program_FigS5.pdf"),
            external_commands=("bash", "bedtools", "curl", "liftOver"),
            r_packages=("ComplexHeatmap", "circlize", "data.table", "dplyr", "plyr", "RColorBrewer", "reshape2"),
        ),
        PlotStageSpec(
            key="edfig5_supplementary",
            name="EDFig5 supplementary workflow",
            summary="Run the supplementary growth-screening, causal-inference, cell-cycle, and autophagy scripts.",
            scripts=(
                "supplementary_figures/EDFig5/1_growth_screening_comparison.R",
                "supplementary_figures/EDFig5/2_causal_inference_regulators.R",
                "supplementary_figures/EDFig5/3_cellCycledistribution.R",
                "supplementary_figures/EDFig5/4_autophagy_qqlot.R",
            ),
            upstream_products=("cNMF spectra/regulation outputs.", "GeneBayes posterior outputs.", "Growth-screening CSV."),
            expected_outputs=("data/fitness_compare/*.txt", "data/programs_causal_relation_summary.txt", "Fig5I.pdf", "*_Fig5J.pdf"),
            r_packages=("cowplot", "data.table", "DescTools", "dplyr", "ggplot2", "ggrepel"),
        ),
        PlotStageSpec(
            key="multi_cell_supplementary",
            name="Multi-cell supplementary workflow",
            summary="Run the FigureS4-12 multi-cell correlation, overlap, regression, and conditional plotting scripts.",
            scripts=(
                "supplementary_figures/FigureS4-12_multiCell_model/1_regulatory_effects_correlation.R",
                "supplementary_figures/FigureS4-12_multiCell_model/2_overlap_program_genes.R",
                "supplementary_figures/FigureS4-12_multiCell_model/3_make_overview_heatmap_burdenCor.R",
                "supplementary_figures/FigureS4-12_multiCell_model/4_stepWise_regression.R",
                "supplementary_figures/FigureS4-12_multiCell_model/5_stepWise_regression_plot_IGF1.R",
                "supplementary_figures/FigureS4-12_multiCell_model/6_stepWise_regression_plot_HbA1c.R",
            ),
            upstream_products=("Multi-cell cNMF roots.", "Multi-cell program association outputs.", "GeneBayes posterior outputs."),
            expected_outputs=("multiCell_model/*.txt", "multiCell_model/*.pdf", "multiCell_model/conditional/*.pdf"),
            r_packages=("ComplexHeatmap", "circlize", "leaps", "pheatmap", "RColorBrewer", "reshape2"),
        ),
        PlotStageSpec(
            key="cnmf_validation_supplementary",
            name="Figure5 validation workflow",
            summary="Run the permutation, leave-one-out, cross-validation, and effect-direction legacy scripts.",
            scripts=(
                "Figure5/1_permutation_test_step1.R",
                "Figure5/2_permutation_test_step2.R",
                "Figure5/3_LeaveOneOut_test_step1.sh",
                "Figure5/4_LeaveOneOut_test_step2.R",
                "Figure5/5_split_crossValidation.R",
                "Figure5/6_predict_effectDirection.R",
                "Figure5/7_crossTraits_effectDirection.R",
            ),
            upstream_products=("cNMF spectra/regulation outputs.", "GeneBayes posterior outputs.", "Figure5 validation work dirs."),
            expected_outputs=("permutation_test/*", "loo_topgene_prediction/*", "FigS7_*.pdf", "*_genes_prediction.txt"),
            external_commands=("bash",),
            r_packages=("ggarchery", "ggplot2", "leaps", "plotrix", "plyr"),
        ),
    ),
)


PLOT_WORKFLOW_SPECS: dict[str, PlotWorkflowSpec] = {
    GENE_BURDEN_WORKFLOW.key: GENE_BURDEN_WORKFLOW,
    PERTURBSEQ_GENE_LEVEL_WORKFLOW.key: PERTURBSEQ_GENE_LEVEL_WORKFLOW,
    PERTURBSEQ_CNMF_WORKFLOW.key: PERTURBSEQ_CNMF_WORKFLOW,
    SUPPLEMENTARY_WORKFLOW.key: SUPPLEMENTARY_WORKFLOW,
}

PLOT_DEFAULT_INPUT_ARTIFACTS: dict[str, dict[str, tuple[str, ...]]] = {
    "gene_burden": {
        "genebayes_posterior_dir": ("genebayes", "posterior"),
    },
    "perturbseq_gene_level": {
        "gene_level_summary_dir": ("perturbseq", "gene_level", "K562GW"),
        "genebayes_posterior_dir": ("genebayes", "posterior"),
    },
    "perturbseq_cnmf": {
        "cnmf_root": ("perturbseq", "cnmf_genomewide", "cNMF"),
        "cnmf_regulation_dir": ("perturbseq", "cnmf_genomewide", "cNMF_regulation", "K562GW"),
        "program_association_dir": (
            "perturbseq",
            "cnmf_genomewide",
            "trait_association",
            "K562GW",
            "ProgramLevel",
        ),
        "genebayes_posterior_dir": ("genebayes", "posterior"),
    },
    "supplementary": {
        "gene_level_summary_dir": ("perturbseq", "gene_level", "K562GW"),
        "cnmf_root": ("perturbseq", "cnmf_genomewide", "cNMF"),
        "cnmf_regulation_dir": ("perturbseq", "cnmf_genomewide", "cNMF_regulation", "K562GW"),
        "program_association_dir": (
            "perturbseq",
            "cnmf_genomewide",
            "trait_association",
            "K562GW",
            "ProgramLevel",
        ),
        "genebayes_posterior_dir": ("genebayes", "posterior"),
    },
}

PLOT_DEFAULT_OUTPUT_ARTIFACTS: dict[str, tuple[str, ...]] = {
    "gene_burden": ("plots", "gene_burden"),
    "perturbseq_gene_level": ("plots", "perturbseq_gene_level"),
    "perturbseq_cnmf": ("plots", "perturbseq_cnmf"),
    "supplementary": ("plots", "supplementary"),
}

DEFAULT_POSTERIOR_TRAIT_FILES: dict[str, str] = {
    "MCH": "Backman_2021_86.per_gene_estimates.tsv",
    "RDW": "Backman_2021_88.per_gene_estimates.tsv",
    "IRF": "Backman_2021_106.per_gene_estimates.tsv",
}


def _resolve_required_path(config: LoadedConfig, value: object, label: str) -> Path:
    path = config.resolve_path(value if isinstance(value, (str, Path)) else None)
    if path is None:
        raise ConfigError(f"{label} must be a non-empty path")
    return path


def _resolve_script_path(config: LoadedConfig, relative_path: str) -> Path:
    _ = config
    return plot_script_path(relative_path)


def _prepare_output_alias_dir(root: Path, legacy_relative_path: str, generated_relative_path: str) -> None:
    link_output_path(root / "generated" / Path(generated_relative_path), root / Path(legacy_relative_path), is_directory=True)


def _prepare_output_alias_file(root: Path, legacy_relative_path: str, generated_relative_path: str) -> None:
    link_output_path(root / "generated" / Path(generated_relative_path), root / Path(legacy_relative_path), is_directory=False)


def _load_plot_workflow_context(config: LoadedConfig, spec: PlotWorkflowSpec) -> ResolvedPlotWorkflow:
    plot = config.workflow("plot")
    workflow = plot.get(spec.key, {})
    if not workflow:
        raise ConfigError(f"workflows.plot.{spec.key} is required")
    if not isinstance(workflow, dict):
        raise ConfigError(f"workflows.plot.{spec.key} must be a mapping")

    inputs = workflow.get("inputs", {})
    outputs = workflow.get("outputs", {})
    parameters = workflow.get("parameters", {})
    if not isinstance(inputs, dict):
        raise ConfigError(f"workflows.plot.{spec.key}.inputs must be a mapping")
    if not isinstance(outputs, dict):
        raise ConfigError(f"workflows.plot.{spec.key}.outputs must be a mapping")
    if not isinstance(parameters, dict):
        raise ConfigError(f"workflows.plot.{spec.key}.parameters must be a mapping")

    default_inputs = PLOT_DEFAULT_INPUT_ARTIFACTS.get(spec.key, {})
    resolved_inputs: dict[str, Path] = {}
    input_rows: list[tuple[str, str, Path]] = []
    for input_key, description in spec.input_descriptions:
        if input_key in inputs:
            path = _resolve_required_path(
                config,
                inputs[input_key],
                f"workflows.plot.{spec.key}.inputs.{input_key}",
            )
        elif input_key in default_inputs:
            path = config.artifact_path(*default_inputs[input_key])
        else:
            raise ConfigError(f"workflows.plot.{spec.key}.inputs.{input_key} is required")
        resolved_inputs[input_key] = path
        input_rows.append((input_key, description, path))

    if "output_dir" in outputs:
        output_dir = _resolve_required_path(
            config,
            outputs.get("output_dir"),
            f"workflows.plot.{spec.key}.outputs.output_dir",
        )
    else:
        output_dir = config.artifact_path(*PLOT_DEFAULT_OUTPUT_ARTIFACTS[spec.key])
    return ResolvedPlotWorkflow(
        config=config,
        spec=spec,
        inputs=resolved_inputs,
        input_rows=tuple(input_rows),
        output_dir=output_dir,
        parameters=dict(parameters),
    )


def _string_list(value: object, label: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (str, Path)):
        return (str(value),)
    if isinstance(value, list) and all(isinstance(item, (str, int, Path)) for item in value):
        return tuple(str(item) for item in value)
    raise ConfigError(f"{label} must be a string or list of strings")


def _int_list(value: object, label: str) -> tuple[int, ...]:
    if value is None:
        return ()
    if isinstance(value, int):
        return (value,)
    if isinstance(value, list) and all(isinstance(item, int) for item in value):
        return tuple(value)
    raise ConfigError(f"{label} must be an integer or list of integers")


def _mapping_list(value: object, label: str) -> tuple[dict[str, Any], ...]:
    if value is None:
        return ()
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise ConfigError(f"{label} must be a list of mappings")
    return tuple(dict(item) for item in value)


def _string_mapping(value: object, label: str) -> dict[str, str]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError(f"{label} must be a mapping")
    if not all(isinstance(key, str) and isinstance(item, (str, Path)) for key, item in value.items()):
        raise ConfigError(f"{label} must map strings to strings")
    return {str(key): str(item) for key, item in value.items()}


def _mapping_list_with_defaults(
    value: object,
    label: str,
    defaults: dict[str, Any],
) -> tuple[dict[str, Any], ...]:
    items = _mapping_list(value, label)
    if not items:
        return (dict(defaults),)
    return tuple({**defaults, **item} for item in items)


def _bool_parameter(value: object, label: str, *, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "yes", "on"}:
            return True
        if normalized in {"0", "false", "no", "off"}:
            return False
    raise ConfigError(f"{label} must be a boolean")


def _mode_parameter(
    value: object,
    label: str,
    *,
    allowed: tuple[str, ...],
    default: str,
) -> str:
    if value is None:
        return default
    if not isinstance(value, str):
        raise ConfigError(f"{label} must be one of: {', '.join(allowed)}")
    normalized = value.strip().lower()
    if normalized not in allowed:
        raise ConfigError(f"{label} must be one of: {', '.join(allowed)}")
    return normalized


def _plot_rscript_prefix(ctx: ResolvedPlotWorkflow) -> list[str]:
    return ctx.config.executable_command(
        "plot_rscript",
        ctx.config.executable_command("rscript", ["Rscript"]),
    )


def _plot_bash_prefix(ctx: ResolvedPlotWorkflow) -> list[str]:
    return ctx.config.executable_command("plot_bash", ["bash"])


def _ldsc_py(ctx: ResolvedPlotWorkflow) -> str:
    return ctx.config.executable("ldsc_py", "ldsc.py")


def _lift_over(ctx: ResolvedPlotWorkflow) -> str:
    return ctx.config.executable("liftover", "liftOver")


def _ldsc_baseline_prefix(ctx: ResolvedPlotWorkflow) -> Path:
    return ctx.input_path("ldsc_baseline_dir") / "baseline."


def _rscript_command(ctx: ResolvedPlotWorkflow, relative_script: str, *args: object) -> list[str]:
    return [
        *_plot_rscript_prefix(ctx),
        str(_resolve_script_path(ctx.config, relative_script)),
        *(str(arg) for arg in args),
    ]


def _bash_command(ctx: ResolvedPlotWorkflow, relative_script: str, *args: object) -> list[str]:
    return [
        *_plot_bash_prefix(ctx),
        str(_resolve_script_path(ctx.config, relative_script)),
        *(str(arg) for arg in args),
    ]


def _link_gene_maps(root: Path, gene_map: Path) -> None:
    link_alias(gene_map, root / "data" / "gencode_v41_gname_gid_ALL_sorted_onlyID")
    link_alias(gene_map, root / "data" / "gencode_v41_gname_gid_ALL_sorted")


def _resolve_optional_parameter_path(
    ctx: ResolvedPlotWorkflow,
    value: object,
    label: str,
    *,
    base_dir: Path | None = None,
    default: Path,
) -> Path:
    if value is None:
        return default.resolve()
    if not isinstance(value, (str, Path)):
        raise ConfigError(f"{label} must be a path string")
    path = Path(value)
    if path.is_absolute():
        return path
    if base_dir is not None and len(path.parts) == 1:
        return (base_dir / path).resolve()
    resolved = ctx.config.resolve_path(path)
    if resolved is None:
        raise ConfigError(f"{label} must be a non-empty path")
    return resolved


def _posterior_trait_file_paths(
    ctx: ResolvedPlotWorkflow,
    parameter_name: str,
    label: str,
) -> dict[str, Path]:
    configured = _string_mapping(ctx.parameter(parameter_name), label)
    posterior_dir = ctx.input_path("genebayes_posterior_dir")
    return {
        trait: _resolve_optional_parameter_path(
            ctx,
            configured.get(trait),
            f"{label}.{trait}",
            base_dir=posterior_dir,
            default=posterior_dir / filename,
        )
        for trait, filename in DEFAULT_POSTERIOR_TRAIT_FILES.items()
    }


def _prepare_legacy_cnmf_root(source_root: Path, target_root: Path) -> None:
    if target_root.exists() or target_root.is_symlink():
        replace_tree(target_root)
    ensure_directories([target_root])
    populate_directory(source_root, target_root)

    aggregate_dir = source_root / "cNMF_all"
    if aggregate_dir.is_dir():
        for child in aggregate_dir.iterdir():
            if child.is_file() and child.name.startswith("cNMF_all."):
                link_alias(child, target_root / child.name)


def _create_static_stage(spec: PlotStageSpec, cwd: Path, commands: tuple[list[str], ...]) -> PlannedStage:
    preview_lines = tuple(format_command_preview(command, cwd=cwd) for command in commands)

    def _action() -> None:
        run_commands(commands, cwd=cwd)

    return PlannedStage(spec=spec, preview_lines=preview_lines, action=_action)


def _parallel_stage_workers(configured_jobs: int | None, *, max_jobs: int) -> int:
    if max_jobs < 1:
        return 1
    if configured_jobs is not None:
        return max(1, min(configured_jobs, max_jobs))

    for env_key in ("SLURM_CPUS_PER_TASK", "PAPER_PIPELINE_PLOT_PARALLEL_JOBS"):
        raw = os.environ.get(env_key)
        if raw is None:
            continue
        try:
            jobs = int(raw)
        except ValueError:
            continue
        if jobs > 0:
            return min(jobs, max_jobs)
    return 1


def _plot_parallel_stage_env() -> dict[str, str]:
    return {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "R_DATATABLE_NUM_THREADS": "1",
    }


def _create_parallel_stage(
    spec: PlotStageSpec,
    cwd: Path,
    *,
    prepare_commands: tuple[list[str], ...] = (),
    parallel_commands: tuple[list[str], ...] = (),
    finalize_commands: tuple[list[str], ...] = (),
    configured_parallel_jobs: int | None = None,
) -> PlannedStage:
    preview_lines: list[str] = []
    preview_lines.extend(format_command_preview(command, cwd=cwd) for command in prepare_commands)
    if parallel_commands:
        parallel_workers = configured_parallel_jobs if configured_parallel_jobs is not None else "auto"
        preview_lines.append(
            f"Run {len(parallel_commands)} program jobs in parallel (max workers: {parallel_workers}; "
            "runtime auto falls back to SLURM_CPUS_PER_TASK or 1)."
        )
        preview_lines.extend(format_command_preview(command, cwd=cwd) for command in parallel_commands)
    preview_lines.extend(format_command_preview(command, cwd=cwd) for command in finalize_commands)

    def _action() -> None:
        run_commands(prepare_commands, cwd=cwd)
        if parallel_commands:
            max_workers = _parallel_stage_workers(configured_parallel_jobs, max_jobs=len(parallel_commands))
            parallel_env = _plot_parallel_stage_env()
            if max_workers <= 1 or len(parallel_commands) == 1:
                run_commands(parallel_commands, cwd=cwd, env=parallel_env)
            else:
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    futures = [
                        executor.submit(run_command, command, cwd=cwd, env=parallel_env)
                        for command in parallel_commands
                    ]
                    for future in futures:
                        future.result()
        run_commands(finalize_commands, cwd=cwd)

    return PlannedStage(spec=spec, preview_lines=tuple(preview_lines), action=_action)


def _create_skipped_stage(spec: PlotStageSpec, message: str) -> PlannedStage:
    def _action() -> None:
        return None

    return PlannedStage(spec=spec, preview_lines=(f"Skip stage: {message}",), action=_action)


def _build_stage_inventory(
    config: LoadedConfig,
    spec: PlotWorkflowSpec,
    output_dir: Path,
    stage_plans: tuple[PlannedStage, ...],
) -> None:
    inventory_path = output_dir / "stage_inventory.tsv"
    stage_plan_map = {stage.spec.key: stage for stage in stage_plans}
    with inventory_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "stage_key",
                "stage_name",
                "scripts",
                "script_exists",
                "summary",
                "upstream_products",
                "expected_outputs",
                "r_packages",
                "external_commands",
                "notes",
                "planned_execution",
            ]
        )
        for stage in spec.stages:
            resolved_scripts = tuple(
                (relative_path, _resolve_script_path(config, relative_path)) for relative_path in stage.scripts
            )
            plan = stage_plan_map.get(stage.key)
            writer.writerow(
                [
                    stage.key,
                    stage.name,
                    "; ".join(relative_path for relative_path, _ in resolved_scripts),
                    "; ".join(
                        f"{relative_path}={'yes' if resolved_path.exists() else 'no'}"
                        for relative_path, resolved_path in resolved_scripts
                    ),
                    stage.summary,
                    "; ".join(stage.upstream_products),
                    "; ".join(stage.expected_outputs),
                    "; ".join(stage.r_packages),
                    "; ".join(stage.external_commands),
                    "; ".join(stage.notes),
                    " || ".join(plan.preview_lines) if plan else "",
                ]
            )


def _write_plot_manifest(ctx: ResolvedPlotWorkflow, stage_plans: tuple[PlannedStage, ...]) -> None:
    manifest_path = ctx.output_dir / "workflow_manifest.md"
    stage_plan_map = {stage.spec.key: stage for stage in stage_plans}
    lines = [
        f"# {ctx.spec.title}",
        "",
        ctx.spec.summary,
        "",
        "This workflow prepares a legacy-compatible workspace under the configured output directory",
        "and executes package-owned plotting scripts vendored under `src/paper_pipeline/plot/scripts/`.",
        "",
        "## Upstream workflows",
    ]
    lines.extend(f"- `{name}`" for name in ctx.spec.upstream_workflows)
    lines.extend(
        [
            "",
            "## Configured inputs",
            "| key | description | resolved path | exists |",
            "| --- | --- | --- | --- |",
        ]
    )
    for key, description, path in ctx.input_rows:
        lines.append(f"| `{key}` | {description} | `{path}` | {'yes' if path.exists() else 'no'} |")

    lines.extend(["", "## Output directory", f"- `{ctx.output_dir}`", "", "## Stage plan"])
    for stage in ctx.spec.stages:
        plan = stage_plan_map.get(stage.key)
        lines.extend(["", f"### {stage.name}", "", stage.summary, "", "Scripts:"])
        for relative_path in stage.scripts:
            script_path = _resolve_script_path(ctx.config, relative_path)
            lines.append(f"- `{relative_path}` (exists: {'yes' if script_path.exists() else 'no'})")
        lines.extend(["", "Expected outputs:"])
        lines.extend(f"- `{item}`" for item in stage.expected_outputs)
        if plan is not None:
            lines.extend(["", "Planned execution:"])
            lines.extend(f"- `{line}`" for line in plan.preview_lines)
    write_text(manifest_path, "\n".join(lines) + "\n")
    _build_stage_inventory(ctx.config, ctx.spec, ctx.output_dir, stage_plans)


def _build_plot_workflow_tasks(
    config: LoadedConfig,
    spec: PlotWorkflowSpec,
    prepare_workspace: Callable[[ResolvedPlotWorkflow], None],
    build_stage_plans: Callable[[ResolvedPlotWorkflow], tuple[PlannedStage, ...]],
) -> list[Task]:
    ctx = _load_plot_workflow_context(config, spec)
    stage_plans = build_stage_plans(ctx)
    manifest_path = ctx.output_dir / "workflow_manifest.md"
    inventory_path = ctx.output_dir / "stage_inventory.tsv"

    return [
        Task(
            name=f"Create {spec.cli_name} output directory",
            preview=f"Create output directory: {ctx.output_dir}",
            action=lambda: ensure_directories([ctx.output_dir]),
        ),
        Task(
            name=f"Prepare {spec.cli_name} legacy workspace",
            preview=f"Link plotting inputs into legacy workspace rooted at {ctx.output_dir}",
            action=lambda: prepare_workspace(ctx),
        ),
        Task(
            name=f"Write {spec.cli_name} manifest and inventory",
            preview=f"Write {manifest_path} and {inventory_path}",
            action=lambda: _write_plot_manifest(ctx, stage_plans),
        ),
        *[
            Task(
                name=f"{spec.cli_name}: {planned_stage.spec.name}",
                preview=" ; ".join(planned_stage.preview_lines),
                action=planned_stage.action,
            )
            for planned_stage in stage_plans
        ],
    ]


def _prepare_gene_burden_workspace(ctx: ResolvedPlotWorkflow) -> None:
    root = ctx.output_dir
    ensure_directories(
        [
            root / "data",
            root / "data" / "LoF",
            root / "data" / "GWAS",
        ]
    )
    _prepare_output_alias_dir(root, "data/GWAS/closest_genes", "gene_burden/closest_genes")
    _prepare_output_alias_dir(root, "clump", "gene_burden/clump")
    _prepare_output_alias_dir(root, "enrichment_result", "gene_burden/enrichment_result")
    link_path(ctx.input_path("genebayes_posterior_dir"), root / "data" / "LoF" / "GeneBayes_posterior")
    link_path(ctx.input_path("raw_burden_dir"), root / "data" / "LoF" / "raw_burden")
    populate_directory(
        ctx.input_path("gwas_dir"),
        root / "data" / "GWAS",
        exclude_names=("closest_genes",),
    )
    hla_genes = ctx.input_path("gwas_dir") / "closest_genes" / "HLA.genes"
    if hla_genes.exists():
        link_alias(hla_genes, root / "data" / "GWAS" / "closest_genes" / "HLA.genes")
    link_path(ctx.input_path("geneset_dir"), root / "data" / "geneset")
    _link_gene_maps(root, ctx.input_path("gene_map"))
    link_alias(ctx.input_path("trait_correspondence"), root / "data" / "Backman_Niele_corresp")


def _prepare_gene_level_workspace(ctx: ResolvedPlotWorkflow) -> None:
    root = ctx.output_dir
    ensure_directories(
        [
            root / "data",
            root / "data" / "LoF",
            root / "data" / "Perturbseq" / "trait_association" / "K562GW",
            root / "data" / "BurdenRegCor",
        ]
    )
    _prepare_output_alias_dir(
        root,
        "data/Perturbseq/trait_association/K562GW/GeneLevel",
        "perturbseq_gene_level/gene_level_trait_association",
    )
    link_path(ctx.input_path("genebayes_posterior_dir"), root / "data" / "LoF" / "GeneBayes_posterior")
    link_path(ctx.input_path("raw_burden_dir"), root / "data" / "LoF" / "raw_burden")
    link_path(ctx.input_path("gene_level_summary_dir"), root / "data" / "Perturbseq" / "geneLevel" / "K562GW")
    link_path(ctx.input_path("gwas_dir"), root / "data" / "GWAS")
    link_path(ctx.input_path("burden_reg_correlation_dir"), root / "data" / "BurdenRegCor" / "GeneLevel")
    link_alias(
        ctx.input_path("gene_level_trait_metadata"),
        root / "data" / "Perturbseq" / "trait_association" / "K562GW" / "GeneLevel" / "Backman_Niele_corresp_all",
    )
    _link_gene_maps(root, ctx.input_path("gene_map"))
    link_alias(ctx.input_path("shet_path"), root / "data" / "shet_10bins.txt")


def _prepare_cnmf_workspace(ctx: ResolvedPlotWorkflow) -> None:
    root = ctx.output_dir
    ensure_directories(
        [
            root / "data",
            root / "data" / "LoF",
            root / "data" / "Perturbseq" / "cNMF",
            root / "data" / "Perturbseq" / "cNMF_regulation",
            root / "data" / "Perturbseq" / "trait_association" / "K562GW",
        ]
    )
    _prepare_output_alias_dir(
        root,
        "data/Perturbseq/trait_association/K562GW/ProgramLevel",
        "perturbseq_cnmf/program_level",
    )
    _prepare_output_alias_dir(root, "program_trans", "perturbseq_cnmf/program_trans")
    link_path(ctx.input_path("genebayes_posterior_dir"), root / "data" / "LoF" / "GeneBayes_posterior")
    _prepare_legacy_cnmf_root(
        ctx.input_path("cnmf_root"),
        root / "data" / "Perturbseq" / "cNMF" / "K562GW",
    )
    link_path(ctx.input_path("cnmf_regulation_dir"), root / "data" / "Perturbseq" / "cNMF_regulation" / "K562GW")
    link_path(ctx.input_path("gwas_dir"), root / "data" / "GWAS")
    link_alias(ctx.input_path("shet_path"), root / "data" / "shet_10bins.txt")
    link_alias(ctx.input_path("trans_eqtl_path"), root / "data" / ctx.input_path("trans_eqtl_path").name)
    _link_gene_maps(root, ctx.input_path("gene_map"))


def _prepare_supplementary_workspace(ctx: ResolvedPlotWorkflow) -> None:
    root = ctx.output_dir
    ensure_directories([root / "data"])
    run_chip_collection = bool(ctx.parameter("run_chip_collection", False))

    populate_directory(
        ctx.input_path("ldsc_dir"),
        root / "data" / "LDSC",
        exclude_names=("ANNOTATION_EUR", "Partitioning", "result_summary", "figure", "genetic_correlation"),
    )
    _prepare_output_alias_dir(root, "data/LDSC/ANNOTATION_EUR", "supplementary/ldsc/ANNOTATION_EUR")
    _prepare_output_alias_dir(root, "data/LDSC/Partitioning", "supplementary/ldsc/Partitioning")
    _prepare_output_alias_dir(root, "data/LDSC/result_summary", "supplementary/ldsc/result_summary")
    _prepare_output_alias_dir(root, "data/LDSC/figure", "supplementary/ldsc/figure")
    _prepare_output_alias_dir(root, "data/LDSC/genetic_correlation", "supplementary/ldsc/genetic_correlation")
    chip_excludes = [
        "chip_score",
        "cNMF_enrichment",
        "chip_score_KOeffects_MWUtest.txt",
        "selected_scores_for_prediction.txt",
    ]
    if run_chip_collection:
        chip_excludes.append("hg19")
    populate_directory(
        ctx.input_path("chip_dir"),
        root / "data" / "TF_ChIP",
        exclude_names=tuple(chip_excludes),
    )
    ensure_directories(
        [
            root / "data" / "TF_ChIP",
            root / "data" / "Perturbseq" / "cNMF",
            root / "data" / "Perturbseq" / "cNMF_regulation",
            root / "data" / "Perturbseq" / "trait_association",
            root / "data" / "Perturbseq" / "metadata",
            root / "data" / "LoF",
        ]
    )
    if run_chip_collection:
        _prepare_output_alias_dir(root, "data/TF_ChIP/hg19", "supplementary/chip/hg19")
    _prepare_output_alias_dir(root, "data/TF_ChIP/chip_score", "supplementary/chip/chip_score")
    _prepare_output_alias_dir(root, "data/TF_ChIP/cNMF_enrichment", "supplementary/chip/cNMF_enrichment")
    _prepare_output_alias_file(
        root,
        "data/TF_ChIP/chip_score_KOeffects_MWUtest.txt",
        "supplementary/chip/chip_score_KOeffects_MWUtest.txt",
    )
    _prepare_output_alias_file(
        root,
        "data/TF_ChIP/selected_scores_for_prediction.txt",
        "supplementary/chip/selected_scores_for_prediction.txt",
    )
    _prepare_output_alias_dir(root, "data/fitness_compare", "supplementary/edfig5/fitness_compare")
    _prepare_output_alias_file(
        root,
        "data/programs_causal_relation_summary.txt",
        "supplementary/edfig5/programs_causal_relation_summary.txt",
    )
    _prepare_output_alias_dir(root, "multiCell_model", "supplementary/multi_cell/multiCell_model")
    _prepare_output_alias_dir(root, "permutation_test", "supplementary/validation/permutation_test")
    _prepare_output_alias_dir(root, "loo_topgene_prediction", "supplementary/validation/loo_topgene_prediction")

    populate_directory(
        ctx.input_path("multicell_cnmf_root"),
        root / "data" / "Perturbseq" / "cNMF",
        exclude_names=("K562GW",),
    )
    populate_directory(
        ctx.input_path("multicell_cnmf_regulation_dir"),
        root / "data" / "Perturbseq" / "cNMF_regulation",
        exclude_names=("K562GW",),
    )
    populate_directory(
        ctx.input_path("multicell_trait_association_dir"),
        root / "data" / "Perturbseq" / "trait_association",
        exclude_names=("K562GW",),
    )

    _prepare_legacy_cnmf_root(
        ctx.input_path("cnmf_root"),
        root / "data" / "Perturbseq" / "cNMF" / "K562GW",
    )
    link_path(ctx.input_path("cnmf_regulation_dir"), root / "data" / "Perturbseq" / "cNMF_regulation" / "K562GW")
    ensure_directories([root / "data" / "Perturbseq" / "trait_association" / "K562GW"])
    link_path(
        ctx.input_path("program_association_dir"),
        root / "data" / "Perturbseq" / "trait_association" / "K562GW" / "ProgramLevel",
    )
    link_path(ctx.input_path("gene_level_summary_dir"), root / "data" / "Perturbseq" / "geneLevel" / "K562GW")
    link_path(ctx.input_path("genebayes_posterior_dir"), root / "data" / "LoF" / "GeneBayes_posterior")
    link_path(ctx.input_path("gwas_dir"), root / "data" / "GWAS")
    link_path(ctx.input_path("geneset_dir"), root / "data" / "geneset")
    link_alias(ctx.input_path("growth_screening_csv"), root / "data" / "Growth_screening_2014.csv")
    link_alias(ctx.input_path("metadata"), root / "data" / "Perturbseq" / "metadata" / "gwps_metadata.csv")
    link_alias(ctx.input_path("shet_path"), root / "data" / "shet_10bins.txt")
    _link_gene_maps(root, ctx.input_path("gene_map"))


def _build_gene_burden_stage_plans(ctx: ResolvedPlotWorkflow) -> tuple[PlannedStage, ...]:
    gwas_traits = _string_list(ctx.parameter("gwas_traits"), "workflows.plot.gene_burden.parameters.gwas_traits")
    if not gwas_traits:
        gwas_traits = ("30050",)
    manhattan_trait = str(ctx.parameter("manhattan_gwas_trait", gwas_traits[0]))
    volcano_burden_file = str(ctx.parameter("volcano_burden_file", "Backman_2021_86_M1_001.summary_statistics.csv"))
    cross_trait_traits = _string_list(
        ctx.parameter("cross_trait_traits"),
        "workflows.plot.gene_burden.parameters.cross_trait_traits",
    ) or ("86", "88")
    enrichment_pairs = _mapping_list(
        ctx.parameter("lof_enrichment_pairs"),
        "workflows.plot.gene_burden.parameters.lof_enrichment_pairs",
    ) or (
        {"burden_trait": "86", "gwas_trait": "30050", "trait_name": "MCH"},
    )

    clumping_commands = tuple(
        command
        for trait in gwas_traits
        for command in (
            _bash_command(ctx, "Figure2/3_GWAS_clumping.sh", trait, ctx.input_path("ld_reference_dir")),
            _bash_command(ctx, "Figure2/4_annotate_closest_genes.sh", trait),
        )
    )

    figure_commands: list[list[str]] = [
        _rscript_command(ctx, "Figure2/1_manhattanplot.R", manhattan_trait),
        _rscript_command(ctx, "Figure2/2_volcanoplot.R", volcano_burden_file),
    ]
    for pair in enrichment_pairs:
        figure_commands.append(
            _rscript_command(ctx, "Figure2/5_pathway_enrichment_GWAS.R", pair["gwas_trait"])
        )
        figure_commands.append(
            _rscript_command(
                ctx,
                "Figure2/6_pathway_enrichment_LoF.R",
                pair["burden_trait"],
                pair["gwas_trait"],
            )
        )
        figure_commands.append(
            _rscript_command(ctx, "Figure2/7_pathway_enrichment_plot.R", pair["trait_name"])
        )
    figure_commands.append(
        _rscript_command(ctx, "Figure2/8_CrossTrait_plot.R", ",".join(cross_trait_traits))
    )

    return (
        _create_static_stage(ctx.spec.stages[0], ctx.output_dir, clumping_commands),
        _create_static_stage(ctx.spec.stages[1], ctx.output_dir, tuple(figure_commands)),
    )


def _build_gene_level_stage_plans(ctx: ResolvedPlotWorkflow) -> tuple[PlannedStage, ...]:
    trait_files = _string_list(
        ctx.parameter("trait_files"),
        "workflows.plot.perturbseq_gene_level.parameters.trait_files",
    ) or ("Backman_2021_86.per_gene_estimates.tsv",)
    enrichment_gwas_trait = str(ctx.parameter("enrichment_gwas_trait", "30050"))
    enrichment_comp = str(ctx.parameter("enrichment_comp", "closest"))
    enrichment_target_gene = str(ctx.parameter("enrichment_target_gene", "HBA1"))
    enrichment_lof_burden_file = str(
        ctx.parameter("enrichment_lof_burden_file", "Backman_2021_86_M1_001.summary_statistics.csv")
    )
    render_trait_file = str(ctx.parameter("render_trait_file", trait_files[0]))

    association_commands = tuple(
        _rscript_command(ctx, "Figure3/1_gamma_beta_cor.R", trait_file)
        for trait_file in trait_files
    )
    enrichment_commands = (
        _rscript_command(
            ctx,
            "Figure3/2_GWAS_enrich.R",
            enrichment_gwas_trait,
            enrichment_comp,
            enrichment_target_gene,
            enrichment_lof_burden_file,
        ),
    )
    rendering_commands = (
        _rscript_command(ctx, "Figure3/3_plot_gamma_vs_cor.R", render_trait_file),
        _rscript_command(ctx, "Figure3/4_plot_QQ.R"),
    )

    return (
        _create_static_stage(ctx.spec.stages[0], ctx.output_dir, association_commands),
        _create_static_stage(ctx.spec.stages[1], ctx.output_dir, enrichment_commands),
        _create_static_stage(ctx.spec.stages[2], ctx.output_dir, rendering_commands),
    )


def _build_cnmf_stage_plans(ctx: ResolvedPlotWorkflow) -> tuple[PlannedStage, ...]:
    k = int(ctx.parameter("k", 60))
    trait_files = _string_list(
        ctx.parameter("trait_files"),
        "workflows.plot.perturbseq_cnmf.parameters.trait_files",
    ) or ("Backman_2021_86.per_gene_estimates.tsv",)
    mode = _mode_parameter(
        ctx.parameter("mode"),
        "workflows.plot.perturbseq_cnmf.parameters.mode",
        allowed=("generic", "legacy"),
        default="generic",
    )
    run_multiple_regression = _bool_parameter(
        ctx.parameter("run_multiple_regression"),
        "workflows.plot.perturbseq_cnmf.parameters.run_multiple_regression",
        default=(mode == "legacy"),
    )
    run_trans_eqtl_follow_up = _bool_parameter(
        ctx.parameter("run_trans_eqtl_follow_up"),
        "workflows.plot.perturbseq_cnmf.parameters.run_trans_eqtl_follow_up",
        default=(mode == "legacy"),
    )
    plot_label_programs = _int_list(
        ctx.parameter("plot_label_programs"),
        "workflows.plot.perturbseq_cnmf.parameters.plot_label_programs",
    ) or (4, 16, 25, 40)
    corregulation_pairs = _mapping_list(
        ctx.parameter("corregulation_pairs"),
        "workflows.plot.perturbseq_cnmf.parameters.corregulation_pairs",
    ) or ({"program_a": "P25", "program_b": "P16"},)
    trans_eqtl_top_n = int(ctx.parameter("trans_eqtl_top_n", 100))
    trans_eqtl_programs = _int_list(
        ctx.parameter("trans_eqtl_programs"),
        "workflows.plot.perturbseq_cnmf.parameters.trans_eqtl_programs",
    ) or tuple(range(1, k + 1))
    posterior_trait_files = (
        _posterior_trait_file_paths(
            ctx,
            "posterior_trait_files",
            "workflows.plot.perturbseq_cnmf.parameters.posterior_trait_files",
        )
        if run_multiple_regression
        else {}
    )

    trans_eqtl_gwas_trait = "30050"
    trans_eqtl_gwas_file: Path | None = None
    trans_eqtl_output_label = "MCH"
    trans_eqtl_parallel_jobs: int | None = None
    trans_eqtl_regulator_file: Path | None = None
    if run_trans_eqtl_follow_up:
        trans_eqtl_gwas_trait = str(ctx.parameter("trans_eqtl_gwas_trait", "30050"))
        trans_eqtl_gwas_file = _resolve_optional_parameter_path(
            ctx,
            ctx.parameter("trans_eqtl_gwas_file"),
            "workflows.plot.perturbseq_cnmf.parameters.trans_eqtl_gwas_file",
            base_dir=ctx.input_path("gwas_dir"),
            default=ctx.input_path("gwas_dir") / f"{trans_eqtl_gwas_trait}_irnt.tsv.gz",
        )
        trans_eqtl_output_label = str(ctx.parameter("trans_eqtl_output_label", "MCH"))
        trans_eqtl_parallel_jobs_param = ctx.parameter("trans_eqtl_parallel_jobs")
        trans_eqtl_parallel_jobs = None if trans_eqtl_parallel_jobs_param is None else int(trans_eqtl_parallel_jobs_param)
        if trans_eqtl_parallel_jobs is not None and trans_eqtl_parallel_jobs < 1:
            raise ConfigError("workflows.plot.perturbseq_cnmf.parameters.trans_eqtl_parallel_jobs must be >= 1")
        trans_eqtl_regulator_file = _resolve_optional_parameter_path(
            ctx,
            ctx.parameter("trans_eqtl_regulator_file"),
            "workflows.plot.perturbseq_cnmf.parameters.trans_eqtl_regulator_file",
            base_dir=ctx.input_path("program_association_dir"),
            default=ctx.input_path("program_association_dir") / f"regulators_enrichment_K{k}_{trait_files[0]}",
        )

    association_commands = tuple(
        _rscript_command(ctx, "Figure4/1_burden_program_regulators_test.R", trait_file, k)
        for trait_file in trait_files
    )

    rendering_commands: list[list[str]] = []
    label_arg = ",".join(str(program) for program in plot_label_programs)
    for trait_file in trait_files:
        rendering_commands.append(
            _rscript_command(ctx, "Figure4/2_plot_burden_program_regulators.R", trait_file, k, label_arg)
        )
    for pair in corregulation_pairs:
        rendering_commands.append(
            _rscript_command(
                ctx,
                "Figure4/3_corregulation_plot.R",
                pair["program_a"],
                pair["program_b"],
                k,
            )
        )
    if run_multiple_regression:
        rendering_commands.append(
            _rscript_command(
                ctx,
                "Figure4/4_multipleRegression.R",
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
            )
        )

    trans_eqtl_prepare_commands: tuple[list[str], ...] = ()
    trans_eqtl_program_commands: tuple[list[str], ...] = ()
    trans_eqtl_finalize_commands: tuple[list[str], ...] = ()
    if run_trans_eqtl_follow_up:
        assert trans_eqtl_gwas_file is not None
        assert trans_eqtl_regulator_file is not None
        trans_eqtl_prepare_commands = (
            _rscript_command(ctx, "Figure4/5_1_transeQTL_prepare.R", k, trans_eqtl_top_n),
        )
        trans_eqtl_program_commands = tuple(
            _rscript_command(
                ctx,
                "Figure4/5_2_transeQTL_GWASalleles.R",
                program,
                trans_eqtl_top_n,
                trans_eqtl_gwas_file,
                trans_eqtl_output_label,
            )
            for program in trans_eqtl_programs
        )
        trans_eqtl_finalize_commands = (
            _rscript_command(
                ctx,
                "Figure4/5_3_plot_transeQTL_vs_Perturbmodel.R",
                k,
                "program_trans/result",
                trans_eqtl_regulator_file,
                trans_eqtl_output_label,
            ),
        )

    trans_eqtl_stage = (
        _create_parallel_stage(
            ctx.spec.stages[2],
            ctx.output_dir,
            prepare_commands=tuple(trans_eqtl_prepare_commands),
            parallel_commands=trans_eqtl_program_commands,
            finalize_commands=trans_eqtl_finalize_commands,
            configured_parallel_jobs=trans_eqtl_parallel_jobs,
        )
        if run_trans_eqtl_follow_up
        else _create_skipped_stage(
            ctx.spec.stages[2],
            "workflows.plot.perturbseq_cnmf.parameters.run_trans_eqtl_follow_up is false; "
            "skip Figure4 trans-eQTL follow-up.",
        )
    )

    return (
        _create_static_stage(ctx.spec.stages[0], ctx.output_dir, association_commands),
        _create_static_stage(ctx.spec.stages[1], ctx.output_dir, tuple(rendering_commands)),
        trans_eqtl_stage,
    )


def _build_ldsc_stage(ctx: ResolvedPlotWorkflow) -> PlannedStage:
    stage = ctx.spec.stages[0]
    ldsc_cwd = ctx.output_dir / "data" / "LDSC"
    annotation_beds = _string_list(
        ctx.parameter("ldsc_annotation_beds"),
        "workflows.plot.supplementary.parameters.ldsc_annotation_beds",
    )
    ldsc_sumstats = _mapping_list(
        ctx.parameter("ldsc_sumstats"),
        "workflows.plot.supplementary.parameters.ldsc_sumstats",
    )
    genetic_pairs = _mapping_list(
        ctx.parameter("genetic_correlation_pairs"),
        "workflows.plot.supplementary.parameters.genetic_correlation_pairs",
    )
    if not ldsc_sumstats:
        return _create_skipped_stage(
            stage,
            "workflows.plot.supplementary.parameters.ldsc_sumstats is empty; "
            "Figure1 LDSC plotting is disabled.",
        )
    if not annotation_beds:
        annotation_beds = tuple(path.stem for path in sorted((ldsc_cwd / "bed").glob("*.bed")))

    preview_lines = []
    preview_lines.extend(
        format_command_preview(
            _bash_command(
                ctx,
                "Figure1/1_annotation_LDSC.sh",
                bed_name,
                ctx.input_path("ldsc_ld_reference_dir"),
                ctx.input_path("ldsc_baseline_dir"),
                _ldsc_py(ctx),
            ),
            cwd=ldsc_cwd,
        )
        for bed_name in annotation_beds
    )
    preview_lines.extend(
        format_command_preview(
            _bash_command(
                ctx,
                "Figure1/2_stratify_LDSC.sh",
                ctx.input_path("ldsc_sumstats_dir") / item["file"],
                item["label"],
                ctx.input_path("ldsc_weights_dir"),
                ctx.input_path("ldsc_frq_dir"),
                _ldsc_baseline_prefix(ctx),
                _ldsc_py(ctx),
            ),
            cwd=ldsc_cwd,
        )
        for item in ldsc_sumstats
    )
    preview_lines.append(
        format_command_preview(
            _bash_command(ctx, "Figure1/3_summarize_stratify_results.sh"),
            cwd=ldsc_cwd,
        )
    )
    preview_lines.extend(
        format_command_preview(
            _bash_command(
                ctx,
                "Figure1/4_geneticCorrelation.sh",
                item["label_a"],
                item["label_b"],
                ctx.input_path("ldsc_sumstats_dir") / item["file_a"],
                ctx.input_path("ldsc_sumstats_dir") / item["file_b"],
                ctx.input_path("ldsc_ld_reference_dir"),
                _ldsc_py(ctx),
            ),
            cwd=ldsc_cwd,
        )
        for item in genetic_pairs
    )
    preview_lines.append(
        format_command_preview(
            _rscript_command(ctx, "Figure1/5_makePlot.R", "data/LDSC"),
            cwd=ctx.output_dir,
        )
    )

    def _action() -> None:
        for bed_name in annotation_beds:
            run_command(
                _bash_command(
                    ctx,
                    "Figure1/1_annotation_LDSC.sh",
                    bed_name,
                    ctx.input_path("ldsc_ld_reference_dir"),
                    ctx.input_path("ldsc_baseline_dir"),
                    _ldsc_py(ctx),
                ),
                cwd=ldsc_cwd,
            )
        for item in ldsc_sumstats:
            run_command(
                _bash_command(
                    ctx,
                    "Figure1/2_stratify_LDSC.sh",
                    ctx.input_path("ldsc_sumstats_dir") / item["file"],
                    item["label"],
                    ctx.input_path("ldsc_weights_dir"),
                    ctx.input_path("ldsc_frq_dir"),
                    _ldsc_baseline_prefix(ctx),
                    _ldsc_py(ctx),
                ),
                cwd=ldsc_cwd,
            )
        run_command(_bash_command(ctx, "Figure1/3_summarize_stratify_results.sh"), cwd=ldsc_cwd)
        for item in genetic_pairs:
            run_command(
                _bash_command(
                    ctx,
                    "Figure1/4_geneticCorrelation.sh",
                    item["label_a"],
                    item["label_b"],
                    ctx.input_path("ldsc_sumstats_dir") / item["file_a"],
                    ctx.input_path("ldsc_sumstats_dir") / item["file_b"],
                    ctx.input_path("ldsc_ld_reference_dir"),
                    _ldsc_py(ctx),
                ),
                cwd=ldsc_cwd,
            )
        run_command(_rscript_command(ctx, "Figure1/5_makePlot.R", "data/LDSC"), cwd=ctx.output_dir)

    return PlannedStage(spec=stage, preview_lines=tuple(preview_lines), action=_action)


def _build_chip_stage(ctx: ResolvedPlotWorkflow) -> PlannedStage:
    stage = ctx.spec.stages[1]
    k = int(ctx.parameter("k", 60))
    chip_programs = _int_list(
        ctx.parameter("chip_programs"),
        "workflows.plot.supplementary.parameters.chip_programs",
    ) or tuple(range(1, k + 1))
    run_collection = bool(ctx.parameter("run_chip_collection", False))
    preview_lines: list[str] = []
    if run_collection:
        preview_lines.append(
            format_command_preview(
                _bash_command(
                    ctx,
                    "supplementary_figures/EDFig4_ChIP/1_collect_ChIP_K562.sh",
                    _lift_over(ctx),
                    ctx.input_path("liftOver_chain"),
                ),
                cwd=ctx.output_dir,
            )
        )
    preview_lines.extend(
        [
            format_command_preview(
                _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/2_calculate_TF_chipScore.R"),
                cwd=ctx.output_dir,
            ),
            format_command_preview(
                _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/3_ChIPscore_program_enrichment.R"),
                cwd=ctx.output_dir,
            ),
            format_command_preview(
                _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/4_select_useful_experiments.R"),
                cwd=ctx.output_dir,
            ),
        ]
    )
    preview_lines.extend(
        format_command_preview(
            _rscript_command(
                ctx,
                "supplementary_figures/EDFig4_ChIP/5_program_enrichment_test.R",
                program,
                k,
            ),
            cwd=ctx.output_dir,
        )
        for program in chip_programs
    )
    preview_lines.append(
        format_command_preview(
            _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/6_plot_program_TF_enrichment.R"),
            cwd=ctx.output_dir,
        )
    )

    def _action() -> None:
        if run_collection:
            run_command(
                _bash_command(
                    ctx,
                    "supplementary_figures/EDFig4_ChIP/1_collect_ChIP_K562.sh",
                    _lift_over(ctx),
                    ctx.input_path("liftOver_chain"),
                ),
                cwd=ctx.output_dir,
            )
        run_command(
            _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/2_calculate_TF_chipScore.R"),
            cwd=ctx.output_dir,
        )
        run_command(
            _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/3_ChIPscore_program_enrichment.R"),
            cwd=ctx.output_dir,
        )
        run_command(
            _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/4_select_useful_experiments.R"),
            cwd=ctx.output_dir,
        )
        for program in chip_programs:
            run_command(
                _rscript_command(
                    ctx,
                    "supplementary_figures/EDFig4_ChIP/5_program_enrichment_test.R",
                    program,
                    k,
                ),
                cwd=ctx.output_dir,
            )
        run_command(
            _rscript_command(ctx, "supplementary_figures/EDFig4_ChIP/6_plot_program_TF_enrichment.R"),
            cwd=ctx.output_dir,
        )

    return PlannedStage(spec=stage, preview_lines=tuple(preview_lines), action=_action)


def _build_edfig5_stage(ctx: ResolvedPlotWorkflow) -> PlannedStage:
    stage = ctx.spec.stages[2]
    k = int(ctx.parameter("k", 60))
    growth_programs = _string_list(
        ctx.parameter("growth_screening_programs"),
        "workflows.plot.supplementary.parameters.growth_screening_programs",
    ) or ("P16",)
    posterior_trait_files = _posterior_trait_file_paths(
        ctx,
        "posterior_trait_files",
        "workflows.plot.supplementary.parameters.posterior_trait_files",
    )
    autophagy_geneset_file = _resolve_optional_parameter_path(
        ctx,
        ctx.parameter("autophagy_geneset_file"),
        "workflows.plot.supplementary.parameters.autophagy_geneset_file",
        base_dir=ctx.input_path("geneset_dir"),
        default=ctx.input_path("geneset_dir") / "Autophagosome_genes.txt",
    )

    commands: list[list[str]] = []
    commands.extend(
        _rscript_command(
            ctx,
            "supplementary_figures/EDFig5/1_growth_screening_comparison.R",
            program,
            k,
        )
        for program in growth_programs
    )
    commands.extend(
        [
            _rscript_command(ctx, "supplementary_figures/EDFig5/2_causal_inference_regulators.R", k),
            _rscript_command(ctx, "supplementary_figures/EDFig5/3_cellCycledistribution.R", k),
            _rscript_command(
                ctx,
                "supplementary_figures/EDFig5/4_autophagy_qqlot.R",
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
                autophagy_geneset_file,
            ),
        ]
    )
    return _create_static_stage(stage, ctx.output_dir, tuple(commands))


def _build_multicell_stage(ctx: ResolvedPlotWorkflow) -> PlannedStage:
    stage = ctx.spec.stages[3]
    trait_files = _string_list(
        ctx.parameter("multicell_trait_files"),
        "workflows.plot.supplementary.parameters.multicell_trait_files",
    ) or ("Backman_2021_133.per_gene_estimates.tsv", "Backman_2021_131.per_gene_estimates.tsv")

    commands: list[list[str]] = [
        _rscript_command(ctx, "supplementary_figures/FigureS4-12_multiCell_model/1_regulatory_effects_correlation.R"),
        _rscript_command(ctx, "supplementary_figures/FigureS4-12_multiCell_model/2_overlap_program_genes.R"),
        _rscript_command(ctx, "supplementary_figures/FigureS4-12_multiCell_model/3_make_overview_heatmap_burdenCor.R"),
    ]
    commands.extend(
        _rscript_command(
            ctx,
            "supplementary_figures/FigureS4-12_multiCell_model/4_stepWise_regression.R",
            trait_file,
        )
        for trait_file in trait_files
    )
    commands.extend(
        [
            _rscript_command(ctx, "supplementary_figures/FigureS4-12_multiCell_model/5_stepWise_regression_plot_IGF1.R"),
            _rscript_command(ctx, "supplementary_figures/FigureS4-12_multiCell_model/6_stepWise_regression_plot_HbA1c.R"),
        ]
    )
    return _create_static_stage(stage, ctx.output_dir, tuple(commands))


def _build_cnmf_validation_stage(ctx: ResolvedPlotWorkflow) -> PlannedStage:
    stage = ctx.spec.stages[4]
    permutation_defaults = {
        "program_n": 5,
        "regulator_n": 3,
        "program_top_def": 200,
        "lof_thresh": 0.1,
        "trait": "MCH",
    }
    permutation_jobs = _mapping_list_with_defaults(
        ctx.parameter("figure5_permutation_jobs"),
        "workflows.plot.supplementary.parameters.figure5_permutation_jobs",
        permutation_defaults,
    )
    leave_one_out_defaults = {
        "program_n": 5,
        "regulator_n": 3,
        "trait": "MCH",
        "program_top_def": 200,
    }
    leave_one_out_jobs = _mapping_list_with_defaults(
        ctx.parameter("leave_one_out_jobs"),
        "workflows.plot.supplementary.parameters.leave_one_out_jobs",
        leave_one_out_defaults,
    )
    effect_direction_traits = _string_list(
        ctx.parameter("effect_direction_traits"),
        "workflows.plot.supplementary.parameters.effect_direction_traits",
    ) or ("MCH", "RDW", "IRF")
    effect_direction_program_top_defs = _int_list(
        ctx.parameter("effect_direction_program_top_defs"),
        "workflows.plot.supplementary.parameters.effect_direction_program_top_defs",
    ) or (100, 200)
    k = int(ctx.parameter("k", 60))
    posterior_trait_files = _posterior_trait_file_paths(
        ctx,
        "posterior_trait_files",
        "workflows.plot.supplementary.parameters.posterior_trait_files",
    )
    leave_one_out_posterior = posterior_trait_files.get(
        str(leave_one_out_jobs[0].get("trait", "MCH")),
        posterior_trait_files["MCH"],
    )

    preview_lines: list[str] = []
    preview_lines.extend(
        format_command_preview(
            _rscript_command(
                ctx,
                "Figure5/1_permutation_test_step1.R",
                job["program_n"],
                job["regulator_n"],
                job["program_top_def"],
                job["lof_thresh"],
                job.get("trait", "MCH"),
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
            ),
            cwd=ctx.output_dir,
        )
        for job in permutation_jobs
    )
    for top_n in sorted({int(job["program_top_def"]) for job in permutation_jobs}):
        matching = next(job for job in permutation_jobs if int(job["program_top_def"]) == top_n)
        example_file = f'{matching.get("trait", "MCH")}_program{matching["program_n"]}_regulator{matching["regulator_n"]}_LOF{matching["lof_thresh"]}.txt'
        preview_lines.append(
            format_command_preview(
                _rscript_command(
                    ctx,
                    "Figure5/2_permutation_test_step2.R",
                    top_n,
                    matching.get("trait", "MCH"),
                    example_file,
                ),
                cwd=ctx.output_dir,
            )
        )
    preview_lines.append(
        format_command_preview(
            _bash_command(ctx, "Figure5/3_LeaveOneOut_test_step1.sh", k, leave_one_out_posterior),
            cwd=ctx.output_dir,
        )
    )
    preview_lines.append(
        "Run Figure5/4_LeaveOneOut_test_step2.R for each generated chunk in `loo_topgene_prediction/gene_list/`."
    )
    preview_lines.append(
        format_command_preview(
            _rscript_command(
                ctx,
                "Figure5/5_split_crossValidation.R",
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
            ),
            cwd=ctx.output_dir,
        )
    )
    preview_lines.extend(
        format_command_preview(
            _rscript_command(
                ctx,
                "Figure5/6_predict_effectDirection.R",
                trait,
                top_n,
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
            ),
            cwd=ctx.output_dir,
        )
        for trait in effect_direction_traits
        for top_n in effect_direction_program_top_defs
    )
    preview_lines.append(
        format_command_preview(
            _rscript_command(
                ctx,
                "Figure5/7_crossTraits_effectDirection.R",
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
            ),
            cwd=ctx.output_dir,
        )
    )

    def _action() -> None:
        for job in permutation_jobs:
            run_command(
                _rscript_command(
                    ctx,
                    "Figure5/1_permutation_test_step1.R",
                    job["program_n"],
                    job["regulator_n"],
                    job["program_top_def"],
                    job["lof_thresh"],
                    job.get("trait", "MCH"),
                    k,
                    posterior_trait_files["MCH"],
                    posterior_trait_files["RDW"],
                    posterior_trait_files["IRF"],
                ),
                cwd=ctx.output_dir,
            )
        for top_n in sorted({int(job["program_top_def"]) for job in permutation_jobs}):
            matching = next(job for job in permutation_jobs if int(job["program_top_def"]) == top_n)
            example_file = f'{matching.get("trait", "MCH")}_program{matching["program_n"]}_regulator{matching["regulator_n"]}_LOF{matching["lof_thresh"]}.txt'
            run_command(
                _rscript_command(
                    ctx,
                    "Figure5/2_permutation_test_step2.R",
                    top_n,
                    matching.get("trait", "MCH"),
                    example_file,
                ),
                cwd=ctx.output_dir,
            )
        run_command(
            _bash_command(ctx, "Figure5/3_LeaveOneOut_test_step1.sh", k, leave_one_out_posterior),
            cwd=ctx.output_dir,
        )
        gene_list_dir = ctx.output_dir / "loo_topgene_prediction" / "gene_list"
        chunk_files = sorted(path.name for path in gene_list_dir.iterdir() if path.is_file() and path.name != "allgenes.txt")
        for job in leave_one_out_jobs:
            for chunk_file in chunk_files:
                run_command(
                    _rscript_command(
                        ctx,
                        "Figure5/4_LeaveOneOut_test_step2.R",
                        chunk_file,
                        job["program_n"],
                        job["regulator_n"],
                        job.get("trait", "MCH"),
                        job["program_top_def"],
                        k,
                        posterior_trait_files["MCH"],
                        posterior_trait_files["RDW"],
                        posterior_trait_files["IRF"],
                    ),
                    cwd=ctx.output_dir,
                )
        run_command(
            _rscript_command(
                ctx,
                "Figure5/5_split_crossValidation.R",
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
            ),
            cwd=ctx.output_dir,
        )
        for trait in effect_direction_traits:
            for top_n in effect_direction_program_top_defs:
                run_command(
                    _rscript_command(
                        ctx,
                        "Figure5/6_predict_effectDirection.R",
                        trait,
                        top_n,
                        k,
                        posterior_trait_files["MCH"],
                        posterior_trait_files["RDW"],
                        posterior_trait_files["IRF"],
                    ),
                    cwd=ctx.output_dir,
                )
        run_command(
            _rscript_command(
                ctx,
                "Figure5/7_crossTraits_effectDirection.R",
                k,
                posterior_trait_files["MCH"],
                posterior_trait_files["RDW"],
                posterior_trait_files["IRF"],
            ),
            cwd=ctx.output_dir,
        )

    return PlannedStage(spec=stage, preview_lines=tuple(preview_lines), action=_action)


def _build_supplementary_stage_plans(ctx: ResolvedPlotWorkflow) -> tuple[PlannedStage, ...]:
    return (
        _build_ldsc_stage(ctx),
        _build_chip_stage(ctx),
        _build_edfig5_stage(ctx),
        _build_multicell_stage(ctx),
        _build_cnmf_validation_stage(ctx),
    )


def build_plot_gene_burden_tasks(config: LoadedConfig) -> list[Task]:
    return _build_plot_workflow_tasks(
        config,
        GENE_BURDEN_WORKFLOW,
        _prepare_gene_burden_workspace,
        _build_gene_burden_stage_plans,
    )


def build_plot_perturbseq_gene_level_tasks(config: LoadedConfig) -> list[Task]:
    return _build_plot_workflow_tasks(
        config,
        PERTURBSEQ_GENE_LEVEL_WORKFLOW,
        _prepare_gene_level_workspace,
        _build_gene_level_stage_plans,
    )


def build_plot_perturbseq_cnmf_tasks(config: LoadedConfig) -> list[Task]:
    return _build_plot_workflow_tasks(
        config,
        PERTURBSEQ_CNMF_WORKFLOW,
        _prepare_cnmf_workspace,
        _build_cnmf_stage_plans,
    )


def build_plot_supplementary_tasks(config: LoadedConfig) -> list[Task]:
    return _build_plot_workflow_tasks(
        config,
        SUPPLEMENTARY_WORKFLOW,
        _prepare_supplementary_workspace,
        _build_supplementary_stage_plans,
    )


PLOT_SUBWORKFLOW_BUILDERS = {
    "gene_burden": build_plot_gene_burden_tasks,
    "perturbseq_gene_level": build_plot_perturbseq_gene_level_tasks,
    "perturbseq_cnmf": build_plot_perturbseq_cnmf_tasks,
    "supplementary": build_plot_supplementary_tasks,
}


def build_plot_tasks(config: LoadedConfig) -> list[Task]:
    plot = config.workflow("plot")
    tasks: list[Task] = []
    for key, builder in PLOT_SUBWORKFLOW_BUILDERS.items():
        if plot.get(key):
            tasks.extend(builder(config))
    if not tasks:
        raise ConfigError("workflows.plot must contain at least one configured subworkflow")
    return tasks
