from __future__ import annotations

from collections.abc import Callable

from .dataprep.config import LoadedConfig
from .dataprep.tasks import Task
from .dataprep.workflows import (
    build_genebayes_cpu_tasks,
    build_genebayes_gpu_tasks,
    build_genebayes_tasks,
    build_perturbseq_cnmf_essential_association_tasks,
    build_perturbseq_cnmf_essential_factorize_tasks,
    build_perturbseq_cnmf_essential_post_tasks,
    build_perturbseq_cnmf_essential_postbase_tasks,
    build_perturbseq_cnmf_essential_postbase_base_tasks,
    build_perturbseq_cnmf_essential_postbase_regulation_tasks,
    build_perturbseq_cnmf_essential_prepare_tasks,
    build_perturbseq_cnmf_essential_tasks,
    build_perturbseq_cnmf_genomewide_association_tasks,
    build_perturbseq_cnmf_genomewide_factorize_tasks,
    build_perturbseq_cnmf_genomewide_post_tasks,
    build_perturbseq_cnmf_genomewide_postbase_tasks,
    build_perturbseq_cnmf_genomewide_postbase_base_tasks,
    build_perturbseq_cnmf_genomewide_postbase_regulation_tasks,
    build_perturbseq_cnmf_genomewide_prepare_tasks,
    build_perturbseq_cnmf_genomewide_tasks,
    build_perturbseq_gene_level_limma_tasks,
    build_perturbseq_gene_level_prepare_tasks,
    build_perturbseq_gene_level_summarize_tasks,
    build_perturbseq_gene_level_tasks,
    build_perturbseq_tasks,
)
from .figures.workflows import (
    build_figures_burden_volcano_tasks,
    build_figures_posterior_volcano_tasks,
    build_figures_cnmf_program_enrichment_tasks,
    build_figures_cnmf_program_top_genes_tasks,
    build_figures_cnmf_tasks,
    build_figures_cross_trait_tasks,
    build_figures_cross_trait_heatmap_tasks,
    build_figures_gene_level_qq_tasks,
    build_figures_gene_level_scatter_tasks,
    build_figures_gwas_locus_zoom_tasks,
    build_figures_gwas_manhattan_tasks,
    build_figures_program_rankings_tasks,
    build_figures_program_heatmap_tasks,
    build_figures_trait_program_gene_panel_tasks,
    build_figures_tasks,
)
from .plot.workflows import (
    build_plot_gene_burden_tasks,
    build_plot_perturbseq_cnmf_tasks,
    build_plot_perturbseq_gene_level_tasks,
    build_plot_supplementary_tasks,
    build_plot_tasks,
)


WorkflowBuilder = Callable[[LoadedConfig], list[Task]]
TOP_LEVEL_WORKFLOWS: tuple[str, ...] = ("genebayes", "perturbseq", "plot", "figures")

WORKFLOW_BUILDERS: dict[str, WorkflowBuilder] = {
    "genebayes": build_genebayes_tasks,
    "genebayes-cpu": build_genebayes_cpu_tasks,
    "genebayes-gpu": build_genebayes_gpu_tasks,
    "perturbseq": build_perturbseq_tasks,
    "perturbseq-gene-level": build_perturbseq_gene_level_tasks,
    "perturbseq-gene-level-prepare": build_perturbseq_gene_level_prepare_tasks,
    "perturbseq-gene-level-limma": build_perturbseq_gene_level_limma_tasks,
    "perturbseq-gene-level-summarize": build_perturbseq_gene_level_summarize_tasks,
    "perturbseq-cnmf-essential": build_perturbseq_cnmf_essential_tasks,
    "perturbseq-cnmf-essential-prepare": build_perturbseq_cnmf_essential_prepare_tasks,
    "perturbseq-cnmf-essential-factorize": build_perturbseq_cnmf_essential_factorize_tasks,
    "perturbseq-cnmf-essential-postbase": build_perturbseq_cnmf_essential_postbase_tasks,
    "perturbseq-cnmf-essential-postbase-base": build_perturbseq_cnmf_essential_postbase_base_tasks,
    "perturbseq-cnmf-essential-postbase-regulation": build_perturbseq_cnmf_essential_postbase_regulation_tasks,
    "perturbseq-cnmf-essential-association": build_perturbseq_cnmf_essential_association_tasks,
    "perturbseq-cnmf-essential-post": build_perturbseq_cnmf_essential_post_tasks,
    "perturbseq-cnmf-genomewide": build_perturbseq_cnmf_genomewide_tasks,
    "perturbseq-cnmf-genomewide-prepare": build_perturbseq_cnmf_genomewide_prepare_tasks,
    "perturbseq-cnmf-genomewide-factorize": build_perturbseq_cnmf_genomewide_factorize_tasks,
    "perturbseq-cnmf-genomewide-postbase": build_perturbseq_cnmf_genomewide_postbase_tasks,
    "perturbseq-cnmf-genomewide-postbase-base": build_perturbseq_cnmf_genomewide_postbase_base_tasks,
    "perturbseq-cnmf-genomewide-postbase-regulation": build_perturbseq_cnmf_genomewide_postbase_regulation_tasks,
    "perturbseq-cnmf-genomewide-association": build_perturbseq_cnmf_genomewide_association_tasks,
    "perturbseq-cnmf-genomewide-post": build_perturbseq_cnmf_genomewide_post_tasks,
    "plot": build_plot_tasks,
    "plot-gene-burden": build_plot_gene_burden_tasks,
    "plot-perturbseq-gene-level": build_plot_perturbseq_gene_level_tasks,
    "plot-perturbseq-cnmf": build_plot_perturbseq_cnmf_tasks,
    "plot-supplementary": build_plot_supplementary_tasks,
    "figures": build_figures_tasks,
    "figures-cnmf": build_figures_cnmf_tasks,
    "figures-cnmf-program-top-genes": build_figures_cnmf_program_top_genes_tasks,
    "figures-cnmf-program-enrichment": build_figures_cnmf_program_enrichment_tasks,
    "figures-burden-volcano": build_figures_burden_volcano_tasks,
    "figures-posterior-volcano": build_figures_posterior_volcano_tasks,
    "figures-cross-trait": build_figures_cross_trait_tasks,
    "figures-cross-trait-heatmap": build_figures_cross_trait_heatmap_tasks,
    "figures-gene-level-qq": build_figures_gene_level_qq_tasks,
    "figures-gene-level-scatter": build_figures_gene_level_scatter_tasks,
    "figures-gwas-locus-zoom": build_figures_gwas_locus_zoom_tasks,
    "figures-gwas-manhattan": build_figures_gwas_manhattan_tasks,
    "figures-program-rankings": build_figures_program_rankings_tasks,
    "figures-program-heatmap": build_figures_program_heatmap_tasks,
    "figures-trait-program-gene-panel": build_figures_trait_program_gene_panel_tasks,
}


__all__ = ["TOP_LEVEL_WORKFLOWS", "WORKFLOW_BUILDERS", "WorkflowBuilder"]
