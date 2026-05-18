from __future__ import annotations

from ..config import ConfigError, LoadedConfig
from ..tasks import Task
from .genebayes import build_genebayes_cpu_tasks, build_genebayes_gpu_tasks, build_genebayes_tasks
from .perturbseq_cnmf import (
    build_perturbseq_cnmf_essential_tasks,
    build_perturbseq_cnmf_essential_prepare_tasks,
    build_perturbseq_cnmf_essential_factorize_tasks,
    build_perturbseq_cnmf_essential_postbase_tasks,
    build_perturbseq_cnmf_essential_postbase_base_tasks,
    build_perturbseq_cnmf_essential_postbase_regulation_tasks,
    build_perturbseq_cnmf_essential_association_tasks,
    build_perturbseq_cnmf_essential_post_tasks,
    build_perturbseq_cnmf_genomewide_tasks,
    build_perturbseq_cnmf_genomewide_prepare_tasks,
    build_perturbseq_cnmf_genomewide_factorize_tasks,
    build_perturbseq_cnmf_genomewide_postbase_tasks,
    build_perturbseq_cnmf_genomewide_postbase_base_tasks,
    build_perturbseq_cnmf_genomewide_postbase_regulation_tasks,
    build_perturbseq_cnmf_genomewide_association_tasks,
    build_perturbseq_cnmf_genomewide_post_tasks,
)
from .perturbseq import (
    build_perturbseq_gene_level_limma_tasks,
    build_perturbseq_gene_level_prepare_tasks,
    build_perturbseq_gene_level_summarize_tasks,
    build_perturbseq_gene_level_tasks,
)


def build_perturbseq_tasks(config: LoadedConfig) -> list[Task]:
    perturbseq = config.workflow("perturbseq")
    tasks: list[Task] = []
    if perturbseq.get("gene_level"):
        tasks.extend(build_perturbseq_gene_level_tasks(config))
    if perturbseq.get("cnmf_essential"):
        tasks.extend(build_perturbseq_cnmf_essential_tasks(config))
    if perturbseq.get("cnmf_genomewide"):
        tasks.extend(build_perturbseq_cnmf_genomewide_tasks(config))
    if not tasks:
        raise ConfigError("workflows.perturbseq must contain at least one configured subworkflow")
    return tasks


__all__ = [
    "build_genebayes_tasks",
    "build_genebayes_cpu_tasks",
    "build_genebayes_gpu_tasks",
    "build_perturbseq_tasks",
    "build_perturbseq_gene_level_tasks",
    "build_perturbseq_gene_level_prepare_tasks",
    "build_perturbseq_gene_level_limma_tasks",
    "build_perturbseq_gene_level_summarize_tasks",
    "build_perturbseq_cnmf_essential_tasks",
    "build_perturbseq_cnmf_essential_prepare_tasks",
    "build_perturbseq_cnmf_essential_factorize_tasks",
    "build_perturbseq_cnmf_essential_postbase_tasks",
    "build_perturbseq_cnmf_essential_postbase_base_tasks",
    "build_perturbseq_cnmf_essential_postbase_regulation_tasks",
    "build_perturbseq_cnmf_essential_association_tasks",
    "build_perturbseq_cnmf_essential_post_tasks",
    "build_perturbseq_cnmf_genomewide_tasks",
    "build_perturbseq_cnmf_genomewide_prepare_tasks",
    "build_perturbseq_cnmf_genomewide_factorize_tasks",
    "build_perturbseq_cnmf_genomewide_postbase_tasks",
    "build_perturbseq_cnmf_genomewide_postbase_base_tasks",
    "build_perturbseq_cnmf_genomewide_postbase_regulation_tasks",
    "build_perturbseq_cnmf_genomewide_association_tasks",
    "build_perturbseq_cnmf_genomewide_post_tasks",
]
