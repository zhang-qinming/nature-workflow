# `generate_plot_batch_scripts.sh` 说明

这份文档对应仓库根目录下的：

- [generate_plot_batch_scripts.sh](../generate_plot_batch_scripts.sh)

目标是根据当前代码生成 plotting 阶段的 Slurm 脚本与配置文件。

## 1. 默认值与建议

脚本内置默认值是集群绝对路径（维护者环境），不是仓库内相对路径：

- `BATCH_ROOT=/gpfs/.../scripts/plot`
- `CONDA_SH=/gpfs/.../miniconda3/etc/profile.d/conda.sh`
- `LDSC_REPO_DIR=/gpfs/.../software/ldsc`

如果你希望输出到当前仓库，建议显式覆盖：

```bash
BATCH_ROOT=./scripts/plot PROJECT_ROOT=. bash generate_plot_batch_scripts.sh
```

注意：脚本每次运行都会删除并重建 `PLOT_ROOT` 下的：

- `traits/`
- `supplementary/`
- `tools/`

因此不要在这些目录里长期手工修改生成物；需要定制请通过环境变量覆盖，或把生成结果复制到别处再改。
`logs/` 目录不会被清空，历史输出会累积。

## 2. 会生成什么

执行：

```bash
bash generate_plot_batch_scripts.sh
```

会在 `PLOT_ROOT`（默认等于 `BATCH_ROOT`）下生成：

- `traits/trait<N>/`
- `supplementary/`
- `tools/`
- `logs/`
- `submit_plot_workflows.sh`

### 2.1 每个 trait 目录

每个 `trait<N>` 会生成：

- `pipeline.plot.yaml`
- `plot_gene_burden.sh`
- `plot_perturbseq_gene_level.sh`
- `plot_perturbseq_cnmf.sh`
- `submit_trait_plots.sh`

`submit_trait_plots.sh` 会依次提交这 3 个主图作业。
其中 `plot_perturbseq_cnmf.sh` 里的 Figure4 trans-eQTL follow-up 会按 program 并发执行。

### 2.2 supplementary 目录

- `pipeline.plot.yaml`
- `plot_supplementary.sh`
- `submit_supplementary.sh`

默认生成的 supplementary 配置里，Figure1 LDSC 相关参数 `ldsc_sumstats` / `genetic_correlation_pairs` 是空列表。
当前代码会在 `ldsc_sumstats` 为空时自动跳过 Figure1 LDSC stage；如果要启用这一部分，需要手动补齐这些参数。

### 2.3 tools 目录

当 `USE_TOOLS_ENV_WRAPPERS=1`（默认）时：

- `ldsc_py_wrapper.sh`
- `liftOver_wrapper.sh`

## 3. 多 trait 生成

```bash
PLOT_TRAIT_NUMS=86,106,131 bash generate_plot_batch_scripts.sh
```

脚本会按 trait 列表生成 `traits/trait86`、`traits/trait106`、`traits/trait131`。

## 4. 生成后如何提交

### 4.1 提交单个 trait 的某一类图

```bash
cd scripts/plot/traits/trait86
sbatch plot_gene_burden.sh
sbatch plot_perturbseq_gene_level.sh
sbatch plot_perturbseq_cnmf.sh
```

### 4.2 提交单个 trait 的整套主图

```bash
cd scripts/plot/traits/trait86
./submit_trait_plots.sh
```

### 4.3 从根目录批量提交所有 trait

```bash
cd scripts/plot
./submit_plot_workflows.sh
```

默认行为：

- 遍历 `traits/trait*` 下的 `submit_trait_plots.sh`
- 提交每个 trait 的 3 个主图作业
- 默认不提交 supplementary

要把 supplementary 一起提交：

```bash
cd scripts/plot
INCLUDE_SUPPLEMENTARY=1 ./submit_plot_workflows.sh
```

## 5. 常用可调参数

### 5.1 路径和环境

- `BATCH_ROOT`
- `PLOT_ROOT`
- `PLOT_TRAITS_ROOT`
- `PLOT_SUPPLEMENTARY_ROOT`
- `PLOT_TOOLS_ROOT`
- `LOGS_ROOT`
- `PLOT_CONFIG_NAME`
- `CONTROL_ENV`
- `PLOT_ENV`
- `TOOLS_ENV`
- `CONDA_SH`

### 5.2 工具与 wrapper

- `USE_TOOLS_ENV_WRAPPERS`
- `LDSC_REPO_DIR`
- `LDSC_PY`
- `LIFTOVER_BIN`

### 5.3 trait 与图参数

- `PLOT_TRAIT_NUMS`
- `PLOT_TRAIT_NUM`
- `PLOT_TRAIT_FILE`
- `PLOT_RAW_BURDEN_FILE`
- `PLOT_GWAS_TRAIT`
- `PLOT_CROSS_TRAIT_TRAITS`
- `PLOT_CNMF_K`
- `PLOT_CNMF_MODE`
- `PLOT_LABEL_PROGRAMS`
- `PLOT_CORREGULATION_PROGRAM_A`
- `PLOT_CORREGULATION_PROGRAM_B`
- `PLOT_ENABLE_MULTIPLE_REGRESSION`
- `PLOT_ENABLE_TRANS_EQTL`
- `PLOT_TRANS_EQTL_TOP_N`
- `PLOT_TRANS_EQTL_MATCH_TABLE`
- `PLOT_CNMF_TRANS_EQTL_PARALLEL_JOBS`
- `PLOT_MCH_POSTERIOR_FILE`
- `PLOT_RDW_POSTERIOR_FILE`
- `PLOT_IRF_POSTERIOR_FILE`
- `PLOT_AUTOPHAGY_GENESET_FILE`
- `SUPPLEMENTARY_GROWTH_PROGRAMS`

其中：

- `PLOT_CNMF_MODE`
  - `generic` 或 `legacy`
  - 默认 `generic`
  - `generic` 适合批量 catalog LoF trait，默认跳过原 paper 的特例 follow-up
  - `legacy` 才默认启用 paper-specific 的 Figure4 follow-up 逻辑
- `PLOT_ENABLE_MULTIPLE_REGRESSION`
  - 默认 `inherit`
  - `inherit` 时跟随 `PLOT_CNMF_MODE`
  - 控制是否运行 `Figure4/4_multipleRegression.R`
- `PLOT_ENABLE_TRANS_EQTL`
  - 默认 `inherit`
  - `inherit` 时，`generic` 默认关闭，`legacy` 仅在当前 trait 匹配到 GWAS 时开启
  - 控制是否运行 Figure4 trans-eQTL follow-up
  - 只有在当前 trait 有对应 variant-level GWAS summary 时才建议开启
- `PLOT_TRANS_EQTL_MATCH_TABLE`
  - 用于提供 trait -> trans-eQTL GWAS 文件的映射；匹配不到时会跳过 follow-up

### 5.4 Slurm 资源

公共参数：

- `PLOT_PARTITION`
- `PLOT_TIME`

分类参数：

- `PLOT_GENE_BURDEN_MEM`
- `PLOT_GENE_BURDEN_CPUS`
- `PLOT_GENE_LEVEL_MEM`
- `PLOT_GENE_LEVEL_CPUS`
- `PLOT_CNMF_MEM`
- `PLOT_CNMF_CPUS`
- `PLOT_SUPPLEMENTARY_MEM`
- `PLOT_SUPPLEMENTARY_CPUS`
- `PLOT_SUPPLEMENTARY_TIME`

## 6. wrapper 机制

默认会在 `tools/` 中生成 wrapper，并在配置中把：

- `executables.ldsc_py`
- `executables.liftover`

指向 wrapper。

wrapper 内部执行：

```bash
conda run --no-capture-output -n <TOOLS_ENV> python <LDSC_REPO_DIR>/ldsc.py
conda run --no-capture-output -n <TOOLS_ENV> liftOver
```

如果不想用 wrapper：

```bash
USE_TOOLS_ENV_WRAPPERS=0
LDSC_PY=/abs/path/to/ldsc.py
LIFTOVER_BIN=/abs/path/to/liftOver
```

## 7. 运行前提

1. 上游结果已准备
- `plot-gene-burden` 依赖 `genebayes`
- `plot-perturbseq-gene-level` 依赖 `perturbseq-gene-level-summarize` 与 `genebayes`
- `plot-perturbseq-cnmf` 依赖 `perturbseq-cnmf-genomewide-postbase`、`perturbseq-cnmf-genomewide-association` 与 `genebayes`
- `plot-supplementary` 还依赖 `LDSC / TF_ChIP / multi-cell` 输入

2. 环境已安装

```bash
mamba env create -f envs/paper-pipeline-plot.yml
conda run -n paper-pipeline-control pip install -e . --no-deps
```

3. 根据需要补装 CRAN 包

```bash
conda run -n paper-pipeline-plot Rscript -e "install.packages(c('topr', 'ggallin', 'ggarchery'), repos='https://cloud.r-project.org')"
```
