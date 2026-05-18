# 环境说明

项目继续使用多套 conda 环境，而不是把所有依赖堆进同一个环境。

当前分工：
- `paper-pipeline-control`
  - 运行主 CLI `paper-pipeline`
  - 负责读取配置并调度各条 workflow
- `paper-pipeline-genebayes`
  - GeneBayes CPU 环境
- `paper-pipeline-genebayes-gpu`
  - GeneBayes GPU 环境
- `paper-pipeline-perturbseq`
  - Perturbseq Python 环境
- `paper-pipeline-r`
  - dataprep 阶段的 R 环境
- `paper-pipeline-plot`
  - plotting 阶段的 R + bash + 外部工具环境
- `paper-pipeline-tools`
  - plotting shell 链里单独调用的外部命令环境
  - 当前用于 `bulik/ldsc` 的 `ldsc.py` 和 `liftOver`

## 创建环境

```bash
mamba env create -f envs/paper-pipeline-control.yml
mamba env create -f envs/paper-pipeline-genebayes.yml
mamba env create -f envs/paper-pipeline-genebayes-gpu.yml
mamba env create -f envs/paper-pipeline-perturbseq.yml
mamba env create -f envs/paper-pipeline-r.yml
mamba env create -f envs/paper-pipeline-plot.yml
mamba env create -f envs/paper-pipeline-tools.yml
```

更新环境：

```bash
mamba env update -f envs/paper-pipeline-control.yml --prune
mamba env update -f envs/paper-pipeline-genebayes.yml --prune
mamba env update -f envs/paper-pipeline-genebayes-gpu.yml --prune
mamba env update -f envs/paper-pipeline-perturbseq.yml --prune
mamba env update -f envs/paper-pipeline-r.yml --prune
mamba env update -f envs/paper-pipeline-plot.yml --prune
mamba env update -f envs/paper-pipeline-tools.yml --prune
```

## 安装当前项目

需要安装 editable package 的环境：

```bash
conda run -n paper-pipeline-control pip install -e . --no-deps
conda run -n paper-pipeline-genebayes pip install -e . --no-deps
conda run -n paper-pipeline-genebayes-gpu pip install -e . --no-deps
conda run -n paper-pipeline-perturbseq pip install -e . --no-deps

```

通常不需要把当前 Python 包安装进：
- `paper-pipeline-r`
- `paper-pipeline-plot`

原因是：
- `paper-pipeline-r` 只负责跑 dataprep 的 R 脚本
- `paper-pipeline-plot` 只负责执行 `plot_bash` / `plot_rscript`
- plotting 的实际脚本路径由 control 环境里安装的 Python 包解析到 `paper_pipeline/plot/scripts/...`

## plotting 环境补充

`envs/paper-pipeline-plot.yml` 当前主要包含：
- plotting 常用 R 包
- `ggrastr`
- `scattermore`
- `ggbreak`
- `MASS`
- `bedtools`
- `plink`
- `curl`

以下两类工具现在推荐单独放进：
- `paper-pipeline-tools`
  - `bulik/ldsc` checkout 所需的 Python 2 依赖
  - `liftOver`

当前通过配置指定：

```yaml
executables:
  ldsc_py: /abs/path/to/scripts/plot/tools/ldsc_py_wrapper.sh
  liftover: /abs/path/to/scripts/plot/tools/liftOver_wrapper.sh
```

`generate_plot_batch_scripts.sh` 现在默认会在 `scripts/plot/tools/` 下生成这两个 wrapper，
内部通过：

```bash
conda run -n paper-pipeline-tools python /abs/path/to/bulik/ldsc/ldsc.py
conda run -n paper-pipeline-tools liftOver
```

其中 `/abs/path/to/bulik/ldsc` 由 `LDSC_REPO_DIR` 控制。生成 plotting 脚本前需要先准备官方仓库 checkout，例如：
```bash
git clone https://github.com/bulik/ldsc /path/to/bulik-ldsc
LDSC_REPO_DIR=/path/to/bulik-ldsc bash generate_plot_batch_scripts.sh
```

如果你不想用 wrapper，也可以把：
- `LDSC_PY`
- `LIFTOVER_BIN`

覆盖成你机器上的真实绝对路径。

以下 3 个 plotting 依赖目前不在 `conda-forge` / `bioconda`，因此没有写进默认
`paper-pipeline-plot.yml`，需要在 Linux 上手工补装 CRAN 包：
- `topr`
  - 影响 `plot-gene-burden`
- `ggallin`
  - 影响 `plot-perturbseq-cnmf`
- `ggarchery`
  - 影响 `plot-supplementary`

可选的手工补包方式：

```bash
conda run -n paper-pipeline-plot Rscript -e "install.packages(c('topr', 'ggallin', 'ggarchery'), repos='https://cloud.r-project.org')"
```

## Linux 上线前最小检查

建议至少确认：

```bash
conda run -n paper-pipeline-control pip install -e . --no-deps
paper-pipeline plan --config configs/paper-pipeline.default.yaml all
paper-pipeline plan --config configs/paper-pipeline.default.yaml plot-gene-burden
paper-pipeline plan --config configs/paper-pipeline.default.yaml plot-perturbseq-gene-level
paper-pipeline plan --config configs/paper-pipeline.default.yaml plot-perturbseq-cnmf
paper-pipeline plan --config configs/paper-pipeline.default.yaml plot-supplementary
paper-pipeline run --config configs/paper-pipeline.default.yaml plot-gene-burden --dry-run
```

如果要启用 Figure1 或 EDFig4 的 shell 链，再额外确认：
- `executables.ldsc_py`
- `executables.liftover`
- `workflows.plot.supplementary.inputs.ldsc_*`
- `workflows.plot.supplementary.inputs.liftOver_chain`

如果你想看从环境创建到整套 workflow 推荐执行顺序，见：
- `docs/全流程运行说明.md`
