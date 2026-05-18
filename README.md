# paper-pipeline

基于 YAML 配置的 workflow 入口与编排，用于运行论文流水线（GeneBayes + Perturb-seq dataprep + plotting/figures）。

## 文档入口

- 环境与安装：[envs/README.md](envs/README.md)
- 配置字段语义：[configs/README.md](configs/README.md)
- Workflow 清单与推荐顺序：[docs/全流程运行说明.md](docs/全流程运行说明.md)
- 输入/中间件/输出与格式约束：[docs/全流程数据描述_按代码顺序.md](docs/全流程数据描述_按代码顺序.md)
- Plotting 准备说明：[docs/plotting_阶段准备说明.md](docs/plotting_阶段准备说明.md)
- Slurm 生成与提交：
  - [docs/一次性阶段Slurm提交说明.md](docs/一次性阶段Slurm提交说明.md)
  - [docs/trait_batch_scripts_说明.md](docs/trait_batch_scripts_说明.md)
  - [docs/plot_batch_scripts_说明.md](docs/plot_batch_scripts_说明.md)

## 常用命令

```bash
paper-pipeline plan --config configs/paper-pipeline.default.yaml <workflow>
paper-pipeline run  --config configs/paper-pipeline.default.yaml <workflow>
paper-pipeline run  --config configs/paper-pipeline.default.yaml <workflow> --dry-run
paper-pipeline run  --config configs/paper-pipeline.default.yaml all
```

`<workflow>` 可为单条 workflow（如 `genebayes`、`perturbseq-gene-level-summarize`、`plot-gene-burden`），也可为聚合入口（`perturbseq`、`plot`、`figures`），以及 `all`（展开配置中定义的顶层聚合：`genebayes/perturbseq/plot/figures`）。

## 生成 Slurm 脚本（集群）

```bash
BATCH_ROOT=./scripts/run_once PROJECT_ROOT=. bash generate_run_once_batch_scripts.sh
BATCH_ROOT=./scripts/traits   PROJECT_ROOT=. bash generate_trait_batch_scripts.sh
BATCH_ROOT=./scripts/plot     PROJECT_ROOT=. bash generate_plot_batch_scripts.sh
BATCH_ROOT=./scripts/figures  PROJECT_ROOT=. bash generate_figure_sbatch.sh
```

这些生成器会覆盖/重建输出目录下的生成物；如需定制请通过环境变量覆盖。

## plotting 备注

plotting workflow 会在 `outputs/plots/...` 下创建 legacy 兼容工作区并挂接输入后运行包内脚本。若运行环境无法创建 symlink 而退化为拷贝目录，重跑前建议先清理对应工作区（见 [docs/plotting_阶段准备说明.md](docs/plotting_阶段准备说明.md)）。
