# plotting 阶段准备说明

实现入口：`src/paper_pipeline/plot/workflows.py`

## 当前状态

plotting 已是可执行 workflow。四条入口：`plot-gene-burden`、`plot-perturbseq-gene-level`、`plot-perturbseq-cnmf`、`plot-supplementary`（以及聚合入口 `plot`）。

执行流程：
1. 在 `outputs/plots/...` 下创建 legacy 兼容工作区
2. 通过软链接/复制挂接输入数据
3. 顺序执行 `src/paper_pipeline/plot/scripts/...` 包内脚本
4. 产出 `workflow_manifest.md` 和 `stage_inventory.tsv`

注意：`code.backup/` 仅作参考，运行入口已切到包内脚本。若因文件系统限制退化为拷贝目录，重跑前建议先删除对应工作区再重跑。

## 四条 workflow 的脚本与依赖

| workflow | 脚本范围 | 主要上游依赖 |
|----------|----------|-------------|
| `plot-gene-burden` | Figure2/* | GeneBayes posterior, raw burden, GWAS, geneset, LD reference |
| `plot-perturbseq-gene-level` | Figure3/* | gene-level limma 汇总, GeneBayes posterior, GWAS, shet |
| `plot-perturbseq-cnmf` | Figure4/* | cNMF genomewide 产物, GeneBayes posterior, GWAS, trans-eQTL |
| `plot-supplementary` | Figure1,5 + EDFig4/5 + multi-cell | LDSC, TF ChIP, K562GW, multi-cell, GeneBayes posterior |

关键参数化覆盖（可通过 config 参数控制，不再硬编码）：

- Figure2：`gwas_traits`, `manhattan_gwas_trait`, `lof_enrichment_pairs`, `cross_trait_traits`
- Figure3：`trait_files`, `enrichment_*`, `render_trait_file`
- Figure4：`k`, `trait_files`, `plot_label_programs`, `corregulation_pairs`, `trans_eqtl_*`, `mode`

`mode: generic` 适合批量 catalog LoF trait（默认禁用 paper-specific follow-up），`mode: legacy` 启用原 paper 的 Figure4 follow-up。

## 环境与工具

- `paper-pipeline-control`：运行 CLI
- `paper-pipeline-plot`：运行 R/bash 脚本
- 外部命令：`ldsc_py` 和 `liftover`（需按部署环境配置为真实路径或 wrapper）

手工补装 CRAN 包（当前仍需要）：
```bash
conda run -n paper-pipeline-plot Rscript -e "install.packages(c('topr','ggallin','ggarchery'), repos='https://cloud.r-project.org')"
```

## 配置最小检查

- `executables.plot_bash` / `executables.plot_rscript` 可用
- `executables.ldsc_py` / `executables.liftover` 路径正确
- `workflows.plot.*.inputs` 路径可达
- `plot-supplementary` 的 LDSC 与 ChIP 输入完整；不需要 ChIP 下载链时保持 `run_chip_collection: false`

## 常用命令

```bash
paper-pipeline plan --config config.yaml plot-gene-burden     # 先看命令展开
paper-pipeline run  --config config.yaml plot-gene-burden     # 实际执行
paper-pipeline run  --config config.yaml plot                # 展开配置中所有已启用的 plot 子项
```

## 相关文档

- 配置参考：`configs/README.md`
- Workflow 顺序：`docs/全流程运行说明.md`
- Slurm 批处理：`docs/plot_batch_scripts_说明.md`
