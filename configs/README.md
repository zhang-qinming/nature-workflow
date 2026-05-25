# 配置说明

建议从 `configs/paper-pipeline.default.yaml` 或 `configs/paper-pipeline.example.yaml` 复制一份自己的 `config.yaml`，再按机器环境和真实数据路径修改。

## 路径解析规则

`project_root` 只影响两件事：
- workflow 子命令的 `cwd`
- `plan` 输出里的命令预览

**关键规则：`inputs.*` / `outputs.*` 的相对路径都相对于配置文件所在目录解析，与 `project_root` 无关。**

`artifact_root` 兼容旧字段名 `output_root` / `outputs_root`（优先级：`artifact_root` > `output_root` > `outputs_root`，默认 `outputs`）。

如果配置文件放在 `configs/` 目录（与仓库模板一致）：

```yaml
project_root: ..
artifact_root: ../outputs
```

如果配置文件放在仓库根目录（例如 `./config.yaml`）：

```yaml
project_root: .
artifact_root: outputs
```

注意：仓库模板配置位于 `configs/`，因此凡是指向仓库根目录 `data/` 的路径都会写成 `../data/...`。若把模板复制到仓库根目录运行，请改成 `data/...` 或绝对路径。

## CLI workflow 名称

支持三种运行方式：
- 具体 workflow：`genebayes`、`perturbseq-gene-level-limma`、`plot-gene-burden` 等
- 聚合 workflow：`perturbseq`、`plot`、`figures`（只展开配置中存在的子项）
- `all`：展开配置里定义的顶层聚合 workflow（`genebayes`、`perturbseq`、`plot`、`figures`）

`all` 不会重复展开子 workflow。如果配置文件里没有定义任何顶层聚合 workflow，`all` 会报错。

## `executables`

`executables` 描述 workflow 运行时使用的命令前缀：

```yaml
executables:
  plot_bash:
    - conda
    - run
    - --no-capture-output
    - -n
    - paper-pipeline-plot
    - bash
  plot_rscript:
    - conda
    - run
    - --no-capture-output
    - -n
    - paper-pipeline-plot
    - Rscript
  ldsc_py: ldsc.py
  liftover: liftOver
```

`ldsc_py` 和 `liftover` 按外部命令处理；批量生成的 plot 配置会覆盖成 wrapper，手动运行时请按环境改成绝对路径或 wrapper。

## cNMF workflow 拆分与索引

`perturbseq-cnmf-essential` 与 `perturbseq-cnmf-genomewide` 都支持三层粒度：

1. 完整链（不带后缀）
2. 后处理快捷链：`*-post`（等于 `postbase + association`）
3. 最细粒度：`*-postbase-base`、`*-postbase-regulation`、`*-association`

硬约束：

- `worker_indices`：完整链要求覆盖所有 worker
- `program_indices`：`*-post`、`*-postbase`、`*-association` 要求覆盖所有 program；只有 `*-postbase-regulation` 支持按 program 子集分片

推荐做法：
- 默认不写 `worker_indices` / `program_indices`，走全覆盖默认值
- worker 分片：在不同配置中写不同 `worker_indices`，只运行 `*-factorize`
- program 分片：在不同配置中写不同 `program_indices`，只运行 `*-postbase-regulation`
- 汇总：先 `*-postbase-base`，最后 `*-association`（全覆盖）

## plotting workflow 配置

plotting 按数据来源拆成四条 workflow。每条都会创建 legacy 兼容工作区、挂接输入、执行包内脚本，产出 `workflow_manifest.md` 和 `stage_inventory.tsv`。

### `plot-gene-burden`（Figure2）

```yaml
workflows:
  plot:
    gene_burden:
      inputs:
        ld_reference_dir: ../data/GWAS/reference    # Figure2 clumping 需要
```

### `plot-perturbseq-gene-level`（Figure3）

```yaml
parameters:
  trait_files: [Backman_2021_86.per_gene_estimates.tsv]
  enrichment_gwas_trait: 30050
  enrichment_comp: closest
  enrichment_target_gene: HBA1
  enrichment_lof_burden_file: Backman_2021_86_M1_001.summary_statistics.csv
  render_trait_file: Backman_2021_86.per_gene_estimates.tsv
```

### `plot-perturbseq-cnmf`（Figure4）

```yaml
parameters:
  mode: generic          # generic（批量 catalog trait）或 legacy
  k: 60
  trait_files: [Backman_2021_86.per_gene_estimates.tsv]
  run_multiple_regression: false
  run_trans_eqtl_follow_up: false
  plot_label_programs: [4, 16, 25, 40]
  corregulation_pairs:
    - program_a: P25
      program_b: P16
```

`mode: generic` 时默认禁用 paper-specific follow-up（multiple regression / trans-eQTL），适合批量 catalog LoF trait。

### `plot-supplementary`（Figure1、Figure5、EDFig4/5、multi-cell）

需要额外输入：LDSC reference、ChIP chain file、multi-cell 根目录等。

```yaml
workflows:
  plot:
    supplementary:
      inputs:
        ldsc_ld_reference_dir: ../data/LDSC/ref/1000G_EUR_Phase3_plink
        ldsc_baseline_dir: ../data/LDSC/ref/baseline_v1.2
        ldsc_weights_dir: ../data/LDSC/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
        ldsc_frq_dir: ../data/LDSC/ref/1000G_Phase3_frq/1000G.EUR.QC.
        ldsc_sumstats_dir: ../data/LDSC/GWAS
        multicell_cnmf_root: ../data/Perturbseq/cNMF
        multicell_cnmf_regulation_dir: ../data/Perturbseq/cNMF_regulation
        multicell_trait_association_dir: ../data/Perturbseq/trait_association
        liftOver_chain: ../data/hg38ToHg19.over.chain.gz
```

LDSC 开关约定：`ldsc_sumstats` 为空时跳过整个 Figure1 LDSC stage；`genetic_correlation_pairs` 为空时只跳过 genetic correlation 子步骤。

## figures workflow 配置

`figures` 是新的可复用画图主线，直接消费 pipeline 产物或表格输入，不依赖 legacy `data/...` workspace。默认配置未启用 `figures`，需要参考 `configs/paper-pipeline.example.yaml` 添加 `workflows.figures`。

所有 figures workflow 共用 `file_id_map`（4 列 TSV：`id1`、`id2`、`path1`、`path2`，分别对应 GWAS 与 LoF）。示例见 `configs/example_file_id_map.tsv`。

完整配置示例：

```yaml
workflows:
  figures:
    cnmf:
      parameters:
        k: 60
        trait_targets:
          - trait_file: Backman_2021_86.per_gene_estimates.tsv
            trait_id: GCST90081631
        plot_label_programs: [4, 16, 25, 40]
        corregulation_pairs: [{program_a: P25, program_b: P16}]

    burden_volcano:
      inputs:
        file_id_map: example_file_id_map.tsv
        geneset_dir: ../data/geneset
        gene_map: ../data/gencode_v41_gname_gid_ALL_sorted_onlyID
      parameters:
        lof_ids: [GCST90081631]

    gwas_manhattan:
      inputs:
        file_id_map: example_file_id_map.tsv
        gene_annotation: ../data/GWAS/genes.protein_coding.v39.gtf
        geneset_dir: ../data/geneset
      parameters:
        gwas_ids: [PA00638]
        sampling_trigger_rows: 100000
        sampling_base_points: 50000
        sampling_fraction: 0.01
        sampling_max_points: 300000
        sampling_seed: 1

    program_heatmap:
      parameters:
        k: 60
        trait_targets:
          - trait_file: Backman_2021_86.per_gene_estimates.tsv
            trait_id: GCST90081631
        metrics: [program_score, regulator_score]

    gene_level_scatter:
      inputs:
        file_id_map: example_file_id_map.tsv
        correlation_dir: ../outputs/plots/perturbseq_gene_level/data/Perturbseq/trait_association/K562GW/GeneLevel
      parameters:
        lof_ids: [GCST90081631]
        highlight_genes: [HBA1]

    cross_trait:
      inputs:
        file_id_map: example_file_id_map.tsv
      parameters:
        lof_pairs:
          - id_x: GCST90081631
            id_y: GCST90084301
            output_id: GCST90081631__GCST90084301

    gene_level_qq:
      inputs:
        file_id_map: example_file_id_map.tsv
        correlation_dir: ../outputs/plots/perturbseq_gene_level/data/Perturbseq/trait_association/K562GW/GeneLevel
      parameters:
        lof_ids: [GCST90081631]
        output_id: GCST90081631

    program_rankings:
      parameters:
        k: 60
        trait_targets:
          - trait_file: Backman_2021_86.per_gene_estimates.tsv
            trait_id: GCST90081631
        top_n: 12

    gwas_locus_zoom:
      inputs:
        file_id_map: example_file_id_map.tsv
        gene_annotation: ../data/GWAS/genes.protein_coding.v39.gtf
      parameters:
        locus_targets:
          - gwas_id: PA00638
            chrom: "1"
            start: 1000000
            end: 1100000
            flank_bp: 250000
            output_id: PA00638__chr1_1000000_1100000

    cross_trait_heatmap:
      inputs:
        file_id_map: example_file_id_map.tsv
      parameters:
        lof_ids: [GCST90081631, GCST90084301]
        method: pearson
```

各 workflow 的输出统一写在 `artifact_root/figures/<workflow>/` 下，按 `tables/`、`plots/`、`meta/` 分类。

## 仍需在 Linux 上确认的内容

- `ldsc_py` 和 `liftover` 的真实可执行路径
- Figure1 的 `ldsc_annotation_beds` / `ldsc_sumstats` / `genetic_correlation_pairs`
- EDFig4 的 `run_chip_collection: true` 是否符合部署方式
- multi-cell supplementary 的根目录路径

## 相关文档

- 环境与安装：`envs/README.md`
- Workflow 清单与推荐顺序：`docs/全流程运行说明.md`
- 数据格式与字段规范：`docs/全流程数据描述_按代码顺序.md`
- Plotting 准备说明：`docs/plotting_阶段准备说明.md`
- Slurm 提交：`docs/一次性阶段Slurm提交说明.md`、`docs/trait_batch_scripts_说明.md`、`docs/plot_batch_scripts_说明.md`
