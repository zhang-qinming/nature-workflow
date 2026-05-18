# 一次性阶段 Slurm 提交说明

这份说明对应：

- [generate_run_once_batch_scripts.sh](../generate_run_once_batch_scripts.sh)

当前仓库中没有 `submit_run_once_stages.sh` 总控脚本；一次性阶段通过生成后的 3 个链式提交器来运行。

## 1. 一次性阶段包含什么

### Perturbseq gene-level

- `perturbseq-gene-level-prepare`
- `perturbseq-gene-level-limma`（按 chunk 分片）
- `perturbseq-gene-level-summarize`

### shared cNMF essential

- `perturbseq-cnmf-essential-prepare`
- `perturbseq-cnmf-essential-factorize`（按 worker 分片）
- `perturbseq-cnmf-essential-postbase-base`
- `perturbseq-cnmf-essential-postbase-regulation`（按 program 分片）

### shared cNMF genomewide

- `perturbseq-cnmf-genomewide-prepare`
- `perturbseq-cnmf-genomewide-factorize`（按 worker 分片）
- `perturbseq-cnmf-genomewide-postbase-base`
- `perturbseq-cnmf-genomewide-postbase-regulation`（按 program 分片）

## 2. 生成后的目录与脚本

建议显式指定输出目录到仓库内：

```bash
BATCH_ROOT=./scripts/run_once PROJECT_ROOT=. bash generate_run_once_batch_scripts.sh
```

注意：生成器每次运行都会清空并重建 `shared/` 目录（见脚本中的 `rm -rf`），因此不要在 `scripts/run_once/shared/` 下长期手工修改生成物。
`logs/` 目录不会被清空，历史输出会累积。

生成后：

- `scripts/run_once/shared/`
  - `pipeline.shared.yaml`
  - `pipeline.gene_level.limma_<start>_<end>.yaml`
  - `pipeline.essential.factorize_<start>_<end>.yaml`
  - `pipeline.essential.postbase_regulation_<start>_<end>.yaml`
  - `pipeline.genomewide.factorize_<start>_<end>.yaml`
  - `pipeline.genomewide.postbase_regulation_<start>_<end>.yaml`
  - `perturbseq_gene_level_prepare.sh`
  - `perturbseq_gene_level_limma_<start>_<end>.sh`
  - `perturbseq_gene_level_summarize.sh`
  - `perturbseq_cnmf_essential_prepare.sh`
  - `perturbseq_cnmf_essential_factorize_*.sh`
  - `perturbseq_cnmf_essential_postbase.sh`
  - `perturbseq_cnmf_essential_postbase_base.sh`
  - `perturbseq_cnmf_essential_postbase_regulation_*.sh`
  - `perturbseq_cnmf_genomewide_prepare.sh`
  - `perturbseq_cnmf_genomewide_factorize_*.sh`
  - `perturbseq_cnmf_genomewide_postbase.sh`
  - `perturbseq_cnmf_genomewide_postbase_base.sh`
  - `perturbseq_cnmf_genomewide_postbase_regulation_*.sh`
- `scripts/run_once/`
  - `submit_gene_level_chain.sh`
  - `submit_essential_chain.sh`
  - `submit_genomewide_chain.sh`
  - `logs/`

## 3. 提交方式

```bash
cd scripts/run_once
./submit_gene_level_chain.sh
./submit_essential_chain.sh
./submit_genomewide_chain.sh
```

其中：

- `submit_gene_level_chain.sh`
  - 先提 `prepare`
  - 再并行提所有 limma 分片
  - 最后提 `summarize`
- `submit_essential_chain.sh`
  - `prepare -> factorize shards -> postbase base -> postbase regulation shards`
- `submit_genomewide_chain.sh`
  - `prepare -> factorize shards -> postbase base -> postbase regulation shards`

说明：

- 生成器会同时产出 `*_postbase.sh`（聚合版）和 `*_postbase_base.sh + *_postbase_regulation_*.sh`（细分版）。
- 当前链式提交器默认使用细分版，便于按 program 分片。
- 如果你已经在其他系统完成 program-level regulation 分片，也可以手动改成提交聚合版 `*_postbase.sh`。

## 4. 常见覆盖项

生成前可调：

```bash
GENE_LEVEL_TOTAL_CHUNKS=198 \
GENE_LEVEL_CHUNKS_PER_JOB=25 \
ESSENTIAL_TOTAL_WORKERS=20 \
ESSENTIAL_WORKERS_PER_JOB=2 \
GENOMEWIDE_TOTAL_WORKERS=100 \
GENOMEWIDE_WORKERS_PER_JOB=2 \
BATCH_ROOT=./scripts/run_once PROJECT_ROOT=. \
bash generate_run_once_batch_scripts.sh
```

## 5. 如何接 trait 阶段

一次性阶段提交后，后续可在 trait 目录提交：

```bash
cd scripts/traits/trait86
ESSENTIAL_POSTBASE_JOBID=<jobid> \
GENOMEWIDE_POSTBASE_JOBID=<jobid> \
./submit_trait_chain.sh
```

如果 shared postbase 已完成，可不传这两个环境变量。
