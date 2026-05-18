# `generate_trait_batch_scripts.sh` 说明

这份说明对应仓库根目录下的：

- [generate_trait_batch_scripts.sh](../generate_trait_batch_scripts.sh)

它负责生成按 trait 运行的脚本（`genebayes-cpu` + `genebayes-gpu` + 两条 association）。

## 1. 默认值与建议

脚本内置默认 `BATCH_ROOT` 是集群绝对路径（维护者环境，`/gpfs/.../scripts/traits`）。

如果希望输出到仓库内，建议：

```bash
BATCH_ROOT=./scripts/traits PROJECT_ROOT=. bash generate_trait_batch_scripts.sh
```

注意：脚本会对每个 `trait<N>/` 目录先清空再生成（见脚本中的 `rm -rf`），重跑会覆盖目录内所有生成物。
当 `GENEBAYES_TRAITS_PER_GPU_JOB > 1` 生成 packed GPU 脚本时，也会清空并重建 `genebayes_packed/`。

## 2. 会生成什么

在 `BATCH_ROOT` 下会生成：

- `trait<N>/`
  - `pipeline.yaml`
  - `genebayes_cpu.sh`
  - `genebayes_gpu.sh`
  - `perturbseq_cnmf_essential_association.sh`
  - `perturbseq_cnmf_genomewide_association.sh`
  - `submit_trait_chain.sh`
- `genebayes_packed/`（可选）
  - `genebayes_pack_<start>_<end>.sh`
- `logs/`

按仓库内推荐路径时即：

- `scripts/traits/trait<N>/...`
- `scripts/traits/genebayes_packed/...`
- `scripts/traits/logs/...`

## 3. 提交顺序建议

先跑一次性阶段，再跑 trait 阶段。

一次性阶段生成与提交见：

- [generate_run_once_batch_scripts.sh](../generate_run_once_batch_scripts.sh)
- [一次性阶段Slurm提交说明.md](./一次性阶段Slurm提交说明.md)

trait 阶段提交示例：

```bash
cd scripts/traits/trait86
./submit_trait_chain.sh
```

## 4. `submit_trait_chain.sh` 的依赖逻辑

该脚本会：

1. 先提交 `genebayes_cpu.sh`
2. CPU 成功后提交 `genebayes_gpu.sh`
3. GPU 成功后再提交
   - `perturbseq_cnmf_essential_association.sh`
   - `perturbseq_cnmf_genomewide_association.sh`

如果 shared postbase 仍在排队，可传：

```bash
ESSENTIAL_POSTBASE_JOBID=<jobid> \
GENOMEWIDE_POSTBASE_JOBID=<jobid> \
./submit_trait_chain.sh
```

如果不传，脚本会默认 shared postbase 已经完成。

## 5. packed GeneBayes 脚本

当 `GENEBAYES_TRAITS_PER_GPU_JOB > 1` 时会生成 `genebayes_packed`。

这些脚本仅用于 `genebayes-gpu`，不会提交 association。
使用前提是对应 trait 的 `genebayes_cpu.sh` 已经完成。

当前默认并发策略：

- A800 节点：并发 3 个 trait（`GENEBAYES_PARALLEL_A800=3`）
- 其他 GPU：并发 2 个 trait（`GENEBAYES_PARALLEL_GPU=2`）

## 6. 常用覆盖项

- trait 范围：`START_TRAIT`、`END_TRAIT`
- 输出目录：`BATCH_ROOT`、`LOGS_ROOT`
- GPU 打包：`GENEBAYES_TRAITS_PER_GPU_JOB`、`GENEBAYES_PACKED_ROOT`
- CPU 环境与资源：`GENEBAYES_CPU_ENV`、`GENEBAYES_CPU_MEM`、`GENEBAYES_CPU_CPUS`、`GENEBAYES_CPU_PARTITION`、`GENEBAYES_CPU_TIME`
- GPU 环境与资源：`GENEBAYES_GPU_ENV`、`GENEBAYES_GPU_MEM`、`GENEBAYES_GPU_CPUS`、`GENEBAYES_GPU_PARTITION`、`GENEBAYES_GPU_GRES`、`GENEBAYES_GPU_TIME`
- 兼容旧变量：`GENEBAYES_ENV`、`GENEBAYES_MEM`、`GENEBAYES_CPUS`、`GENEBAYES_PARTITION`、`GENEBAYES_GRES`、`GENEBAYES_TIME`
- `GENEBAYES_ENV` 现在只作为 `GENEBAYES_CPU_ENV` 的兼容别名；GPU 环境请显式使用 `GENEBAYES_GPU_ENV`
