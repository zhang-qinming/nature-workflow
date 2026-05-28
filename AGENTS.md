# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Project overview

Config-driven workflow orchestration for a computational biology paper pipeline (GeneBayes LoF prediction + Perturb-seq cNMF/gene-level analysis + plotting/figures). Python 3.10+ CLI (`paper-pipeline`) reads a YAML config and dispatches shell commands or Python actions across multiple conda environments.

**This is not a git repository.** There is no version control history to reference.

## Setup

See `envs/README.md` for conda environment creation and editable install instructions. Key points:

- `paper-pipeline-r` and `paper-pipeline-plot` do NOT need the Python package installed
- 3 CRAN packages must be installed manually into `paper-pipeline-plot`: `topr`, `ggallin`, `ggarchery`

```bash
conda run -n paper-pipeline-plot Rscript -e "install.packages(c('topr', 'ggallin', 'ggarchery'), repos='https://cloud.r-project.org')"
```

## Commands

```bash
paper-pipeline plan --config configs/paper-pipeline.default.yaml <workflow>    # Print steps
paper-pipeline run  --config configs/paper-pipeline.default.yaml <workflow>    # Execute
paper-pipeline run  --config configs/paper-pipeline.default.yaml <workflow> --dry-run
paper-pipeline run  --config configs/paper-pipeline.default.yaml all           # All top-level workflows
```

See `docs/全流程运行说明.md` for the full workflow inventory and recommended execution order.

## Slurm batch scripts

```bash
BATCH_ROOT=./scripts/run_once PROJECT_ROOT=. bash generate_run_once_batch_scripts.sh
BATCH_ROOT=./scripts/traits   PROJECT_ROOT=. bash generate_trait_batch_scripts.sh
BATCH_ROOT=./scripts/plot     PROJECT_ROOT=. bash generate_plot_batch_scripts.sh
BATCH_ROOT=./scripts/figures  PROJECT_ROOT=. bash generate_figure_sbatch.sh
```

See `docs/一次性阶段Slurm提交说明.md`, `docs/trait_batch_scripts_说明.md`, `docs/plot_batch_scripts_说明.md`.

## Architecture

### Entrypoint and dispatch

- CLI: `src/paper_pipeline/cli.py` → delegates to `src/paper_pipeline/dataprep/cli.py`
- Subcommands: `plan` (print tasks), `run` (execute with optional `--dry-run`)
- `all` resolves to configured top-level aggregates: `genebayes`, `perturbseq`, `plot`, `figures`
- Workflow registry: `src/paper_pipeline/workflows.py` — `WORKFLOW_BUILDERS` dict maps names to builder functions

### Core abstractions (`src/paper_pipeline/dataprep/`)

- **`config.py`**: `load_config(path)` → `LoadedConfig` dataclass. All `inputs.*` / `outputs.*` paths resolve relative to the **config file's directory**. Provides `executable_command()`, `artifact_path()`, `resolve_path_or_artifact()`. Raise `ConfigError` for config problems.
- **`tasks.py`**: `Task` dataclass — `name`, `preview`, and either a `command: list[str]` or an `action: Callable[[], None]`.
- **`runner.py`**: `print_tasks()` / `run_tasks()` — sequential `subprocess.run()` with `check=True`.

Every workflow builder is `(LoadedConfig) -> list[Task]`. The `executables` config section maps logical names to conda `run -n <env>` prefixes.

### Source tree

```
src/paper_pipeline/
  dataprep/          Core config, task system, runner, scientific workflows
    workflows/         GeneBayes + Perturbseq workflow builders
    genebayes/         Feature prediction, merge, training
    perturbseq/        cNMF pipeline, gene-level limma, R scripts (bundled as package data)
  plot/              Legacy paper reproduction plotting (symlink-heavy workspaces)
    workflows.py       Builders: gene_burden, gene_level, cnmf, supplementary
    scripts/           R/shell scripts by figure number
  figures/           New reusable figure workflows (explicit inputs, deterministic outputs)
    workflows.py       Builders: cnmf, burden_volcano, gwas_manhattan, cross_trait, etc.
    scripts/           Standalone R rendering scripts
```

### Design boundary: `plot` vs `figures`

- **`plot-*`**: Legacy paper reproduction. Creates symlink-heavy workspaces under `outputs/plots/`, expects specific file layouts. Uses `code.backup/` as frozen reference only — execution is from `src/paper_pipeline/plot/scripts/`.
- **`figures-*`**: Newer reusable rendering. Takes explicit file/path inputs, no implicit `data/...` assumptions, writes structured outputs under `tables/`, `plots/`, `meta/` with deterministic IDs. **Prefer adding new figure types here.**

### Config resolution rules (critical)

- All `inputs.*` / `outputs.*` relative paths resolve relative to the **config file's directory**, NOT `project_root`
- `project_root` only affects task `cwd` and command preview display
- `artifact_root` defaults to `outputs` (falls back to legacy keys `output_root` / `outputs_root`)
- Default config at `configs/paper-pipeline.default.yaml` sets `project_root: ..` and `artifact_root: ../outputs` because it lives in `configs/`

### cNMF workflow splitting

Both `essential` and `genomewide` cNMF chains support three granularity levels:
1. Full chain — requires coverage of all workers and programs
2. Post shortcut (`*-post` = postbase + association)
3. Fine-grained: `*-postbase-base`, `*-postbase-regulation` (supports program sharding), `*-association`

Recommended cluster sharding order:
`prepare → factorize(shard by worker_indices) → postbase-base → postbase-regulation(shard by program_indices) → association`

### Slurm generator env var overrides

Shell generators use environment variables, not manual edits. Key overrides:

- `LDSC_REPO_DIR`: path to `bulik/ldsc` checkout (for `generate_plot_batch_scripts.sh`)
- `LDSC_PY` / `LIFTOVER_BIN`: override binary paths directly
- `BATCH_ROOT` / `PROJECT_ROOT`: output and project directories

### Root-level files

- `code.backup/` — frozen reference copy of original paper scripts; not executed
- `code.tree` — generated tree snapshot of the repo
- `paper_pipeline.txt` — likely pip freeze or dependency listing

## Further reading

- Workflow inventory and execution order: `docs/全流程运行说明.md`
- Data formats, field specs, and input/output examples: `docs/全流程数据描述_按代码顺序.md`
- Config field reference: `configs/README.md`
- Conda environment setup: `envs/README.md`
- Plotting preparation: `docs/plotting_阶段准备说明.md`

## Code conventions

- Use `ConfigError` for config parse errors
- Workflow builders return `list[Task]`; keep tasks explicit and composable
- Dataclasses for resolved workflow configs
- Shell generators use env var overrides, not manual edits
- `figures` workflows: explicit file/path inputs, deterministic outputs — no implicit `data/...` assumptions
