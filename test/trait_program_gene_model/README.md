# Trait-Program-Gene Model Test Workflow

This directory runs a paper-style Trait-Program-Gene model workflow for each LoF trait.

Inputs are intentionally scattered across existing pipeline outputs and reference folders. Outputs are written only under `OUTPUT_DIR`; the input folders are not modified.

By default this version uses trait-specific automatic model selection:

- `TPGM_PROGRAM_N=auto`: select cNMF programs whose shet-matched burden enrichment passes BH-FDR `TPGM_PROGRAM_FDR_THRESHOLD` (default `0.05`).
- `TPGM_REGULATOR_N=auto`: use `leaps::regsubsets` and BIC to choose the regulator model size, then remove the `shet` covariate and keep the selected cNMF programs. The search is capped by `TPGM_REGULATOR_MAX_N` (default `8`) for runtime and interpretability.
- `TPGM_MAX_GENES_PER_SIDE=all`: do not truncate graph genes in the output tables. Set this to an integer only when you intentionally want a display-size subset.

Set `TPGM_PROGRAM_N` or `TPGM_REGULATOR_N` to an integer to force a fixed number. The sibling directory `trait_program_gene_model_5program_3regulator` keeps the fixed 5+3 wrapper.

## Required Inputs

- `FILE_ID_MAP`: TSV with `id1`, `id2`, `path1`, `path2`. The workflow iterates over `id2`.
- `TPGM_POSTERIOR_DIR`: directory containing `{file_label(path2)}.per_gene_estimates.tsv`.
- `TPGM_GENE_MAP`: two-column ENSG-to-symbol map.
- `TPGM_SHET_PATH`: table with at least `ensg`, `shet_BIN`; `shet` is used if present.
- `TPGM_SPECTRA_PATH`: cNMF gene spectra score matrix, for example `cNMF_all.gene_spectra_score.k_60.dt_0_5.txt`.
- `TPGM_REGULATION_DIR`: directory containing `K60_program{index}_perturb_effects.txt`.

Defaults are provided in `config.base.yaml`, but every input above can be overridden with an environment variable.

## Output Location

Set `OUTPUT_DIR` to control outputs. Default:

```bash
/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs/trait_program_gene_model
```

The workflow writes:

- `tables/{trait}_long.tsv`: graph genes for display; discordant genes are retained with parenthesized `gene_label`.
- `tables/{trait}_concordant_long.tsv`: concordant-only subset of the graph table.
- `tables/{trait}_programs.tsv`: selected programs, similar to the current panel program table.
- `tables/{trait}_gene_predictions.tsv`: all tested graph genes with concordant/discordant flags.
- `tables/{trait}_program_rank.tsv`: program-burden ranking, including raw `P` and BH-FDR `q_value`.
- `tables/{trait}_regulator_coefficients.tsv`: regulator model coefficients.
- `tables/{trait}_permutation.tsv`: permutation summary when `TPGM_PERMUTATION_ITERATIONS > 0`.
- `plots/{trait}.pdf` and `plots/{trait}.png`.

## Typical Run

```bash
OUTPUT_DIR=/path/to/output \
BATCH_ROOT=/path/to/scripts \
FILE_ID_MAP=/path/to/path.file_id_map.tsv \
TPGM_POSTERIOR_DIR=/path/to/genebayes/posterior \
TPGM_GENE_MAP=/path/to/gencode_v41_gname_gid_ALL_sorted_onlyID \
TPGM_SHET_PATH=/path/to/shet_10bins.txt \
TPGM_SPECTRA_PATH=/path/to/cNMF_all.gene_spectra_score.k_60.dt_0_5.txt \
TPGM_REGULATION_DIR=/path/to/cNMF_regulation/K562GW \
bash test/trait_program_gene_model/1_generate.sh

bash test/trait_program_gene_model/2_submit.sh
bash test/trait_program_gene_model/3_check.sh
```
