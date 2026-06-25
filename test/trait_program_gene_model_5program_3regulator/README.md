# Trait-Program-Gene Model, 5 Program + 3 Regulator

This is a fixed-default wrapper around `test/trait_program_gene_model`.

It uses the same paper-style data processing, inputs, and output schema, but defaults to:

- `TPGM_PROGRAM_N=5`
- `TPGM_REGULATOR_N=3`
- `TASK_NAME=trait_program_gene_model_5program_3regulator`
- `JOB_NAME=tpgm5p3r`
- `OUTPUT_DIR=${OUTPUT_ROOT}/trait_program_gene_model_5program_3regulator`
- `BATCH_ROOT=/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/trait_program_gene_model_5program_3regulator`

Inputs remain scattered across the same source locations as the base model workflow and can be overridden with:

```bash
FILE_ID_MAP=/path/to/path.file_id_map.tsv
TPGM_POSTERIOR_DIR=/path/to/genebayes/posterior
TPGM_GENE_MAP=/path/to/gencode_v41_gname_gid_ALL_sorted_onlyID
TPGM_SHET_PATH=/path/to/shet_10bins.txt
TPGM_SPECTRA_PATH=/path/to/cNMF_all.gene_spectra_score.k_60.dt_0_5.txt
TPGM_REGULATION_DIR=/path/to/cNMF_regulation/K562GW
```

Typical run:

```bash
bash test/trait_program_gene_model_5program_3regulator/1_generate.sh
bash test/trait_program_gene_model_5program_3regulator/2_submit.sh
bash test/trait_program_gene_model_5program_3regulator/3_check.sh
```
