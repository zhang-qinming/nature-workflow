#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/../trait_program_gene_model"

export TASK_NAME="${TASK_NAME:-trait_program_gene_model_5program_3regulator}"
export JOB_NAME="${JOB_NAME:-tpgm5p3r}"
export OUTPUT_ROOT="${OUTPUT_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/outputs}"
export OUTPUT_DIR="${OUTPUT_DIR:-${OUTPUT_ROOT}/trait_program_gene_model_5program_3regulator}"
export RUN_ROOT="${RUN_ROOT:-${OUTPUT_ROOT}/${TASK_NAME}}"
export BATCH_ROOT="${BATCH_ROOT:-/gpfs/chencao/qinminzhang/workflow/catalog_lof/figure_all/scripts/${TASK_NAME}}"
export TPGM_PROGRAM_N=5
export TPGM_REGULATOR_N=3

bash "${BASE_DIR}/4_rerun.sh" "$@"
echo "Check fixed 5+3 wrapper: bash test/trait_program_gene_model_5program_3regulator/3_check.sh"
