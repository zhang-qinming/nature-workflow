#!/bin/bash
set -euo pipefail

t1=${1:?first trait label is required}
t2=${2:?second trait label is required}
sumstats_a=${3:?first LDSC sumstats file is required}
sumstats_b=${4:?second LDSC sumstats file is required}
lddir=${5:?LD reference directory is required}
ldsc_py=${6:?ldsc.py executable is required}

mkdir -p genetic_correlation

"${ldsc_py}" \
  --rg "${sumstats_a},${sumstats_b}" \
  --ref-ld-chr "${lddir}/" \
  --w-ld-chr "${lddir}/" \
  --out "genetic_correlation/${t1}_${t2}"
