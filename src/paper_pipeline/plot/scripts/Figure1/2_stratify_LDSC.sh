#!/bin/bash
set -euo pipefail

gwas_sumstats=${1:?LDSC sumstats file is required}
label=${2:?output label is required}
weights=${3:?weights reference prefix is required}
frqfile=${4:?frequency reference prefix is required}
baseline=${5:?baseline reference prefix is required}
ldsc_py=${6:?ldsc.py executable is required}

odir="Partitioning/EUR/${label}"
data_list=$(ls ANNOTATION_EUR)

mkdir -p "$odir"

for data in $data_list; do
  "${ldsc_py}" \
    --h2 "$gwas_sumstats" \
    --ref-ld-chr "ANNOTATION_EUR/${data}/chr@,${baseline}" \
    --frqfile-chr "$frqfile" \
    --w-ld-chr "$weights" \
    --overlap-annot --print-cov --print-coefficients --print-delete-vals \
    --out "${odir}/${data}"
done
