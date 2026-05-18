#!/bin/bash

GWAS=$1 ##loop over GWAS
ODIR="Partitioning/EUR/${GWAS}"

DATA_list=$( ls ANNOTATION_EUR ) ##annotations 

weights="path to weights reference"
frqfile="path to frequency reference"
baseline="path to baseline reference"

mkdir -p $ODIR

for DATA in $DATA_list; do
 ~/software/ldsc/ldsc.py   \
   --h2          $IFILE    \
   --ref-ld-chr  ANNOTATION_EUR/${DATA}/chr@,$baseline \
   --frqfile-chr $frqfile  \
   --w-ld-chr    $weights  \
   --overlap-annot     --print-cov  --print-coefficients  --print-delete-vals    \
   --out         $ODIR/${DATA}
done
