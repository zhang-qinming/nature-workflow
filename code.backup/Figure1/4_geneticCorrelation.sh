#!/bin/bash

T1=$1 ##TRAIT1
T2=$2 ##TRAIT2

dd="path to GWAS files"
lddir="path to LD reference"

 ~/software/ldsc/ldsc.py   \
--rg ${dd}/${T1}_irnt.sumstats.gz,${dd}/${T2}_irnt.sumstats.gz \
--ref-ld-chr ${lddir}/ \
--w-ld-chr ${lddir}/ \
--out genetic_correlation/${T1}_${T2}

