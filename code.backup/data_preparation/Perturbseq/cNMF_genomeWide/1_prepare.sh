#!/bin/bash

K=$1 ##30 60 90 120
module load python/3.9
module load cairo

##step1 normalize and prepare run parameters
cnmf prepare --output-dir ./cNMF  --name cNMF_K$K -c data/Perturbseq/K562_gwps_raw_singlecell_01.h5ad -k $K  --n-iter 100 --seed 14 --numgenes 2000 --total-workers 200
