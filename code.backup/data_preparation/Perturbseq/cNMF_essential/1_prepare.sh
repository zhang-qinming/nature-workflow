#!/bin/bash
#SBATCH
cd data/Perturbseq/

CT=$1
module load python/3.9
module load cairo
##step1 normalize and prepare run parameters
cnmf prepare --output-dir ./cNMF/${CT}  --name test1 -c filtered_data/${CT}.h5ad -k 60  --n-iter 100 --seed 14 --numgenes 2000 --total-workers 20
