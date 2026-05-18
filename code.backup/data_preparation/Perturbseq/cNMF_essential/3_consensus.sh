#!/bin/bash
#SBATCH 
module load python/3.9
module load cairo
cd data/Perturbseq/

CT=$1

cnmf combine --output-dir ./cNMF/${CT} --name test1

cnmf consensus --output-dir ./cNMF/${CT}  --name test1 --components 60 --local-density-threshold 0.4 --show-clustering
