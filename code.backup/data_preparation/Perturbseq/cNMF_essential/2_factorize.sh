#!/bin/bash
#SBATCH

##run with array jobs. --array=0-19
cd data/Perturbseq

module load python/3.9
module load cairo
SEED=${SLURM_ARRAY_TASK_ID}
CT=$1

cnmf factorize --output-dir ./cNMF/${CT}  --name test1 --worker-index $SEED --total-workers 20
