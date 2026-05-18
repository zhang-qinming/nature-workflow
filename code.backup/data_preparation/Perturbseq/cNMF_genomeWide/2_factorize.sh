#!/bin/bash
module load python/3.9
module load cairo
K=$1 ##30 60 90 120
SEED=$2 ##0-199

cnmf factorize --output-dir ./cNMF --name cNMF_K$K --worker-index $SEED --total-workers 200
