#!/bin/bash

##before this, move files in cNMF/cNMF_K$K/cnmf_tmp to cNMF/cNMF_all/cnmf_tmp and rename files.

module load python/3.9
module load cairo

cnmf prepare --output-dir ./cNMF  --name cNMF_all -c /oak/stanford/groups/pritch/users/mineto/reference/perturb-seq/Weissman/K562_gwps_raw_singlecell_01.h5ad -k 30 60 90 120  --n-iter 100 --seed 14 --numgenes 2000 --total-workers 200
