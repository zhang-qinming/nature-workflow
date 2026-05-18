#!/bin/bash
module load python/3.9
module load cairo

K=$1 #30 60 90 120
T=$2 #0.5

cnmf consensus --output-dir ./cNMF --name cNMF_all --components $K --local-density-threshold $T --show-clustering
