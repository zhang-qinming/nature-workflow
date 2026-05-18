#!/bin/bash
module load python/3.9
module load cairo

cnmf combine --output-dir ./cNMF --name cNMF_all

cnmf k_selection_plot --output-dir ./cNMF --name cNMF_all

