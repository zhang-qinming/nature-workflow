import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3  
sc.logging.print_header()
adata = sc.read_h5ad("data/Perturbseq/K562_gwps_raw_singlecell_01.h5ad")
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=500)

adata.to_df().to_csv("data/Perturbseq/K562gwps_raw_count.csv")