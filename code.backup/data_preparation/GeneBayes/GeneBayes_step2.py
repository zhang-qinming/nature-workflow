import pandas as pd
import sys


for i in [86,88,106] :
    df1 = pd.read_csv(f"pred_beta_Backman_2021_{i}.tsv", sep='\t')
    df2 = pd.read_csv(f"pred_magnitude_Backman_2021_{i}.tsv", sep='\t')
    df2["pred_magnitude"] = df2["pred_beta"]
    df2 = df2.drop(columns=["pred_beta"])
    df = df1.merge(df2, on="ensg")
    df.to_csv(f"features.Backman_2021_{i}.tsv", sep='\t', index=None)
