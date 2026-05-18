from shaphypetune import BoostSearch, BoostRFE, BoostRFA, BoostBoruta

import sys
import scipy
import pickle
import pandas as pd
import xgboost as xgb
from sklearn.metrics import mean_squared_error
import numpy as np
from sklearn.model_selection import train_test_split
from hyperopt import STATUS_OK, Trials, fmin, hp, tpe
import time
import gffutils
from statsmodels.stats.multitest import fdrcorrection

from sklearn.datasets import dump_svmlight_file

import seaborn as sns
import matplotlib.pyplot as plt

TRAIT = sys.argv[1]

### s_het ###
s_het = pd.read_csv("s_het.tsv", sep='\t')
s_het = s_het[["hgnc","ensg","chrom","post_mean"]]#,"post_lower_95","post_upper_95"]]

### protein embeddings ###
embed = pd.read_csv("reduced_embeddings.tsv", sep=' ')

gex = pd.read_csv("celltype_nTPM.tsv", sep='\t')


##GeneFormer##
GF=pd.read_csv("GeneFormer_cellclassifier_mixture2.tsv", sep='\t')


### MERGING ###
features = s_het
print("s_het", s_het.columns, s_het.shape)
features = features.merge(embed, left_on="hgnc", right_on="gene", how="outer")
features = features.merge(gex, on="ensg", how="outer")
features = features.merge(GF, on="ensg", how="outer")
trait = pd.read_csv(f"../traits/{TRAIT}.summary_statistics.csv")
data = features.merge(trait, on="ensg")

data = data.drop_duplicates(subset=["ensg"])

X = data.drop([col for col in ["beta","chrom","ensg","hgnc","gene","gene_x", "gene_y", "standard_error","z"] if col in data.columns], axis=1)
genes = data["ensg"]
# to limit numerical issues
y = (data["beta"]**2)*10e2-(data["standard_error"]**2)*10e2
y = y.rank(method='dense')
weights = 1./(data["standard_error"]**2)
weights[:] = 1
space = {'eta': 0.1, 'max_depth': 4, 'min_child_weight': 3, 'n_estimators': 1000, 'reg_alpha': 64, 'reg_lambda': 0, 'subsample': 0.5}


def spearman_corr(y_hat, dtrain):
    y_true = dtrain.get_label()
    return "corr", -scipy.stats.spearmanr(y_true, y_hat, nan_policy="raise")[0]

chroms = [f"chr{i}" for i in range(1,23)]+["chrX"]

print("X.shape", X.shape)

s_het_corr = []
val_corr = []
train_corr = []
genes_all, preds_all = [], 0
truth_all = []

print(X)
print(y)

for chrom in chroms:
    print(chrom)

    X_train1 = X.loc[~data["chrom"].isin([chrom])]
    y_train1 = y.loc[~data["chrom"].isin([chrom])]
    X_test1 = X.loc[data["chrom"].isin([chrom])]
    y_test1 = y.loc[data["chrom"].isin([chrom])]

    w_train1 = weights.loc[~data["chrom"].isin([chrom])]
    w_test1 = weights.loc[data["chrom"].isin([chrom])]

    corr = scipy.stats.spearmanr(X_test1["post_mean"],y_test1)
    print("correlation with s_het", corr)
    s_het_corr.append(corr[0])

    clf = xgb.XGBRegressor(
            gpu_id=0,
            tree_method="gpu_hist",
            max_depth= space['max_depth'],
            reg_alpha = space['reg_alpha'],
            reg_lambda = space['reg_lambda'],
            min_child_weight = space['min_child_weight'],
            eta = space['eta'],
            subsample = space['subsample'],
            n_estimators = space['n_estimators'])
    
    model = BoostRFE(
        clf, step=1, importance_type='shap_importances', n_jobs=-1, verbose=0, param_grid=space, n_iter=150, min_features_to_select=X.shape[1]
    )

    trials = Trials()
    model.fit(X_train1, y_train1, eval_set=[(X_test1, y_test1)],
              verbose=False, early_stopping_rounds=20, trials=trials,
              sample_weight=w_train1, sample_weight_eval_set=[w_test1],
              eval_metric=spearman_corr)
    clf = model.estimator_
    preds_val = clf.predict(X_test1)
    preds_train = clf.predict(X_train1)
    corr = scipy.stats.spearmanr(preds_val, y_test1)
    val_corr.append(corr[0])
    print("val",corr)
    corr = scipy.stats.spearmanr(preds_train, y_train1)
    train_corr.append(corr[0])
    print("train",corr)

    genes_all += genes.loc[data["chrom"].isin([chrom])].tolist()
    preds_all += clf.predict(X) #preds_val.tolist()
    truth_all += y_test1.tolist()
    print()

truth_all = y
genes_all = genes

print("s_het_corr")
for x in s_het_corr:
    print(x)
print("val_corr")
for x in val_corr:
    print(x)
print("train corr")
for x in train_corr:
    print(x)

print(scipy.stats.spearmanr(preds_all, truth_all))

df_out = pd.DataFrame({"ensg":genes_all,"pred_beta":preds_all})
df_out.to_csv(f"pred_magnitude_{TRAIT}.tsv", sep='\t', index=False)

df_out = pd.DataFrame({"s_het_corr":s_het_corr,"val_corr":val_corr,"train_corr":train_corr})
df_out.to_csv(f"train_stats.magnitude_{TRAIT}.tsv", sep='\t', index=False)
