from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import xgboost as xgb
from hyperopt import Trials
from shaphypetune import BoostRFE


_SPACE = {
    "eta": 0.1,
    "max_depth": 4,
    "min_child_weight": 3,
    "n_estimators": 1000,
    "reg_alpha": 64,
    "reg_lambda": 0,
    "subsample": 0.5,
}
_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


@dataclass(slots=True)
class FeaturePredictionArgs:
    trait_label: str
    s_het_path: Path
    embeddings_path: Path
    celltype_ntpm_path: Path
    geneformer_path: Path
    trait_path: Path
    output_dir: Path
    prefer_gpu: bool = True


def _build_parser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("trait_label")
    parser.add_argument("--s-het", dest="s_het_path", default="s_het.tsv")
    parser.add_argument("--embeddings", dest="embeddings_path", default="reduced_embeddings.tsv")
    parser.add_argument("--celltype-ntpm", dest="celltype_ntpm_path", default="celltype_nTPM.tsv")
    parser.add_argument(
        "--geneformer-mixture",
        dest="geneformer_path",
        default="GeneFormer_cellclassifier_mixture2.tsv",
    )
    parser.add_argument("--traits-dir", dest="traits_dir", default="../traits")
    parser.add_argument("--trait-file", dest="trait_file", default=None)
    parser.add_argument("--output-dir", dest="output_dir", default=".")
    parser.add_argument(
        "--prefer-gpu",
        dest="prefer_gpu",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use the GPU XGBoost backend when available in the environment.",
    )
    return parser


def _parse_args(argv: list[str] | None, description: str) -> FeaturePredictionArgs:
    parser = _build_parser(description)
    args = parser.parse_args(argv)
    trait_path = (
        Path(args.trait_file)
        if args.trait_file is not None
        else Path(args.traits_dir) / f"{args.trait_label}.summary_statistics.csv"
    )
    return FeaturePredictionArgs(
        trait_label=args.trait_label,
        s_het_path=Path(args.s_het_path),
        embeddings_path=Path(args.embeddings_path),
        celltype_ntpm_path=Path(args.celltype_ntpm_path),
        geneformer_path=Path(args.geneformer_path),
        trait_path=trait_path,
        output_dir=Path(args.output_dir),
        prefer_gpu=bool(args.prefer_gpu),
    )


def _load_training_data(args: FeaturePredictionArgs) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    s_het = pd.read_csv(args.s_het_path, sep="\t")[["hgnc", "ensg", "chrom", "post_mean"]]
    embed = pd.read_csv(args.embeddings_path, sep=" ")
    gex = pd.read_csv(args.celltype_ntpm_path, sep="\t")
    geneformer = pd.read_csv(args.geneformer_path, sep="\t")
    trait = pd.read_csv(args.trait_path)

    features = s_het.merge(embed, left_on="hgnc", right_on="gene", how="outer")
    features = features.merge(gex, on="ensg", how="outer")
    features = features.merge(geneformer, on="ensg", how="outer")
    data = features.merge(trait, on="ensg").drop_duplicates(subset=["ensg"])
    x_frame = data.drop(
        [
            col
            for col in ["beta", "chrom", "ensg", "hgnc", "gene", "gene_x", "gene_y", "standard_error", "z"]
            if col in data.columns
        ],
        axis=1,
    )
    genes = data["ensg"]
    return data, x_frame, genes


def _spearman_corr(y_hat, dtrain):
    y_true = dtrain.get_label()
    return "corr", -scipy.stats.spearmanr(y_true, y_hat, nan_policy="raise")[0]


def _build_regressor(prefer_gpu: bool) -> xgb.XGBRegressor:
    if prefer_gpu:
        return xgb.XGBRegressor(
            gpu_id=0,
            tree_method="gpu_hist",
            max_depth=_SPACE["max_depth"],
            reg_alpha=_SPACE["reg_alpha"],
            reg_lambda=_SPACE["reg_lambda"],
            min_child_weight=_SPACE["min_child_weight"],
            eta=_SPACE["eta"],
            subsample=_SPACE["subsample"],
            n_estimators=_SPACE["n_estimators"],
        )
    return xgb.XGBRegressor(
        tree_method="hist",
        max_depth=_SPACE["max_depth"],
        reg_alpha=_SPACE["reg_alpha"],
        reg_lambda=_SPACE["reg_lambda"],
        min_child_weight=_SPACE["min_child_weight"],
        eta=_SPACE["eta"],
        subsample=_SPACE["subsample"],
        n_estimators=_SPACE["n_estimators"],
    )


def _fit_feature_model(
    *,
    x_train: pd.DataFrame,
    y_train: pd.Series,
    x_test: pd.DataFrame,
    y_test: pd.Series,
    w_train: pd.Series,
    w_test: pd.Series,
    prefer_gpu: bool,
):
    model = BoostRFE(
        _build_regressor(prefer_gpu),
        step=1,
        importance_type="shap_importances",
        n_jobs=-1,
        verbose=0,
        param_grid=_SPACE,
        n_iter=150,
        min_features_to_select=x_train.shape[1],
    )
    trials = Trials()
    model.fit(
        x_train,
        y_train,
        eval_set=[(x_test, y_test)],
        verbose=False,
        early_stopping_rounds=20,
        trials=trials,
        sample_weight=w_train,
        sample_weight_eval_set=[w_test],
        eval_metric=_spearman_corr,
    )
    return model.estimator_


def _predict_features(
    args: FeaturePredictionArgs,
    *,
    data: pd.DataFrame,
    x_frame: pd.DataFrame,
    genes: pd.Series,
    y: pd.Series,
    weights: pd.Series,
    prediction_name: str,
    stats_name: str,
) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("X.shape", x_frame.shape)

    s_het_corr: list[float] = []
    val_corr: list[float] = []
    train_corr: list[float] = []
    preds_all = np.zeros(len(x_frame), dtype=float)

    print(x_frame)
    print(y)

    for chrom in _CHROMS:
        print(chrom)

        mask = data["chrom"].isin([chrom])
        x_train = x_frame.loc[~mask]
        y_train = y.loc[~mask]
        x_test = x_frame.loc[mask]
        y_test = y.loc[mask]
        w_train = weights.loc[~mask]
        w_test = weights.loc[mask]

        corr = scipy.stats.spearmanr(x_test["post_mean"], y_test)
        print("correlation with s_het", corr)
        s_het_corr.append(corr[0])

        try:
            estimator = _fit_feature_model(
                x_train=x_train,
                y_train=y_train,
                x_test=x_test,
                y_test=y_test,
                w_train=w_train,
                w_test=w_test,
                prefer_gpu=args.prefer_gpu,
            )
        except xgb.core.XGBoostError:
            if not args.prefer_gpu:
                raise
            print("GPU XGBoost backend unavailable, retrying on CPU.")
            estimator = _fit_feature_model(
                x_train=x_train,
                y_train=y_train,
                x_test=x_test,
                y_test=y_test,
                w_train=w_train,
                w_test=w_test,
                prefer_gpu=False,
            )
        preds_val = estimator.predict(x_test)
        preds_train = estimator.predict(x_train)

        corr = scipy.stats.spearmanr(preds_val, y_test)
        val_corr.append(corr[0])
        print("val", corr)

        corr = scipy.stats.spearmanr(preds_train, y_train)
        train_corr.append(corr[0])
        print("train", corr)

        preds_all += estimator.predict(x_frame)
        print()

    print("s_het_corr")
    for value in s_het_corr:
        print(value)
    print("val_corr")
    for value in val_corr:
        print(value)
    print("train corr")
    for value in train_corr:
        print(value)

    print(scipy.stats.spearmanr(preds_all, y))

    pd.DataFrame({"ensg": genes, "pred_beta": preds_all}).to_csv(
        args.output_dir / prediction_name,
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        {"s_het_corr": s_het_corr, "val_corr": val_corr, "train_corr": train_corr}
    ).to_csv(
        args.output_dir / stats_name,
        sep="\t",
        index=False,
    )


def predict_signed_features(args: FeaturePredictionArgs) -> None:
    data, x_frame, genes = _load_training_data(args)
    y = data["beta"] * 10e2
    weights = 1.0 / (data["standard_error"] ** 2)
    _predict_features(
        args,
        data=data,
        x_frame=x_frame,
        genes=genes,
        y=y,
        weights=weights,
        prediction_name=f"pred_beta_{args.trait_label}.tsv",
        stats_name=f"train_stats.beta_{args.trait_label}.tsv",
    )


def predict_magnitude_features(args: FeaturePredictionArgs) -> None:
    data, x_frame, genes = _load_training_data(args)
    y = (data["beta"] ** 2) * 10e2 - (data["standard_error"] ** 2) * 10e2
    y = y.rank(method="dense")
    weights = 1.0 / (data["standard_error"] ** 2)
    weights[:] = 1
    _predict_features(
        args,
        data=data,
        x_frame=x_frame,
        genes=genes,
        y=y,
        weights=weights,
        prediction_name=f"pred_magnitude_{args.trait_label}.tsv",
        stats_name=f"train_stats.magnitude_{args.trait_label}.tsv",
    )


def main_signed(argv: list[str] | None = None) -> None:
    predict_signed_features(_parse_args(argv, "Prepare signed GeneBayes feature predictions."))


def main_magnitude(argv: list[str] | None = None) -> None:
    predict_magnitude_features(
        _parse_args(argv, "Prepare magnitude-style GeneBayes feature predictions.")
    )
