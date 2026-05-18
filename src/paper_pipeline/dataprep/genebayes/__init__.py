from .feature_prediction import (
    FeaturePredictionArgs,
    main_magnitude,
    main_signed,
    predict_magnitude_features,
    predict_signed_features,
)
from .merge import main as main_merge, merge_feature_predictions


def main_training(*args, **kwargs):
    from .training import main

    return main(*args, **kwargs)


def train_gene_bayes_model(*args, **kwargs):
    from .training import train_gene_bayes_model as _train_gene_bayes_model

    return _train_gene_bayes_model(*args, **kwargs)

__all__ = [
    "FeaturePredictionArgs",
    "main_magnitude",
    "main_merge",
    "main_signed",
    "main_training",
    "merge_feature_predictions",
    "predict_magnitude_features",
    "predict_signed_features",
    "train_gene_bayes_model",
]
