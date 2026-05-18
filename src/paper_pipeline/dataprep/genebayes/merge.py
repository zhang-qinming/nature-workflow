from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def merge_feature_predictions(
    trait_labels: list[str],
    prediction_dir: Path,
    feature_dir: Path,
) -> None:
    feature_dir.mkdir(parents=True, exist_ok=True)
    for trait_label in trait_labels:
        signed_path = prediction_dir / f"pred_beta_{trait_label}.tsv"
        magnitude_path = prediction_dir / f"pred_magnitude_{trait_label}.tsv"
        if not signed_path.exists():
            raise FileNotFoundError(signed_path)
        if not magnitude_path.exists():
            raise FileNotFoundError(magnitude_path)

        signed = pd.read_csv(signed_path, sep="\t")
        magnitude = pd.read_csv(magnitude_path, sep="\t").rename(
            columns={"pred_beta": "pred_magnitude"}
        )
        merged = signed.merge(magnitude, on="ensg")
        merged.to_csv(
            feature_dir / f"features.{trait_label}.tsv",
            sep="\t",
            index=False,
        )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Merge GeneBayes feature predictions.")
    parser.add_argument(
        "--traits",
        nargs="+",
        required=True,
        help="Trait labels such as Backman_2021_86.",
    )
    parser.add_argument("--prediction-dir", default=".")
    parser.add_argument("--feature-dir", default=".")
    args = parser.parse_args(argv)

    merge_feature_predictions(
        trait_labels=[str(value) for value in args.traits],
        prediction_dir=Path(args.prediction_dir),
        feature_dir=Path(args.feature_dir),
    )
