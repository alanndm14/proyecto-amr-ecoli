#!/usr/bin/env python3
"""Evaluate a weighted multi-seed MLP as the fourth AMR model."""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.utils.class_weight import compute_sample_weight

RANDOM_STATE = 42
REPO_ROOT = Path(__file__).resolve().parents[1]
SEEDS = [1, 2, 3, 4, 5, 42, 99]
CONFIGS = {"MLP-16": (16,), "MLP-32": (32,), "MLP-16-8": (16, 8)}
DATASETS = {
    "ciprofloxacin": "dataset_ciprofloxacin_model_ready.csv",
    "cefotaxime": "dataset_cefotaxime_model_ready.csv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data" / "processed")
    parser.add_argument("--outdir", type=Path, default=REPO_ROOT / "results" / "mlp_reproduced")
    parser.add_argument(
        "--reference-metrics",
        type=Path,
        default=REPO_ROOT / "results" / "optimized" / "metricas_holdout_todos_los_modelos.csv",
    )
    return parser.parse_args()


def build_mlp(hidden_layers: tuple[int, ...], seed: int) -> MLPClassifier:
    return MLPClassifier(
        hidden_layer_sizes=hidden_layers,
        activation="relu",
        solver="adam",
        alpha=1e-4,
        learning_rate_init=1e-3,
        batch_size=512,
        early_stopping=True,
        validation_fraction=0.15,
        n_iter_no_change=5,
        tol=1e-3,
        max_iter=70,
        random_state=seed,
    )


def metrics(y_true, y_pred, y_score) -> dict:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_score),
        "f1_resistant": f1_score(y_true, y_pred, zero_division=0),
        "recall_resistant": recall_score(y_true, y_pred, zero_division=0),
        "precision_resistant": precision_score(y_true, y_pred, zero_division=0),
        "specificity_susceptible": tn / (tn + fp),
        "false_susceptible": int(fn),
        "false_resistant": int(fp),
        "true_resistant": int(tp),
        "true_susceptible": int(tn),
        "confusion_matrix": [[int(tn), int(fp)], [int(fn), int(tp)]],
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    ensemble_rows: list[dict] = []

    for antibiotic, filename in DATASETS.items():
        df = pd.read_csv(args.data_dir / filename, low_memory=False)
        feature_cols = [c for c in df.columns if c.startswith("amr__")]
        X = df[feature_cols].astype(np.float32).values
        y = df["y"].astype(int).values
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.20, stratify=y, random_state=RANDOM_STATE
        )
        weights = compute_sample_weight("balanced", y_train)

        for candidate_name, hidden in CONFIGS.items():
            probabilities = []
            for seed in SEEDS:
                model = build_mlp(hidden, seed)
                model.fit(X_train, y_train, sample_weight=weights)
                score = model.predict_proba(X_test)[:, 1]
                row = metrics(y_test, (score >= 0.5).astype(int), score)
                row.update(
                    antibiotic=antibiotic,
                    model="mlp_weighted",
                    model_label="MLP ponderado",
                    candidate_name=candidate_name,
                    hidden_layer_sizes=str(hidden),
                    seed=seed,
                    n_iter=int(model.n_iter_),
                )
                rows.append(row)
                probabilities.append(score)

            score = np.mean(probabilities, axis=0)
            row = metrics(y_test, (score >= 0.5).astype(int), score)
            row.update(
                antibiotic=antibiotic,
                model="mlp_weighted_ensemble",
                model_label="MLP ponderado (promedio 7 semillas)",
                candidate_name=candidate_name,
                hidden_layer_sizes=str(hidden),
                seed="ensemble_7",
                n_seeds=len(SEEDS),
            )
            ensemble_rows.append(row)

    runs = pd.DataFrame(rows)
    ensembles = pd.DataFrame(ensemble_rows)
    summary = (
        runs.groupby(["antibiotic", "candidate_name", "hidden_layer_sizes"])
        .agg(
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            balanced_accuracy_std=("balanced_accuracy", "std"),
            roc_auc_mean=("roc_auc", "mean"),
            roc_auc_std=("roc_auc", "std"),
            f1_resistant_mean=("f1_resistant", "mean"),
            f1_resistant_std=("f1_resistant", "std"),
            recall_resistant_mean=("recall_resistant", "mean"),
            recall_resistant_std=("recall_resistant", "std"),
            false_susceptible_mean=("false_susceptible", "mean"),
            false_resistant_mean=("false_resistant", "mean"),
            n_seeds=("seed", "nunique"),
        )
        .reset_index()
        .sort_values(["antibiotic", "balanced_accuracy_mean"], ascending=[True, False])
    )

    runs.to_csv(args.outdir / "mlp_holdout_por_semilla.csv", index=False)
    summary.to_csv(args.outdir / "mlp_holdout_resumen_semillas.csv", index=False)
    ensembles.to_csv(args.outdir / "mlp_holdout_ensemble_semillas.csv", index=False)

    reference = pd.read_csv(args.reference_metrics)
    best = ensembles.sort_values(
        ["antibiotic", "balanced_accuracy"], ascending=[True, False]
    ).groupby("antibiotic").head(1)
    combined = pd.concat(
        [reference, best.assign(selected_for_final_report=False)], ignore_index=True, sort=False
    )
    combined.to_csv(args.outdir / "comparacion_holdout_modelos_previos_mas_mlp.csv", index=False)
    shutil.make_archive(str(args.outdir), "zip", args.outdir)
    print(f"Resultados MLP: {args.outdir}")


if __name__ == "__main__":
    main()
