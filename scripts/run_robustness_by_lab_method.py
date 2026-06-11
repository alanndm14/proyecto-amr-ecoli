#!/usr/bin/env python3
"""Evalúa la robustez de la regresión logística por método de laboratorio."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
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
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RANDOM_STATE = 42
REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data" / "processed")
    parser.add_argument(
        "--outdir", type=Path, default=REPO_ROOT / "results" / "robustness_method_reproduced"
    )
    return parser.parse_args()


def make_model():
    return make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(
            max_iter=1000,
            class_weight="balanced",
            random_state=RANDOM_STATE,
            solver="liblinear",
            penalty="l1",
            C=0.1,
        ),
    )


def normalize_method(value) -> str:
    return "No especificado" if pd.isna(value) or not str(value).strip() else str(value).strip()


def metrics(y_true, y_pred, y_score=None) -> dict:
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "n": int(len(y_true)),
        "n_susceptible": int((y_true == 0).sum()),
        "n_resistant": int((y_true == 1).sum()),
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_score) if y_score is not None else np.nan,
        "recall_resistant": recall_score(y_true, y_pred, zero_division=0),
        "specificity_susceptible": tn / (tn + fp),
        "precision_resistant": precision_score(y_true, y_pred, zero_division=0),
        "f1_resistant": f1_score(y_true, y_pred, zero_division=0),
        "false_susceptible": int(fn),
        "false_resistant": int(fp),
        "true_susceptible": int(tn),
        "true_resistant": int(tp),
    }


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    datasets = {
        "ciprofloxacin": args.data_dir / "dataset_ciprofloxacin_model_ready.csv",
        "cefotaxime": args.data_dir / "dataset_cefotaxime_model_ready.csv",
    }
    holdout_records: list[dict] = []
    leave_out_records: list[dict] = []

    for antibiotic, path in datasets.items():
        df = pd.read_csv(path, low_memory=False)
        df["method_clean"] = df["laboratory_typing_method"].apply(normalize_method)
        features = [c for c in df.columns if c.startswith("amr__")]
        X = df[features].astype(np.int8)
        y = df["y"].astype(int)
        idx = np.arange(len(df))

        train_idx, test_idx = train_test_split(
            idx, test_size=0.20, stratify=y, random_state=RANDOM_STATE
        )
        model = make_model()
        model.fit(X.iloc[train_idx], y.iloc[train_idx])
        test = df.iloc[test_idx].copy()
        test["y_pred"] = model.predict(X.iloc[test_idx])
        test["y_score"] = model.predict_proba(X.iloc[test_idx])[:, 1]

        overall = metrics(test["y"], test["y_pred"], test["y_score"])
        overall.update(
            antibiotic=antibiotic,
            evaluation="holdout_subgroup",
            laboratory_typing_method="TOTAL_HOLDOUT",
        )
        holdout_records.append(overall)
        for method, group in test.groupby("method_clean"):
            if len(group) >= 100 and group["y"].nunique() == 2:
                row = metrics(group["y"], group["y_pred"], group["y_score"])
                row.update(
                    antibiotic=antibiotic,
                    evaluation="holdout_subgroup",
                    laboratory_typing_method=method,
                )
                holdout_records.append(row)

        counts = df.groupby("method_clean")["y"].agg(["size", "sum"])
        for method, values in counts.iterrows():
            n = int(values["size"])
            n_resistant = int(values["sum"])
            if method == "No especificado" or n < 500 or n_resistant < 20 or n - n_resistant < 20:
                continue
            test_mask = df["method_clean"].eq(method)
            train_mask = ~test_mask
            model = make_model()
            model.fit(X.loc[train_mask], y.loc[train_mask])
            score = model.predict_proba(X.loc[test_mask])[:, 1]
            row = metrics(y.loc[test_mask], (score >= 0.5).astype(int), score)
            row.update(
                antibiotic=antibiotic,
                evaluation="leave_method_out",
                laboratory_typing_method=method,
            )
            leave_out_records.append(row)

    columns = [
        "evaluation", "antibiotic", "laboratory_typing_method", "n", "n_susceptible",
        "n_resistant", "accuracy", "balanced_accuracy", "roc_auc", "recall_resistant",
        "specificity_susceptible", "precision_resistant", "f1_resistant",
        "false_susceptible", "false_resistant", "true_susceptible", "true_resistant",
    ]
    holdout = pd.DataFrame(holdout_records)[columns]
    leave_out = pd.DataFrame(leave_out_records)[columns]
    combined = pd.concat([holdout, leave_out], ignore_index=True)
    combined.to_csv(args.outdir / "robustez_metodo_laboratorio_todos.csv", index=False)
    holdout.to_csv(args.outdir / "robustez_holdout_por_metodo_laboratorio.csv", index=False)
    leave_out.to_csv(args.outdir / "robustez_leave_method_out.csv", index=False)

    plot = leave_out.sort_values(["antibiotic", "balanced_accuracy"])
    labels = plot["antibiotic"] + " - " + plot["laboratory_typing_method"]
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.45 * len(plot) + 1.5)))
    ax.barh(labels, plot["balanced_accuracy"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Balanced accuracy")
    ax.set_title("Leave-method-out: regresión logística balanceada")
    fig.tight_layout()
    fig.savefig(args.outdir / "figura_robustez_leave_method_out.png", dpi=250)
    plt.close(fig)
    print(f"Resultados de robustez: {args.outdir}")


if __name__ == "__main__":
    main()
