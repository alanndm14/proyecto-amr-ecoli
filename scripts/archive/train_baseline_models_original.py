#!/usr/bin/env python3
"""
Entrena modelos base para resistencia a antibióticos con genes AMR 0/1.

Uso:
python train_baseline_models.py --input dataset_ciprofloxacin_model_ready.csv --antibiotic ciprofloxacin
python train_baseline_models.py --input dataset_cefotaxime_model_ready.csv --antibiotic cefotaxime
"""

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
    classification_report,
    confusion_matrix,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="CSV model_ready")
    parser.add_argument("--antibiotic", required=True, help="Nombre del antibiótico")
    parser.add_argument("--outdir", default="model_outputs")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    df = pd.read_csv(args.input)
    feature_cols = [c for c in df.columns if c.startswith("amr__")]

    if "y" not in df.columns:
        raise ValueError("No encontré la columna y.")
    if not feature_cols:
        raise ValueError("No encontré columnas de genes AMR que empiecen con amr__.")

    X = df[feature_cols]
    y = df["y"].astype(int)

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.20,
        random_state=42,
        stratify=y
    )

    models = {
        "logistic_regression_balanced": make_pipeline(
            StandardScaler(with_mean=False),
            LogisticRegression(
                max_iter=5000,
                class_weight="balanced",
                solver="liblinear",
                random_state=42
            )
        ),
        "random_forest_balanced": RandomForestClassifier(
            n_estimators=500,
            max_depth=None,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1
        ),
    }

    rows = []
    for name, model in models.items():
        model.fit(X_train, y_train)
        pred = model.predict(X_test)

        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(X_test)[:, 1]
        else:
            proba = None

        metrics = {
            "antibiotic": args.antibiotic,
            "model": name,
            "n_samples": len(df),
            "n_features": len(feature_cols),
            "accuracy": accuracy_score(y_test, pred),
            "balanced_accuracy": balanced_accuracy_score(y_test, pred),
            "f1_resistant": f1_score(y_test, pred, pos_label=1),
            "roc_auc": roc_auc_score(y_test, proba) if proba is not None else None,
            "confusion_matrix": confusion_matrix(y_test, pred).tolist(),
            "classification_report": classification_report(y_test, pred, target_names=["Susceptible", "Resistant"], output_dict=True),
        }
        rows.append(metrics)

    simple_rows = []
    for m in rows:
        simple_rows.append({
            "antibiotic": m["antibiotic"],
            "model": m["model"],
            "n_samples": m["n_samples"],
            "n_features": m["n_features"],
            "accuracy": round(m["accuracy"], 4),
            "balanced_accuracy": round(m["balanced_accuracy"], 4),
            "f1_resistant": round(m["f1_resistant"], 4),
            "roc_auc": round(m["roc_auc"], 4) if m["roc_auc"] is not None else None,
            "confusion_matrix": json.dumps(m["confusion_matrix"])
        })

    pd.DataFrame(simple_rows).to_csv(outdir / f"metrics_{args.antibiotic}.csv", index=False)
    with open(outdir / f"metrics_{args.antibiotic}_full.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)

    print(pd.DataFrame(simple_rows).to_string(index=False))


if __name__ == "__main__":
    main()
