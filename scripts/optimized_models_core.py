#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Análisis final optimizado de resistencia antimicrobiana en Escherichia coli.

Modelos comparados:
1. Regresión logística balanceada
2. Random forest balanceado
3. XGBoost ponderado

Qué hace:
- Lee dataset_ciprofloxacin_model_ready.csv y dataset_cefotaxime_model_ready.csv
- Separa holdout estratificado del 20%
- Ajusta hiperparámetros con RandomizedSearchCV
- Usa RepeatedStratifiedKFold con 5 folds x 2 repeticiones
- Selecciona por balanced accuracy
- Evalúa en holdout
- Calcula intervalos de confianza bootstrap
- Genera figuras en español
- Evalúa estabilidad de genes importantes
- Genera reporte final en Markdown y ZIP de resultados

Uso en terminal:
python generate_final_optimized_amr_analysis.py

Uso en Jupyter:
!python generate_final_optimized_amr_analysis.py

Dependencias:
pip install pandas numpy matplotlib scikit-learn xgboost tabulate joblib
"""

from __future__ import annotations

import argparse
import json
import re
import time
import warnings
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import RandomizedSearchCV, RepeatedStratifiedKFold, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    from xgboost import XGBClassifier
except ImportError as exc:
    raise ImportError(
        "No encontré xgboost. Instálalo con:\n\n"
        "pip install xgboost\n\n"
        "o desde Jupyter:\n\n"
        "%pip install xgboost\n"
    ) from exc

try:
    import joblib
except ImportError:
    joblib = None


# -----------------------------------------------------------------------------
# Configuración general
# -----------------------------------------------------------------------------

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 120

RANDOM_STATE = 42

DATASETS = {
    "ciprofloxacin": Path("dataset_ciprofloxacin_model_ready.csv"),
    "cefotaxime": Path("dataset_cefotaxime_model_ready.csv"),
}

SPANISH_MODEL_NAMES = {
    "logistic_regression_balanced": "Regresión logística balanceada",
    "random_forest_balanced": "Bosque aleatorio balanceado",
    "xgboost_weighted": "XGBoost ponderado",
}

METRIC_LABELS_ES = {
    "accuracy": "Exactitud",
    "balanced_accuracy": "Exactitud balanceada",
    "f1_resistant": "F1 de resistentes",
    "roc_auc": "ROC-AUC",
    "recall_resistant": "Sensibilidad de resistentes",
    "precision_resistant": "Precisión de resistentes",
    "specificity_susceptible": "Especificidad de susceptibles",
}


# -----------------------------------------------------------------------------
# Utilidades de formato y métricas
# -----------------------------------------------------------------------------


def fmt_sig(x: Any, sig: int = 4) -> str:
    """Formatea con cifras significativas sin redondeo excesivo."""
    if x is None:
        return "NA"
    try:
        if np.isnan(x):
            return "NA"
    except TypeError:
        pass
    if x == 0:
        return "0"
    return f"{x:.{sig}g}"


def fmt_pct(x: float, sig: int = 3) -> str:
    return f"{100 * x:.{sig}g}%"


def clean_gene_name(col: str) -> str:
    return re.sub(r"^amr__", "", col)


def specificity_score(y_true: Iterable[int], y_pred: Iterable[int]) -> float:
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    return tn / (tn + fp) if (tn + fp) else np.nan


def bootstrap_metric_ci(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_proba: np.ndarray,
    metric_name: str,
    n_boot: int = 500,
    seed: int = RANDOM_STATE,
) -> Tuple[float, float]:
    """Intervalo de confianza percentil bootstrap para métricas en holdout."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    values: List[float] = []

    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yt = y_true[idx]
        yp = y_pred[idx]
        pr = y_proba[idx]

        if metric_name == "accuracy":
            val = accuracy_score(yt, yp)
        elif metric_name == "balanced_accuracy":
            val = balanced_accuracy_score(yt, yp)
        elif metric_name == "f1_resistant":
            val = f1_score(yt, yp, pos_label=1, zero_division=0)
        elif metric_name == "recall_resistant":
            val = recall_score(yt, yp, pos_label=1, zero_division=0)
        elif metric_name == "precision_resistant":
            val = precision_score(yt, yp, pos_label=1, zero_division=0)
        elif metric_name == "specificity_susceptible":
            val = specificity_score(yt, yp)
        elif metric_name == "roc_auc":
            if len(np.unique(yt)) < 2:
                continue
            val = roc_auc_score(yt, pr)
        else:
            raise ValueError(f"Métrica no reconocida: {metric_name}")

        values.append(val)

    arr = np.asarray(values, dtype=float)
    if len(arr) == 0:
        return np.nan, np.nan
    lo, hi = np.nanpercentile(arr, [2.5, 97.5])
    return float(lo), float(hi)


def evaluate_predictions(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_proba: np.ndarray,
    n_boot: int,
) -> Tuple[Dict[str, Any], np.ndarray]:
    metrics: Dict[str, Any] = {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "f1_resistant": f1_score(y_true, y_pred, pos_label=1, zero_division=0),
        "roc_auc": roc_auc_score(y_true, y_proba),
        "recall_resistant": recall_score(y_true, y_pred, pos_label=1, zero_division=0),
        "precision_resistant": precision_score(y_true, y_pred, pos_label=1, zero_division=0),
        "specificity_susceptible": specificity_score(y_true, y_pred),
    }

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()

    metrics.update(
        {
            "true_susceptible": int(tn),
            "false_resistant": int(fp),
            "false_susceptible": int(fn),
            "true_resistant": int(tp),
            "confusion_matrix": json.dumps(cm.tolist()),
        }
    )

    for metric_name in [
        "accuracy",
        "balanced_accuracy",
        "f1_resistant",
        "roc_auc",
        "recall_resistant",
        "precision_resistant",
        "specificity_susceptible",
    ]:
        lo, hi = bootstrap_metric_ci(
            y_true,
            y_pred,
            y_proba,
            metric_name=metric_name,
            n_boot=n_boot,
            seed=RANDOM_STATE,
        )
        metrics[f"{metric_name}_ci95_low"] = lo
        metrics[f"{metric_name}_ci95_high"] = hi

    return metrics, cm


def round_for_export(df: pd.DataFrame, ndigits: int = 6) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].round(ndigits)
    return out


def md_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except Exception:
        return df.to_string(index=False)


# -----------------------------------------------------------------------------
# Figuras en español
# -----------------------------------------------------------------------------


def plot_class_balance(antibiotic: str, y: pd.Series, outdir: Path) -> None:
    counts = (
        pd.Series(y)
        .map({0: "Susceptible", 1: "Resistente"})
        .value_counts()
        .reindex(["Susceptible", "Resistente"])
    )

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    bars = ax.bar(counts.index, counts.values)

    ax.set_title(f"Distribución de clases para {antibiotic}")
    ax.set_xlabel("Fenotipo")
    ax.set_ylabel("Número de genomas")

    total = counts.sum()
    for bar, val in zip(bars, counts.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val,
            f"{val}\n({fmt_pct(val / total)})",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(outdir / f"figura_01_balance_clases_{antibiotic}.png", dpi=260)
    plt.close(fig)


def plot_model_comparison(antibiotic: str, comp_df: pd.DataFrame, outdir: Path) -> None:
    data = comp_df.copy()
    data = data.sort_values("cv_balanced_accuracy_mean", ascending=False)

    labels = data["model_label"].tolist()
    means = data["cv_balanced_accuracy_mean"].to_numpy()
    sds = data["cv_balanced_accuracy_std"].to_numpy()

    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    bars = ax.bar(labels, means, yerr=sds, capsize=5)

    ax.set_ylim(0, 1)
    ax.set_ylabel("Exactitud balanceada en validación cruzada")
    ax.set_xlabel("Modelo")
    ax.set_title(f"Comparación de modelos para {antibiotic}")
    ax.tick_params(axis="x", rotation=15)

    for bar, mean, sd in zip(bars, means, sds):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            mean,
            f"{fmt_sig(mean)} ± {fmt_sig(sd, 2)}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(outdir / f"figura_02_comparacion_modelos_{antibiotic}.png", dpi=260)
    plt.close(fig)


def plot_metric_panel(antibiotic: str, holdout_df: pd.DataFrame, outdir: Path) -> None:
    data = holdout_df.copy()
    metrics = ["balanced_accuracy", "recall_resistant", "f1_resistant", "roc_auc"]
    metric_labels = [METRIC_LABELS_ES[m] for m in metrics]

    x = np.arange(len(data))
    width = 0.18

    fig, ax = plt.subplots(figsize=(10.5, 5.8))

    for i, metric in enumerate(metrics):
        vals = data[metric].to_numpy()
        ax.bar(x + (i - 1.5) * width, vals, width, label=metric_labels[i])

    ax.set_xticks(x)
    ax.set_xticklabels(data["model_label"], rotation=15, ha="right")
    ax.set_ylim(0, 1)
    ax.set_xlabel("Modelo")
    ax.set_ylabel("Valor de la métrica")
    ax.set_title(f"Métricas en conjunto retenido para {antibiotic}")
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / f"figura_03_metricas_holdout_{antibiotic}.png", dpi=260)
    plt.close(fig)


def plot_confusion_matrix_es(antibiotic: str, model_label: str, cm: np.ndarray, outdir: Path) -> None:
    row_sums = cm.sum(axis=1, keepdims=True)
    norm = np.divide(cm, row_sums, out=np.zeros_like(cm, dtype=float), where=row_sums != 0)

    fig, ax = plt.subplots(figsize=(6.0, 5.4))
    im = ax.imshow(norm, interpolation="nearest")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Proporción por clase real")

    classes = ["Susceptible", "Resistente"]
    ax.set(
        xticks=np.arange(2),
        yticks=np.arange(2),
        xticklabels=classes,
        yticklabels=classes,
        xlabel="Clase predicha",
        ylabel="Clase real",
        title=f"Matriz de confusión\n{antibiotic} — {model_label}",
    )

    thresh = norm.max() / 2 if norm.size else 0.5
    for i in range(2):
        for j in range(2):
            text = f"{cm[i, j]}\n{fmt_pct(norm[i, j])}"
            ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                color="white" if norm[i, j] > thresh else "black",
                fontsize=10,
            )

    fig.tight_layout()
    fig.savefig(outdir / f"figura_04_matriz_confusion_{antibiotic}.png", dpi=260)
    plt.close(fig)


def plot_roc_es(
    antibiotic: str,
    model_label: str,
    y_true: np.ndarray,
    y_proba: np.ndarray,
    auc_value: float,
    outdir: Path,
) -> None:
    fpr, tpr, _ = roc_curve(y_true, y_proba)

    fig, ax = plt.subplots(figsize=(6.0, 5.4))
    ax.plot(fpr, tpr, label=f"Curva ROC (AUC = {fmt_sig(auc_value)})")
    ax.plot([0, 1], [0, 1], linestyle="--", label="Azar")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Tasa de falsos positivos")
    ax.set_ylabel("Tasa de verdaderos positivos")
    ax.set_title(f"Curva ROC\n{antibiotic} — {model_label}")
    ax.legend(loc="lower right")

    fig.tight_layout()
    fig.savefig(outdir / f"figura_05_curva_roc_{antibiotic}.png", dpi=260)
    plt.close(fig)


def plot_stable_genes(
    antibiotic: str,
    model_label: str,
    stability: pd.DataFrame,
    total_folds: int,
    outdir: Path,
    n: int = 15,
) -> None:
    top = stability.head(n).iloc[::-1].copy()

    fig, ax = plt.subplots(figsize=(9.0, 6.8))
    bars = ax.barh(top["gene"], top["selection_frequency"])

    ax.set_xlabel("Número de veces que el gen apareció entre los más importantes")
    ax.set_ylabel("Gen AMR")
    ax.set_title(f"Genes más estables\n{antibiotic} — {model_label}")
    ax.set_xlim(0, total_folds * 1.08)

    for bar, freq, frac in zip(top.index.map(lambda i: bars[list(top.index).index(i)]), top["selection_frequency"], top["selection_fraction"]):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2,
            f" {int(freq)}/{total_folds} ({fmt_pct(frac)})",
            va="center",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(outdir / f"figura_06_genes_estables_{antibiotic}.png", dpi=260)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Modelos e hiperparámetros solicitados
# -----------------------------------------------------------------------------


def make_candidate_models(scale_pos_weight: float) -> Dict[str, Dict[str, Any]]:
    """Modelos con los hiperparámetros ampliados solicitados."""
    return {
        "logistic_regression_balanced": {
            "label": SPANISH_MODEL_NAMES["logistic_regression_balanced"],
            "n_iter": 40,
            "estimator": make_pipeline(
                StandardScaler(with_mean=False),
                LogisticRegression(
                    max_iter=8000,
                    class_weight="balanced",
                    random_state=RANDOM_STATE,
                ),
            ),
            "param_dist": [
                {
                    "logisticregression__solver": ["liblinear"],
                    "logisticregression__penalty": ["l1", "l2"],
                    "logisticregression__C": [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100],
                },
                {
                    "logisticregression__solver": ["saga"],
                    "logisticregression__penalty": ["l1", "l2"],
                    "logisticregression__C": [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100],
                },
                {
                    "logisticregression__solver": ["saga"],
                    "logisticregression__penalty": ["elasticnet"],
                    "logisticregression__C": [0.01, 0.03, 0.1, 0.3, 1, 3, 10],
                    "logisticregression__l1_ratio": [0.15, 0.3, 0.5, 0.7, 0.85],
                },
            ],
        },
        "random_forest_balanced": {
            "label": SPANISH_MODEL_NAMES["random_forest_balanced"],
            "n_iter": 50,
            "estimator": RandomForestClassifier(
                random_state=RANDOM_STATE,
                n_jobs=1,
            ),
            "param_dist": {
                "n_estimators": [300, 500, 800],
                "class_weight": ["balanced", "balanced_subsample"],
                "max_depth": [None, 10, 15, 20, 25, 35, 50],
                "min_samples_leaf": [1, 2, 3, 4, 6, 8],
                "min_samples_split": [2, 4, 6, 10],
                "max_features": ["sqrt", "log2", 0.25, 0.5],
                "bootstrap": [True],
                "max_samples": [None, 0.7, 0.85],
            },
        },
        "xgboost_weighted": {
            "label": SPANISH_MODEL_NAMES["xgboost_weighted"],
            "n_iter": 60,
            "estimator": XGBClassifier(
                objective="binary:logistic",
                eval_metric="logloss",
                tree_method="hist",
                scale_pos_weight=scale_pos_weight,
                random_state=RANDOM_STATE,
                n_jobs=1,
            ),
            "param_dist": {
                "n_estimators": [300, 500, 800, 1200],
                "learning_rate": [0.01, 0.03, 0.05, 0.08, 0.1],
                "max_depth": [2, 3, 4, 5, 6],
                "min_child_weight": [1, 2, 5, 10],
                "subsample": [0.7, 0.85, 1.0],
                "colsample_bytree": [0.7, 0.85, 1.0],
                "gamma": [0, 0.1, 0.5, 1],
                "reg_alpha": [0, 0.01, 0.1, 1],
                "reg_lambda": [0.5, 1, 2, 5, 10],
            },
        },
    }


# -----------------------------------------------------------------------------
# Importancia y estabilidad de genes
# -----------------------------------------------------------------------------


def extract_importance(model: Any, feature_cols: List[str], model_family: str) -> pd.DataFrame:
    if model_family == "logistic_regression_balanced":
        vals = model.named_steps["logisticregression"].coef_[0]
        abs_vals = np.abs(vals)
        importance_type = "coeficiente_logistico_absoluto"
    elif model_family == "random_forest_balanced":
        vals = model.feature_importances_
        abs_vals = vals
        importance_type = "importancia_gini_bosque_aleatorio"
    elif model_family == "xgboost_weighted":
        vals = model.feature_importances_
        abs_vals = vals
        importance_type = "importancia_xgboost_normalizada"
    else:
        raise ValueError(f"Modelo no reconocido para importancia: {model_family}")

    return (
        pd.DataFrame(
            {
                "feature": feature_cols,
                "gene": [clean_gene_name(c) for c in feature_cols],
                "importance": vals,
                "abs_importance": abs_vals,
                "importance_type": importance_type,
            }
        )
        .sort_values("abs_importance", ascending=False)
        .reset_index(drop=True)
    )


def summarize_gene_stability(
    estimator: Any,
    model_family: str,
    X: pd.DataFrame,
    y: pd.Series,
    feature_cols: List[str],
    antibiotic: str,
    n_splits: int,
    n_repeats: int,
    top_n: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, int]:
    """Reentrena el modelo seleccionado en folds repetidos y mide estabilidad del top genes."""
    cv = RepeatedStratifiedKFold(
        n_splits=n_splits,
        n_repeats=n_repeats,
        random_state=RANDOM_STATE,
    )

    records = []
    total_folds = n_splits * n_repeats

    for fold, (tr, te) in enumerate(cv.split(X, y), start=1):
        model = clone(estimator)
        model.fit(X.iloc[tr], y.iloc[tr])

        top = extract_importance(model, feature_cols, model_family).head(top_n).copy()
        top["fold"] = fold
        top["rank"] = np.arange(1, len(top) + 1)
        top["antibiotic"] = antibiotic
        top["model"] = model_family
        top["model_label"] = SPANISH_MODEL_NAMES[model_family]
        records.append(top)

    raw = pd.concat(records, ignore_index=True)

    summary = (
        raw.groupby(["antibiotic", "model", "model_label", "feature", "gene"], as_index=False)
        .agg(
            selection_frequency=("fold", "nunique"),
            mean_rank=("rank", "mean"),
            mean_abs_importance=("abs_importance", "mean"),
            mean_signed_importance=("importance", "mean"),
        )
        .sort_values(
            ["selection_frequency", "mean_rank", "mean_abs_importance"],
            ascending=[False, True, False],
        )
        .reset_index(drop=True)
    )

    summary["selection_fraction"] = summary["selection_frequency"] / total_folds
    return raw, summary, total_folds


# -----------------------------------------------------------------------------
# Análisis principal
# -----------------------------------------------------------------------------


def run_analysis(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    all_audit = []
    all_model_comparison = []
    all_search_results = []
    all_holdout = []
    all_feature_importances = []
    all_stability = []
    all_stability_raw = []
    selected_rows = []

    start_global = time.time()

    for antibiotic, path in DATASETS.items():
        if not path.exists():
            raise FileNotFoundError(
                f"No encontré {path}. Coloca este script en la misma carpeta que los archivos dataset_*_model_ready.csv."
            )

        print("=" * 80)
        print(f"Analizando {antibiotic}")
        print("=" * 80)

        df = pd.read_csv(path, low_memory=False)
        feature_cols = [c for c in df.columns if c.startswith("amr__")]
        if "y" not in df.columns:
            raise ValueError(f"El archivo {path} no tiene columna 'y'.")
        if not feature_cols:
            raise ValueError(f"No encontré columnas que empiecen con 'amr__' en {path}.")

        X = df[feature_cols].astype(np.int8)
        y = df["y"].astype(int)

        counts = y.value_counts().to_dict()
        n_susceptible = int(counts.get(0, 0))
        n_resistant = int(counts.get(1, 0))
        resistant_fraction = n_resistant / len(df)

        plot_class_balance(antibiotic, y, outdir)

        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=args.test_size,
            stratify=y,
            random_state=RANDOM_STATE,
        )

        train_counts = y_train.value_counts().to_dict()
        train_sus = int(train_counts.get(0, 0))
        train_res = int(train_counts.get(1, 0))
        scale_pos_weight = train_sus / train_res

        all_audit.append(
            {
                "antibiotic": antibiotic,
                "n_samples": len(df),
                "n_features": len(feature_cols),
                "n_susceptible": n_susceptible,
                "n_resistant": n_resistant,
                "resistant_fraction": resistant_fraction,
                "n_train": len(X_train),
                "n_test": len(X_test),
                "train_susceptible": train_sus,
                "train_resistant": train_res,
                "scale_pos_weight_xgboost_train": scale_pos_weight,
                "duplicated_genome_ids": int(df["genome_id"].duplicated().sum()) if "genome_id" in df.columns else np.nan,
            }
        )

        candidate_models = make_candidate_models(scale_pos_weight=scale_pos_weight)
        candidate_models["logistic_regression_balanced"]["n_iter"] = args.n_iter_logreg
        candidate_models["random_forest_balanced"]["n_iter"] = args.n_iter_rf
        candidate_models["xgboost_weighted"]["n_iter"] = args.n_iter_xgb

        cv = RepeatedStratifiedKFold(
            n_splits=args.cv_splits,
            n_repeats=args.cv_repeats,
            random_state=RANDOM_STATE,
        )

        rows_antibiotic = []
        holdout_candidates = []

        for model_family, spec in candidate_models.items():
            print(f"\nAjustando {spec['label']} ({antibiotic})")
            print(f"n_iter = {spec['n_iter']}")

            search = RandomizedSearchCV(
                estimator=spec["estimator"],
                param_distributions=spec["param_dist"],
                n_iter=spec["n_iter"],
                scoring={
                    "balanced_accuracy": "balanced_accuracy",
                    "roc_auc": "roc_auc",
                    "f1": "f1",
                    "recall": "recall",
                },
                refit="balanced_accuracy",
                cv=cv,
                n_jobs=args.n_jobs_search,
                verbose=args.verbose,
                random_state=RANDOM_STATE,
                pre_dispatch="n_jobs",
                return_train_score=False,
                error_score=np.nan,
            )

            t0 = time.time()
            if args.n_jobs_search != 1 and joblib is not None:
                with joblib.parallel_backend("threading"):
                    search.fit(X_train, y_train)
            else:
                search.fit(X_train, y_train)
            elapsed = time.time() - t0

            best = search.best_estimator_
            best_idx = search.best_index_

            cv_row = {
                "antibiotic": antibiotic,
                "model": model_family,
                "model_label": spec["label"],
                "n_iter": spec["n_iter"],
                "cv_splits": args.cv_splits,
                "cv_repeats": args.cv_repeats,
                "fits_requested": spec["n_iter"] * args.cv_splits * args.cv_repeats,
                "search_elapsed_seconds": elapsed,
                "cv_balanced_accuracy_mean": search.cv_results_["mean_test_balanced_accuracy"][best_idx],
                "cv_balanced_accuracy_std": search.cv_results_["std_test_balanced_accuracy"][best_idx],
                "cv_f1_resistant_mean": search.cv_results_["mean_test_f1"][best_idx],
                "cv_f1_resistant_std": search.cv_results_["std_test_f1"][best_idx],
                "cv_roc_auc_mean": search.cv_results_["mean_test_roc_auc"][best_idx],
                "cv_roc_auc_std": search.cv_results_["std_test_roc_auc"][best_idx],
                "cv_recall_resistant_mean": search.cv_results_["mean_test_recall"][best_idx],
                "cv_recall_resistant_std": search.cv_results_["std_test_recall"][best_idx],
                "best_params": json.dumps(search.best_params_, ensure_ascii=False),
            }
            rows_antibiotic.append(cv_row)
            all_model_comparison.append(cv_row)

            cv_results = pd.DataFrame(search.cv_results_)
            cv_results.insert(0, "antibiotic", antibiotic)
            cv_results.insert(1, "model", model_family)
            cv_results.insert(2, "model_label", spec["label"])
            all_search_results.append(cv_results)

            y_pred = best.predict(X_test)
            y_proba = best.predict_proba(X_test)[:, 1]

            metrics, cm = evaluate_predictions(
                np.asarray(y_test),
                np.asarray(y_pred),
                np.asarray(y_proba),
                n_boot=args.n_boot,
            )
            metrics.update(
                {
                    "antibiotic": antibiotic,
                    "model": model_family,
                    "model_label": spec["label"],
                    "n_train": len(X_train),
                    "n_test": len(X_test),
                    "n_features": len(feature_cols),
                    "best_params": json.dumps(search.best_params_, ensure_ascii=False),
                }
            )
            holdout_candidates.append(metrics)

        comp_df = (
            pd.DataFrame(rows_antibiotic)
            .sort_values(
                ["cv_balanced_accuracy_mean", "cv_roc_auc_mean", "cv_recall_resistant_mean"],
                ascending=[False, False, False],
            )
            .reset_index(drop=True)
        )
        plot_model_comparison(antibiotic, comp_df, outdir)

        holdout_df_antibiotic = pd.DataFrame(holdout_candidates)
        plot_metric_panel(antibiotic, holdout_df_antibiotic, outdir)

        selected_family = comp_df.loc[0, "model"]
        selected_label = comp_df.loc[0, "model_label"]
        selected_params = json.loads(comp_df.loc[0, "best_params"])

        holdout_selected = holdout_df_antibiotic[holdout_df_antibiotic["model"] == selected_family].iloc[0].to_dict()
        selected_rows.append(
            {
                **comp_df.loc[0].to_dict(),
                "holdout_accuracy": holdout_selected["accuracy"],
                "holdout_balanced_accuracy": holdout_selected["balanced_accuracy"],
                "holdout_balanced_accuracy_ci95_low": holdout_selected["balanced_accuracy_ci95_low"],
                "holdout_balanced_accuracy_ci95_high": holdout_selected["balanced_accuracy_ci95_high"],
                "holdout_roc_auc": holdout_selected["roc_auc"],
                "holdout_recall_resistant": holdout_selected["recall_resistant"],
                "holdout_specificity_susceptible": holdout_selected["specificity_susceptible"],
                "holdout_false_susceptible": holdout_selected["false_susceptible"],
                "holdout_false_resistant": holdout_selected["false_resistant"],
                "selection_rule": "Mayor exactitud balanceada media en validación cruzada repetida",
            }
        )

        for row in holdout_candidates:
            row["selected_for_final_report"] = row["model"] == selected_family
            all_holdout.append(row)

        # Reentrenar el modelo seleccionado para figuras finales y estabilidad.
        spec_selected = candidate_models[selected_family]
        selected_estimator = clone(spec_selected["estimator"]).set_params(**selected_params)
        selected_estimator.fit(X_train, y_train)

        y_pred_sel = selected_estimator.predict(X_test)
        y_proba_sel = selected_estimator.predict_proba(X_test)[:, 1]
        _, cm_sel = evaluate_predictions(
            np.asarray(y_test),
            np.asarray(y_pred_sel),
            np.asarray(y_proba_sel),
            n_boot=args.n_boot,
        )

        plot_confusion_matrix_es(antibiotic, selected_label, cm_sel, outdir)
        plot_roc_es(
            antibiotic,
            selected_label,
            np.asarray(y_test),
            np.asarray(y_proba_sel),
            roc_auc_score(y_test, y_proba_sel),
            outdir,
        )

        # Modelo final sobre todos los datos solo para importancias finales.
        final_all_data = clone(spec_selected["estimator"]).set_params(**selected_params)
        final_all_data.fit(X, y)
        imp = extract_importance(final_all_data, feature_cols, selected_family)
        imp["antibiotic"] = antibiotic
        imp["model"] = selected_family
        imp["model_label"] = selected_label
        all_feature_importances.append(imp)

        if args.save_models and joblib is not None:
            model_path = outdir / f"modelo_final_{antibiotic}_{selected_family}.joblib"
            joblib.dump(final_all_data, model_path)

        raw_stability, stability, total_folds = summarize_gene_stability(
            estimator=final_all_data,
            model_family=selected_family,
            X=X,
            y=y,
            feature_cols=feature_cols,
            antibiotic=antibiotic,
            n_splits=args.stability_splits,
            n_repeats=args.stability_repeats,
            top_n=args.stability_top_n,
        )
        all_stability_raw.append(raw_stability)
        all_stability.append(stability)
        plot_stable_genes(antibiotic, selected_label, stability, total_folds, outdir, n=15)

    # ---------------------------------------------------------------------
    # Guardado de tablas
    # ---------------------------------------------------------------------

    audit_df = pd.DataFrame(all_audit)
    comparison_df = pd.DataFrame(all_model_comparison)
    search_results_df = pd.concat(all_search_results, ignore_index=True)
    holdout_df = pd.DataFrame(all_holdout)
    selected_df = pd.DataFrame(selected_rows)
    feature_importances_df = pd.concat(all_feature_importances, ignore_index=True)
    stability_df = pd.concat(all_stability, ignore_index=True)
    stability_raw_df = pd.concat(all_stability_raw, ignore_index=True)

    top_genes = []
    for antibiotic in DATASETS:
        top_genes.append(stability_df[stability_df["antibiotic"] == antibiotic].head(20))
    top_genes_df = pd.concat(top_genes, ignore_index=True)

    tables = {
        "resumen_auditoria_datos.csv": audit_df,
        "comparacion_modelos_validacion_cruzada.csv": comparison_df,
        "resultados_completos_randomized_search.csv": search_results_df,
        "metricas_holdout_todos_los_modelos.csv": holdout_df,
        "modelos_seleccionados_resumen.csv": selected_df,
        "importancias_variables_modelos_seleccionados.csv": feature_importances_df,
        "estabilidad_genes_modelos_seleccionados.csv": stability_df,
        "estabilidad_genes_detalle_folds.csv": stability_raw_df,
        "top_genes_estables_modelos_seleccionados.csv": top_genes_df,
    }

    for filename, table in tables.items():
        round_for_export(table).to_csv(outdir / filename, index=False)

    # ---------------------------------------------------------------------
    # Reporte final
    # ---------------------------------------------------------------------

    selected_sections = []
    for _, row in selected_df.iterrows():
        anti = row["antibiotic"]
        model = row["model"]
        genes = stability_df[(stability_df["antibiotic"] == anti) & (stability_df["model"] == model)].head(12).copy()
        genes = genes[["gene", "selection_frequency", "selection_fraction", "mean_rank", "mean_abs_importance"]]

        selected_sections.append(
            f"### {anti} — {row['model_label']}\n\n"
            f"El modelo seleccionado fue **{row['model_label']}**, con exactitud balanceada media de "
            f"{fmt_sig(row['cv_balanced_accuracy_mean'])} en validación cruzada repetida "
            f"(DE = {fmt_sig(row['cv_balanced_accuracy_std'], 2)}). En el conjunto retenido, obtuvo "
            f"exactitud balanceada de {fmt_sig(row['holdout_balanced_accuracy'])} "
            f"(IC95%: {fmt_sig(row['holdout_balanced_accuracy_ci95_low'])}–{fmt_sig(row['holdout_balanced_accuracy_ci95_high'])}), "
            f"ROC-AUC de {fmt_sig(row['holdout_roc_auc'])}, sensibilidad para resistentes de "
            f"{fmt_sig(row['holdout_recall_resistant'])} y especificidad para susceptibles de "
            f"{fmt_sig(row['holdout_specificity_susceptible'])}. Los falsos susceptibles fueron "
            f"{int(row['holdout_false_susceptible'])}.\n\n"
            f"Genes más estables del modelo seleccionado:\n\n{md_table(round_for_export(genes))}\n"
        )

    report = f"""
# Resultados finales optimizados del análisis AMR en *Escherichia coli*

## Enfoque metodológico

Se compararon tres modelos de clasificación supervisada: regresión logística balanceada, bosque aleatorio balanceado y XGBoost ponderado. Para cada antibiótico se separó un conjunto retenido del {fmt_pct(args.test_size)} de los genomas, estratificado por fenotipo. La selección de hiperparámetros se realizó únicamente sobre el conjunto de entrenamiento mediante `RandomizedSearchCV` con validación cruzada repetida y estratificada de {args.cv_splits} particiones × {args.cv_repeats} repeticiones.

La métrica principal de selección fue la exactitud balanceada, porque los datos presentan más genomas susceptibles que resistentes. También se guardaron ROC-AUC, F1 y sensibilidad para resistentes, ya que un falso susceptible es el error más delicado en resistencia antimicrobiana.

## Auditoría de datos

{md_table(round_for_export(audit_df))}

## Comparación de modelos en validación cruzada

{md_table(round_for_export(comparison_df[['antibiotic','model_label','n_iter','fits_requested','cv_balanced_accuracy_mean','cv_balanced_accuracy_std','cv_f1_resistant_mean','cv_roc_auc_mean','cv_recall_resistant_mean','best_params']]))}

## Evaluación en conjunto retenido

{md_table(round_for_export(holdout_df[['antibiotic','model_label','selected_for_final_report','accuracy','balanced_accuracy','balanced_accuracy_ci95_low','balanced_accuracy_ci95_high','f1_resistant','roc_auc','recall_resistant','specificity_susceptible','false_susceptible','false_resistant']]))}

## Modelos seleccionados

{md_table(round_for_export(selected_df[['antibiotic','model_label','cv_balanced_accuracy_mean','cv_balanced_accuracy_std','holdout_balanced_accuracy','holdout_roc_auc','holdout_recall_resistant','holdout_specificity_susceptible','holdout_false_susceptible','holdout_false_resistant','selection_rule']]))}

## Interpretación por antibiótico

{chr(10).join(selected_sections)}

## Conclusión general

La comparación optimizada muestra si la señal de genes AMR permite clasificar de manera consistente genomas resistentes y susceptibles a ciprofloxacin y cefotaxime. El uso de búsqueda aleatoria controlada, validación cruzada repetida, conjunto retenido independiente, intervalos de confianza bootstrap y estabilidad de genes permite presentar estos resultados como análisis final del modelo base optimizado, no como una corrida preliminar.

## Limitaciones

Aunque el análisis es más robusto, todavía debe mencionarse que no incluye una cohorte externa completamente independiente, no incorpora explícitamente mutaciones puntuales no anotadas y no modela estructura poblacional o linajes bacterianos.
""".strip()

    (outdir / "REPORTE_FINAL_OPTIMIZADO.md").write_text(report, encoding="utf-8")

    elapsed_global = time.time() - start_global
    summary = f"""
Análisis final optimizado completado.

Carpeta de salida: {outdir}
Tiempo total: {elapsed_global:.1f} segundos ({elapsed_global/60:.2f} minutos)

Archivos clave:
- REPORTE_FINAL_OPTIMIZADO.md
- modelos_seleccionados_resumen.csv
- comparacion_modelos_validacion_cruzada.csv
- metricas_holdout_todos_los_modelos.csv
- top_genes_estables_modelos_seleccionados.csv
- figuras PNG en español
""".strip()
    (outdir / "RESUMEN_EJECUCION.txt").write_text(summary, encoding="utf-8")

    zip_path = Path(f"{outdir.name}.zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for f in outdir.rglob("*"):
            zf.write(f, f.relative_to(Path(".")))

    print("\n" + summary)
    print(f"ZIP creado: {zip_path.resolve()}")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Análisis final optimizado AMR E. coli con LogReg, RF y XGBoost."
    )
    parser.add_argument("--outdir", default="ecoli_amr_resultados_finales_optimizados")
    parser.add_argument("--test_size", type=float, default=0.20)

    # Parámetros solicitados para búsqueda robusta.
    parser.add_argument("--cv_splits", type=int, default=5)
    parser.add_argument("--cv_repeats", type=int, default=2)
    parser.add_argument("--n_jobs_search", type=int, default=4)
    parser.add_argument("--n_iter_logreg", type=int, default=40)
    parser.add_argument("--n_iter_rf", type=int, default=50)
    parser.add_argument("--n_iter_xgb", type=int, default=60)
    parser.add_argument("--verbose", type=int, default=2)

    # Estabilidad de genes e intervalos de confianza.
    parser.add_argument("--n_boot", type=int, default=500)
    parser.add_argument("--stability_splits", type=int, default=5)
    parser.add_argument("--stability_repeats", type=int, default=3)
    parser.add_argument("--stability_top_n", type=int, default=20)

    # Guardar modelos finales entrenados en todos los datos.
    parser.add_argument("--save_models", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_analysis(args)
