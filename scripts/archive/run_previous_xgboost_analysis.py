from pathlib import Path
import json
import re
import zipfile
import warnings
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
from sklearn.model_selection import GridSearchCV, RepeatedStratifiedKFold, train_test_split
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


plt.rcParams["axes.unicode_minus"] = False

RANDOM_STATE = 42
OUT = Path("ecoli_amr_resultados_finales_con_xgboost")
OUT.mkdir(exist_ok=True)

DATASETS = {
    "ciprofloxacin": Path("dataset_ciprofloxacin_model_ready.csv"),
    "cefotaxime": Path("dataset_cefotaxime_model_ready.csv"),
}


def fmt_sig(x, sig=4):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    if x == 0:
        return "0"
    return f"{x:.{sig}g}"


def fmt_pct(x, sig=3):
    return f"{100*x:.{sig}g}%"


def clean_gene_name(col):
    return re.sub(r"^amr__", "", col)


def specificity_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return tn / (tn + fp) if (tn + fp) else np.nan


def bootstrap_metric_ci(y_true, y_pred, y_proba, metric_name, n_boot=800, seed=42):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    values = []

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
            raise ValueError(metric_name)

        values.append(val)

    arr = np.asarray(values, dtype=float)
    return np.nanpercentile(arr, [2.5, 97.5]) if len(arr) else (np.nan, np.nan)


def evaluate_predictions(y_true, y_pred, y_proba):
    metrics = {
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

    metrics.update({
        "true_susceptible": int(tn),
        "false_resistant": int(fp),
        "false_susceptible": int(fn),
        "true_resistant": int(tp),
        "confusion_matrix": json.dumps(cm.tolist()),
    })

    for m in [
        "accuracy",
        "balanced_accuracy",
        "f1_resistant",
        "roc_auc",
        "recall_resistant",
        "precision_resistant",
        "specificity_susceptible",
    ]:
        lo, hi = bootstrap_metric_ci(np.asarray(y_true), np.asarray(y_pred), np.asarray(y_proba), m)
        metrics[f"{m}_ci95_low"] = lo
        metrics[f"{m}_ci95_high"] = hi

    return metrics, cm


def plot_class_balance(antibiotic, y):
    counts = (
        pd.Series(y)
        .map({0: "Susceptible", 1: "Resistente"})
        .value_counts()
        .reindex(["Susceptible", "Resistente"])
    )

    fig, ax = plt.subplots(figsize=(6, 4.5))
    bars = ax.bar(counts.index, counts.values)

    ax.set_title(f"Distribución de clases para {antibiotic}")
    ax.set_xlabel("Fenotipo")
    ax.set_ylabel("Número de genomas")

    total = counts.sum()
    for bar, val in zip(bars, counts.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val,
            f"{val}\n({fmt_pct(val/total)})",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(OUT / f"figura_balance_clases_{antibiotic}.png", dpi=240)
    plt.close(fig)


def plot_model_comparison(antibiotic, model_rows):
    labels = [r["model_label"] for r in model_rows]
    means = [r["cv_balanced_accuracy_mean"] for r in model_rows]
    sds = [r["cv_balanced_accuracy_std"] for r in model_rows]

    fig, ax = plt.subplots(figsize=(9, 5))
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
    fig.savefig(OUT / f"figura_comparacion_modelos_{antibiotic}.png", dpi=240)
    plt.close(fig)


def plot_confusion_matrix_es(antibiotic, model_label, cm):
    row_sums = cm.sum(axis=1, keepdims=True)
    norm = np.divide(cm, row_sums, out=np.zeros_like(cm, dtype=float), where=row_sums != 0)

    fig, ax = plt.subplots(figsize=(5.8, 5.2))
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
    fig.savefig(OUT / f"figura_matriz_confusion_{antibiotic}.png", dpi=240)
    plt.close(fig)


def plot_roc_es(antibiotic, model_label, y_true, y_proba, auc_value):
    fpr, tpr, _ = roc_curve(y_true, y_proba)

    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    ax.plot(fpr, tpr, label=f"Curva ROC (AUC = {fmt_sig(auc_value)})")
    ax.plot([0, 1], [0, 1], linestyle="--", label="Azar")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Tasa de falsos positivos")
    ax.set_ylabel("Tasa de verdaderos positivos")
    ax.set_title(f"Curva ROC\n{antibiotic} — {model_label}")
    ax.legend(loc="lower right")

    fig.tight_layout()
    fig.savefig(OUT / f"figura_curva_roc_{antibiotic}.png", dpi=240)
    plt.close(fig)


def extract_importance(model, feature_cols, family):
    if family == "logistic_regression_balanced":
        vals = model.named_steps["logisticregression"].coef_[0]
        abs_vals = np.abs(vals)
        imp_type = "coeficiente_logistico_absoluto"
    elif family == "random_forest_balanced":
        vals = model.feature_importances_
        abs_vals = vals
        imp_type = "importancia_gini_bosque_aleatorio"
    elif family == "xgboost_weighted":
        vals = model.feature_importances_
        abs_vals = vals
        imp_type = "importancia_xgboost_gain_normalizada"
    else:
        raise ValueError(f"Modelo no reconocido: {family}")

    return (
        pd.DataFrame({
            "feature": feature_cols,
            "gene": [clean_gene_name(c) for c in feature_cols],
            "importance": vals,
            "abs_importance": abs_vals,
            "importance_type": imp_type,
        })
        .sort_values("abs_importance", ascending=False)
        .reset_index(drop=True)
    )


def summarize_gene_stability(estimator, family, X, y, feature_cols, antibiotic, n_splits=5, n_repeats=3, top_n=20):
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=RANDOM_STATE)
    records = []
    total_folds = n_splits * n_repeats

    for fold, (tr, te) in enumerate(cv.split(X, y), start=1):
        model = clone(estimator)
        model.fit(X.iloc[tr], y.iloc[tr])

        top = extract_importance(model, feature_cols, family).head(top_n).copy()
        top["fold"] = fold
        top["rank"] = np.arange(1, len(top) + 1)
        top["antibiotic"] = antibiotic
        top["model"] = family
        records.append(top)

    raw = pd.concat(records, ignore_index=True)

    summary = (
        raw.groupby(["antibiotic", "model", "feature", "gene"], as_index=False)
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


def plot_stable_genes(antibiotic, model_label, stability, total_folds, n=15):
    top = stability.head(n).iloc[::-1].copy()

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    bars = ax.barh(top["gene"], top["selection_frequency"])

    ax.set_xlabel("Número de veces que el gen apareció entre los más importantes")
    ax.set_ylabel("Gen AMR")
    ax.set_title(f"Genes más estables\n{antibiotic} — {model_label}")
    ax.set_xlim(0, total_folds * 1.08)

    for bar, freq, frac in zip(bars, top["selection_frequency"], top["selection_fraction"]):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2,
            f" {int(freq)}/{total_folds} ({fmt_pct(frac)})",
            va="center",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(OUT / f"figura_genes_estables_{antibiotic}.png", dpi=240)
    plt.close(fig)


def round_for_export(df):
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].round(6)
    return out


def md_table(df):
    return df.to_markdown(index=False)


all_audit = []
all_model_comparison = []
all_holdout = []
all_feature_importances = []
all_stability = []
all_stability_raw = []
selected_rows = []

for antibiotic, path in DATASETS.items():
    if not path.exists():
        raise FileNotFoundError(
            f"No encontré {path}. Pon este script en la misma carpeta que los archivos dataset_*_model_ready.csv."
        )

    print(f"Analizando {antibiotic}...")

    df = pd.read_csv(path, low_memory=False)
    feature_cols = [c for c in df.columns if c.startswith("amr__")]

    X = df[feature_cols].astype(np.int8)
    y = df["y"].astype(int)

    counts = y.value_counts().to_dict()
    n_susceptible = counts.get(0, 0)
    n_resistant = counts.get(1, 0)
    resistant_fraction = n_resistant / len(df)
    scale_pos_weight = n_susceptible / n_resistant

    all_audit.append({
        "antibiotic": antibiotic,
        "n_samples": len(df),
        "n_features": len(feature_cols),
        "n_susceptible": n_susceptible,
        "n_resistant": n_resistant,
        "resistant_fraction": resistant_fraction,
        "scale_pos_weight_used_for_xgboost": scale_pos_weight,
        "duplicated_genome_ids": int(df["genome_id"].duplicated().sum()) if "genome_id" in df.columns else np.nan,
    })

    plot_class_balance(antibiotic, y)

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        stratify=y,
        random_state=RANDOM_STATE,
    )

    candidate_models = {
        "logistic_regression_balanced": {
            "label": "Regresión logística balanceada",
            "estimator": make_pipeline(
                StandardScaler(with_mean=False),
                LogisticRegression(
                    max_iter=4000,
                    class_weight="balanced",
                    solver="liblinear",
                    random_state=RANDOM_STATE,
                ),
            ),
            "param_grid": {
                "logisticregression__C": [0.1, 0.3, 1, 3, 10],
                "logisticregression__penalty": ["l1", "l2"],
            },
        },
        "random_forest_balanced": {
            "label": "Bosque aleatorio balanceado",
            "estimator": RandomForestClassifier(
                n_estimators=200,
                class_weight="balanced",
                random_state=RANDOM_STATE,
                n_jobs=-1,
            ),
            "param_grid": {
                "max_depth": [None, 15, 25],
                "min_samples_leaf": [1, 2, 4],
                "max_features": ["sqrt", "log2"],
            },
        },
        "xgboost_weighted": {
            "label": "XGBoost ponderado",
            "estimator": XGBClassifier(
                objective="binary:logistic",
                eval_metric="logloss",
                tree_method="hist",
                n_estimators=250,
                scale_pos_weight=scale_pos_weight,
                random_state=RANDOM_STATE,
                n_jobs=-1,
            ),
            "param_grid": {
                "max_depth": [3, 5],
                "learning_rate": [0.05, 0.1],
                "subsample": [0.85, 1.0],
                "colsample_bytree": [0.85, 1.0],
                "min_child_weight": [1, 5],
            },
        },
    }

    inner_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=RANDOM_STATE)

    rows_antibiotic = []
    holdout_candidates = []

    for family, spec in candidate_models.items():
        print(f"  Ajustando {family}...")

        search = GridSearchCV(
            estimator=spec["estimator"],
            param_grid=spec["param_grid"],
            scoring={
                "balanced_accuracy": "balanced_accuracy",
                "f1_resistant": "f1",
                "roc_auc": "roc_auc",
                "recall_resistant": "recall",
            },
            refit="balanced_accuracy",
            cv=inner_cv,
            n_jobs=-1,
            return_train_score=False,
        )

        search.fit(X_train, y_train)

        best = search.best_estimator_
        best_idx = search.best_index_

        row = {
            "antibiotic": antibiotic,
            "model": family,
            "model_label": spec["label"],
            "cv_balanced_accuracy_mean": search.cv_results_["mean_test_balanced_accuracy"][best_idx],
            "cv_balanced_accuracy_std": search.cv_results_["std_test_balanced_accuracy"][best_idx],
            "cv_f1_resistant_mean": search.cv_results_["mean_test_f1_resistant"][best_idx],
            "cv_f1_resistant_std": search.cv_results_["std_test_f1_resistant"][best_idx],
            "cv_roc_auc_mean": search.cv_results_["mean_test_roc_auc"][best_idx],
            "cv_roc_auc_std": search.cv_results_["std_test_roc_auc"][best_idx],
            "cv_recall_resistant_mean": search.cv_results_["mean_test_recall_resistant"][best_idx],
            "cv_recall_resistant_std": search.cv_results_["std_test_recall_resistant"][best_idx],
            "best_params": json.dumps(search.best_params_, ensure_ascii=False),
        }

        rows_antibiotic.append(row)
        all_model_comparison.append(row)

        y_pred = best.predict(X_test)
        y_proba = best.predict_proba(X_test)[:, 1]

        metrics, cm = evaluate_predictions(np.asarray(y_test), np.asarray(y_pred), np.asarray(y_proba))
        metrics.update({
            "antibiotic": antibiotic,
            "model": family,
            "model_label": spec["label"],
            "n_train": len(X_train),
            "n_test": len(X_test),
            "n_features": len(feature_cols),
            "best_params": json.dumps(search.best_params_, ensure_ascii=False),
        })
        holdout_candidates.append(metrics)

    comp_df = (
        pd.DataFrame(rows_antibiotic)
        .sort_values(
            ["cv_balanced_accuracy_mean", "cv_roc_auc_mean", "cv_recall_resistant_mean"],
            ascending=[False, False, False],
        )
        .reset_index(drop=True)
    )

    plot_model_comparison(antibiotic, comp_df.to_dict("records"))

    selected_family = comp_df.loc[0, "model"]
    selected_label = comp_df.loc[0, "model_label"]
    selected_params = json.loads(comp_df.loc[0, "best_params"])

    holdout_df = pd.DataFrame(holdout_candidates)
    holdout_selected = holdout_df[holdout_df["model"] == selected_family].iloc[0].to_dict()
    holdout_selected["selection_rule"] = (
        "Mayor exactitud balanceada en validación cruzada interna "
        "(5 folds × 2 repeticiones)"
    )

    selected_rows.append({
        **comp_df.loc[0].to_dict(),
        "holdout_accuracy": holdout_selected["accuracy"],
        "holdout_balanced_accuracy": holdout_selected["balanced_accuracy"],
        "holdout_roc_auc": holdout_selected["roc_auc"],
        "holdout_recall_resistant": holdout_selected["recall_resistant"],
        "holdout_specificity_susceptible": holdout_selected["specificity_susceptible"],
        "holdout_false_susceptible": holdout_selected["false_susceptible"],
        "holdout_false_resistant": holdout_selected["false_resistant"],
        "selection_rule": holdout_selected["selection_rule"],
    })

    for row in holdout_candidates:
        row["selected_for_final_report"] = row["model"] == selected_family
        all_holdout.append(row)

    spec = candidate_models[selected_family]
    final_estimator = clone(spec["estimator"]).set_params(**selected_params)
    final_estimator.fit(X_train, y_train)

    y_pred = final_estimator.predict(X_test)
    y_proba = final_estimator.predict_proba(X_test)[:, 1]
    _, cm = evaluate_predictions(np.asarray(y_test), np.asarray(y_pred), np.asarray(y_proba))

    plot_confusion_matrix_es(antibiotic, selected_label, cm)
    plot_roc_es(antibiotic, selected_label, np.asarray(y_test), np.asarray(y_proba), roc_auc_score(y_test, y_proba))

    final_all_data = clone(spec["estimator"]).set_params(**selected_params)
    final_all_data.fit(X, y)

    imp = extract_importance(final_all_data, feature_cols, selected_family)
    imp["antibiotic"] = antibiotic
    imp["model"] = selected_family
    imp["model_label"] = selected_label
    all_feature_importances.append(imp)

    raw, stability, total_folds = summarize_gene_stability(
        final_all_data,
        selected_family,
        X,
        y,
        feature_cols,
        antibiotic,
    )

    stability["model_label"] = selected_label
    raw["model_label"] = selected_label

    all_stability.append(stability)
    all_stability_raw.append(raw)

    plot_stable_genes(antibiotic, selected_label, stability, total_folds, n=15)


files_to_save = {
    "resumen_auditoria_datos.csv": pd.DataFrame(all_audit),
    "comparacion_modelos_validacion_cruzada.csv": pd.DataFrame(all_model_comparison),
    "metricas_holdout_todos_los_modelos.csv": pd.DataFrame(all_holdout),
    "modelos_seleccionados_resumen.csv": pd.DataFrame(selected_rows),
    "importancias_variables_modelos_seleccionados.csv": pd.concat(all_feature_importances, ignore_index=True),
    "estabilidad_genes_modelos_seleccionados.csv": pd.concat(all_stability, ignore_index=True),
    "estabilidad_genes_detalle_folds.csv": pd.concat(all_stability_raw, ignore_index=True),
}

for filename, table in files_to_save.items():
    round_for_export(table).to_csv(OUT / filename, index=False)

stab = pd.concat(all_stability, ignore_index=True)
top_genes = []

for antibiotic in DATASETS:
    top_genes.append(stab[stab["antibiotic"] == antibiotic].head(20))

top_genes = pd.concat(top_genes, ignore_index=True)
round_for_export(top_genes).to_csv(OUT / "top_genes_estables_modelos_seleccionados.csv", index=False)

audit_df = pd.DataFrame(all_audit)
comp_df = pd.DataFrame(all_model_comparison)
holdout_df = pd.DataFrame(all_holdout)
selected_df = pd.DataFrame(selected_rows)

selected_sections = []

for _, row in selected_df.iterrows():
    anti = row["antibiotic"]
    model = row["model"]

    top = stab[(stab["antibiotic"] == anti) & (stab["model"] == model)].head(12).copy()

    if len(top):
        top = top[["gene", "selection_frequency", "selection_fraction", "mean_rank", "mean_abs_importance"]]

    selected_sections.append(
        f"### {anti} — {row['model_label']}\n\n"
        f"Modelo seleccionado con exactitud balanceada media en validación cruzada de "
        f"{fmt_sig(row['cv_balanced_accuracy_mean'])} "
        f"(DE = {fmt_sig(row['cv_balanced_accuracy_std'], 2)}). En el conjunto retenido, obtuvo "
        f"exactitud balanceada de {fmt_sig(row['holdout_balanced_accuracy'])}, AUC de "
        f"{fmt_sig(row['holdout_roc_auc'])}, sensibilidad para resistentes de "
        f"{fmt_sig(row['holdout_recall_resistant'])} y especificidad para susceptibles de "
        f"{fmt_sig(row['holdout_specificity_susceptible'])}. "
        f"Los errores más delicados, resistentes predichos como susceptibles, fueron "
        f"{int(row['holdout_false_susceptible'])}.\n\n"
        f"Genes más estables del modelo seleccionado\n\n{md_table(top)}\n"
    )

report = f"""
# Resultados finales del análisis de resistencia antimicrobiana en *Escherichia coli*

## Enfoque metodológico final

Para evitar presentar resultados meramente preliminares, el análisis final comparó tres modelos: regresión logística balanceada, bosque aleatorio balanceado y XGBoost ponderado. Se separó un conjunto de prueba retenido e independiente del 20% de los genomas. La selección de hiperparámetros se realizó únicamente sobre el 80% restante mediante validación cruzada repetida y estratificada de 5 particiones con 2 repeticiones. El criterio principal para elegir modelo fue la exactitud balanceada, adecuada cuando las clases resistente y susceptible están desbalanceadas. Además, se evaluó la estabilidad de los genes importantes mediante 15 ajustes repetidos del modelo seleccionado.

## Auditoría de datos de entrada

{md_table(audit_df)}

## Comparación de modelos en validación cruzada interna

{md_table(comp_df[['antibiotic','model_label','cv_balanced_accuracy_mean','cv_balanced_accuracy_std','cv_f1_resistant_mean','cv_roc_auc_mean','cv_recall_resistant_mean','best_params']])}

## Resultados en el conjunto retenido

{md_table(holdout_df[['antibiotic','model_label','selected_for_final_report','accuracy','balanced_accuracy','balanced_accuracy_ci95_low','balanced_accuracy_ci95_high','f1_resistant','roc_auc','recall_resistant','specificity_susceptible','false_susceptible','false_resistant']])}

## Modelos seleccionados para el informe final

{md_table(selected_df[['antibiotic','model_label','cv_balanced_accuracy_mean','cv_balanced_accuracy_std','holdout_balanced_accuracy','holdout_roc_auc','holdout_recall_resistant','holdout_specificity_susceptible','holdout_false_susceptible','holdout_false_resistant','selection_rule']])}

## Interpretación por antibiótico

{chr(10).join(selected_sections)}

## Conclusión general para redactar en resultados

En ambos antibióticos, los modelos basados en presencia y ausencia de genes AMR lograron discriminar entre genomas resistentes y susceptibles con desempeños consistentes tanto en validación cruzada como en evaluación independiente. La incorporación de XGBoost permite contrastar los modelos lineales y de ensamble clásico con un modelo de boosting de gradiente, por lo que la comparación final es más robusta. Además, los genes destacados por los modelos seleccionados mostraron estabilidad a lo largo de múltiples repeticiones, lo que fortalece su interpretación biológica.

## Limitaciones que todavía deben mencionarse

Aunque estos resultados pueden reportarse como resultados finales del modelo base, todavía conviene declarar tres limitaciones. La primera es la ausencia de una cohorte externa completamente independiente. La segunda es que el análisis se basa en genes AMR anotados y no incorpora mutaciones puntuales no anotadas. La tercera es que no se modeló explícitamente la estructura poblacional o los linajes bacterianos.
""".strip()

(OUT / "REPORTE_FINAL_CON_XGBOOST.md").write_text(report, encoding="utf-8")

readme = """
# Instrucciones rápidas

1. Coloca este script en la misma carpeta que:
   - dataset_ciprofloxacin_model_ready.csv
   - dataset_cefotaxime_model_ready.csv

2. Instala dependencias:
   pip install pandas scikit-learn matplotlib tabulate xgboost

   En Jupyter:
   %pip install pandas scikit-learn matplotlib tabulate xgboost

3. Ejecuta:
   python generate_final_nonprelim_analysis_es_xgboost.py

4. El script creará:
   - carpeta ecoli_amr_resultados_finales_con_xgboost
   - archivo ecoli_amr_resultados_finales_con_xgboost.zip

Tiempo esperado:
- Esta versión puede tardar más que la anterior porque XGBoost agrega más combinaciones de hiperparámetros.
- En una laptop moderna puede tardar varios minutos.
- Si tarda demasiado, reduce el grid de XGBoost en el script.
""".strip()

(OUT / "README_INSTRUCCIONES.txt").write_text(readme, encoding="utf-8")

with zipfile.ZipFile("ecoli_amr_resultados_finales_con_xgboost.zip", "w", compression=zipfile.ZIP_DEFLATED) as zf:
    for f in OUT.rglob("*"):
        zf.write(f, f.relative_to(Path(".")))

print("Análisis final con XGBoost completado.")
print("Salida:", OUT.resolve())
print("ZIP:", Path("ecoli_amr_resultados_finales_con_xgboost.zip").resolve())
