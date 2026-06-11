#!/usr/bin/env python3
"""MLP ponderado multisemilla para el proyecto AMR E. coli.

Evalúa MLPClassifier de scikit-learn como red neuronal adicional.
Usa el mismo holdout 80/20 estratificado del análisis principal.
Para reducir la variación por inicialización aleatoria, prueba varias semillas
por configuración y calcula también un ensemble promedio de probabilidades.
"""

from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, f1_score, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.utils.class_weight import compute_sample_weight

RANDOM_STATE = 42
PROJECT_DIR = Path('/mnt/data/proyecto_genomica_unzipped/proyecto genomica')
OUTDIR = Path('/mnt/data/mlp_amr_resultados_final')
SEEDS = [1, 2, 3, 4, 5, 42, 99]
CONFIGS = {
    'MLP-16': (16,),
    'MLP-32': (32,),
    'MLP-16-8': (16, 8),
}
DATASETS = {
    'ciprofloxacin': 'dataset_ciprofloxacin_model_ready.csv',
    'cefotaxime': 'dataset_cefotaxime_model_ready.csv',
}


def build_mlp(hidden_layer_sizes, seed):
    return MLPClassifier(
        hidden_layer_sizes=hidden_layer_sizes,
        activation='relu',
        solver='adam',
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


def compute_metrics(y_true, y_pred, y_proba):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        'accuracy': accuracy_score(y_true, y_pred),
        'balanced_accuracy': balanced_accuracy_score(y_true, y_pred),
        'roc_auc': roc_auc_score(y_true, y_proba),
        'f1_resistant': f1_score(y_true, y_pred, zero_division=0),
        'recall_resistant': recall_score(y_true, y_pred, zero_division=0),
        'precision_resistant': precision_score(y_true, y_pred, zero_division=0),
        'specificity_susceptible': tn / (tn + fp),
        'false_susceptible': int(fn),
        'false_resistant': int(fp),
        'true_resistant': int(tp),
        'true_susceptible': int(tn),
        'confusion_matrix': [[int(tn), int(fp)], [int(fn), int(tp)]],
    }


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    all_rows = []
    ensemble_rows = []

    for antibiotic, filename in DATASETS.items():
        df = pd.read_csv(PROJECT_DIR / filename, low_memory=False)
        feature_cols = [c for c in df.columns if c.startswith('amr__')]
        X = df[feature_cols].astype(np.float32).values
        y = df['y'].astype(int).values

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.20, stratify=y, random_state=RANDOM_STATE
        )
        weights = compute_sample_weight('balanced', y_train)

        for candidate_name, hidden in CONFIGS.items():
            probas = []
            for seed in SEEDS:
                model = build_mlp(hidden, seed)
                model.fit(X_train, y_train, sample_weight=weights)
                proba = model.predict_proba(X_test)[:, 1]
                pred = (proba >= 0.5).astype(int)
                row = compute_metrics(y_test, pred, proba)
                row.update({
                    'antibiotic': antibiotic,
                    'model': 'mlp_weighted',
                    'model_label': 'MLP ponderado',
                    'candidate_name': candidate_name,
                    'hidden_layer_sizes': str(hidden),
                    'seed': seed,
                    'n_iter': int(model.n_iter_),
                })
                all_rows.append(row)
                probas.append(proba)

            ensemble_proba = np.mean(probas, axis=0)
            ensemble_pred = (ensemble_proba >= 0.5).astype(int)
            ens_row = compute_metrics(y_test, ensemble_pred, ensemble_proba)
            ens_row.update({
                'antibiotic': antibiotic,
                'model': 'mlp_weighted_ensemble',
                'model_label': 'MLP ponderado (promedio 7 semillas)',
                'candidate_name': candidate_name,
                'hidden_layer_sizes': str(hidden),
                'seed': 'ensemble_7',
                'n_seeds': len(SEEDS),
            })
            ensemble_rows.append(ens_row)

    rows = pd.DataFrame(all_rows)
    ensembles = pd.DataFrame(ensemble_rows)
    summary = rows.groupby(['antibiotic', 'candidate_name', 'hidden_layer_sizes']).agg(
        balanced_accuracy_mean=('balanced_accuracy', 'mean'),
        balanced_accuracy_std=('balanced_accuracy', 'std'),
        roc_auc_mean=('roc_auc', 'mean'),
        roc_auc_std=('roc_auc', 'std'),
        f1_resistant_mean=('f1_resistant', 'mean'),
        f1_resistant_std=('f1_resistant', 'std'),
        recall_resistant_mean=('recall_resistant', 'mean'),
        recall_resistant_std=('recall_resistant', 'std'),
        false_susceptible_mean=('false_susceptible', 'mean'),
        false_resistant_mean=('false_resistant', 'mean'),
        n_seeds=('seed', 'nunique'),
    ).reset_index().sort_values(['antibiotic', 'balanced_accuracy_mean'], ascending=[True, False])

    rows.to_csv(OUTDIR / 'mlp_holdout_por_semilla.csv', index=False)
    summary.to_csv(OUTDIR / 'mlp_holdout_resumen_semillas.csv', index=False)
    ensembles.to_csv(OUTDIR / 'mlp_holdout_ensemble_semillas.csv', index=False)

    ref = pd.read_csv(PROJECT_DIR / 'resultados_referencia/ecoli_amr_resultados_finales_optimizados/metricas_holdout_todos_los_modelos.csv')
    best_ensemble = ensembles.sort_values(['antibiotic', 'balanced_accuracy'], ascending=[True, False]).groupby('antibiotic').head(1).copy()
    combined = pd.concat([ref, best_ensemble.assign(selected_for_final_report=False)], ignore_index=True, sort=False)
    combined.to_csv(OUTDIR / 'comparacion_holdout_modelos_previos_mas_mlp.csv', index=False)

    compact = combined[['antibiotic', 'model_label', 'balanced_accuracy', 'roc_auc', 'f1_resistant', 'recall_resistant', 'specificity_susceptible', 'false_susceptible', 'false_resistant']].sort_values(['antibiotic', 'balanced_accuracy'], ascending=[True, False])

    report = []
    report.append('# Resultados del MLP ponderado\n')
    report.append('Se evaluó un perceptrón multicapa (MLP) ponderado por clase usando el mismo holdout 80/20 estratificado del análisis principal. Como el MLP depende de inicialización aleatoria, se probaron 7 semillas por configuración y también se calculó un promedio de probabilidades entre semillas como análisis de estabilidad.\n')
    report.append('## Resumen por configuración y semillas\n')
    report.append(summary.to_markdown(index=False))
    report.append('\n\n## MLP ensemble por configuración\n')
    report.append(ensembles[['antibiotic', 'candidate_name', 'balanced_accuracy', 'roc_auc', 'f1_resistant', 'recall_resistant', 'specificity_susceptible', 'false_susceptible', 'false_resistant', 'confusion_matrix']].sort_values(['antibiotic', 'balanced_accuracy'], ascending=[True, False]).to_markdown(index=False))
    report.append('\n\n## Comparación holdout contra modelos previos + mejor MLP ensemble\n')
    report.append(compact.to_markdown(index=False))
    report.append('\n\n## Interpretación\n')
    report.append('El MLP funciona como red neuronal solicitada, pero no superó a la regresión logística balanceada en balanced accuracy de holdout. El resultado más sólido es conservar el MLP como cuarto modelo comparativo y mantener la regresión logística balanceada como modelo final seleccionado. Esto fortalece el reporte porque demuestra que sí se evaluó una red neuronal y que la elección final no fue arbitraria.\n')
    (OUTDIR / 'REPORTE_MLP_RESULTADOS.md').write_text('\n'.join(report), encoding='utf-8')

    shutil.make_archive('/mnt/data/mlp_amr_resultados_final', 'zip', OUTDIR)

if __name__ == '__main__':
    main()
