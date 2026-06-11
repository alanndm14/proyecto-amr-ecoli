from pathlib import Path
import json, re, zipfile, warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
                             recall_score, precision_score, confusion_matrix, ConfusionMatrixDisplay,
                             RocCurveDisplay, make_scorer)
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV, cross_validate
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

OUT = Path('ecoli_amr_robust_final_analysis')
OUT.mkdir(exist_ok=True)
DATASETS = {
    'ciprofloxacin': Path('dataset_ciprofloxacin_model_ready.csv'),
    'cefotaxime': Path('dataset_cefotaxime_model_ready.csv')
}
RANDOM_STATE = 42


def clean_gene_name(col):
    return re.sub(r'^amr__', '', col)


def specificity_score(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return tn / (tn + fp) if (tn + fp) else np.nan


def extract_importance(model, feature_cols, model_name):
    if model_name == 'logistic_regression_balanced':
        vals = model.named_steps['logisticregression'].coef_[0]
        abs_vals = np.abs(vals)
        itype = 'absolute_logistic_coefficient'
    else:
        vals = model.feature_importances_
        abs_vals = vals
        itype = 'random_forest_gini_importance'
    return pd.DataFrame({
        'feature': feature_cols,
        'gene': [clean_gene_name(c) for c in feature_cols],
        'importance': vals,
        'abs_importance': abs_vals,
        'importance_type': itype
    }).sort_values('abs_importance', ascending=False).reset_index(drop=True)


def plot_confusion(antibiotic, model_name, cm):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ConfusionMatrixDisplay(cm, display_labels=['Susceptible', 'Resistant']).plot(ax=ax, values_format='d')
    ax.set_title(f'{antibiotic} - {model_name}')
    fig.tight_layout()
    fig.savefig(OUT / f'final_confusion_matrix_{antibiotic}_{model_name}.png', dpi=220)
    plt.close(fig)


def plot_roc(antibiotic, model_name, y_test, proba):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    RocCurveDisplay.from_predictions(y_test, proba, ax=ax)
    ax.set_title(f'ROC - {antibiotic} - {model_name}')
    fig.tight_layout()
    fig.savefig(OUT / f'final_roc_curve_{antibiotic}_{model_name}.png', dpi=220)
    plt.close(fig)


def plot_stable_genes(antibiotic, model_name, stability, n=20):
    top = stability.head(n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 7))
    ax.barh(top['gene'], top['selection_frequency'])
    ax.set_xlabel('Frecuencia en folds')
    ax.set_ylabel('Gen AMR')
    ax.set_title(f'Genes estables - {antibiotic} - {model_name}')
    fig.tight_layout()
    fig.savefig(OUT / f'stable_top{n}_genes_{antibiotic}_{model_name}.png', dpi=220)
    plt.close(fig)


def plot_cv(antibiotic, cvsum):
    sub = cvsum[cvsum['antibiotic'] == antibiotic].copy()
    sub['label'] = sub['model'].str.replace('_', ' ')
    metrics = ['balanced_accuracy_mean', 'f1_resistant_mean', 'recall_resistant_mean', 'roc_auc_mean']
    labels = ['Balanced accuracy', 'F1 resistant', 'Recall resistant', 'ROC-AUC']
    x = np.arange(len(sub))
    width = 0.18
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for i, (metric, label) in enumerate(zip(metrics, labels)):
        ax.bar(x + (i - 1.5) * width, sub[metric], width, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(sub['label'], rotation=20, ha='right')
    ax.set_ylim(0, 1)
    ax.set_ylabel('Promedio en validación cruzada')
    ax.set_title(f'Validación cruzada - {antibiotic}')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / f'cv_metrics_{antibiotic}.png', dpi=220)
    plt.close(fig)


def gene_stability(estimator, model_name, X, y, feature_cols, antibiotic, top_n=20):
    cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE)
    records = []
    for fold, (tr, te) in enumerate(cv.split(X, y), start=1):
        m = clone(estimator)
        m.fit(X.iloc[tr], y.iloc[tr])
        imp = extract_importance(m, feature_cols, model_name).head(top_n).copy()
        imp['fold'] = fold
        imp['rank'] = np.arange(1, len(imp) + 1)
        imp['antibiotic'] = antibiotic
        imp['model'] = model_name
        records.append(imp)
    raw = pd.concat(records, ignore_index=True)
    stability = (raw.groupby(['antibiotic', 'model', 'feature', 'gene'], as_index=False)
                 .agg(selection_frequency=('fold', 'nunique'),
                      mean_rank=('rank', 'mean'),
                      mean_abs_importance=('abs_importance', 'mean'),
                      mean_signed_importance=('importance', 'mean'))
                 .sort_values(['selection_frequency', 'mean_rank', 'mean_abs_importance'], ascending=[False, True, False])
                 .reset_index(drop=True))
    return raw, stability


scoring = {
    'accuracy': 'accuracy',
    'balanced_accuracy': 'balanced_accuracy',
    'f1_resistant': make_scorer(f1_score, pos_label=1),
    'roc_auc': 'roc_auc',
    'recall_resistant': make_scorer(recall_score, pos_label=1),
    'precision_resistant': make_scorer(precision_score, pos_label=1, zero_division=0),
    'specificity_susceptible': make_scorer(specificity_score)
}

all_audit, all_tuning, all_holdout, all_cv, all_imp, all_stab, all_stab_raw = [], [], [], [], [], [], []

for antibiotic, path in DATASETS.items():
    if not path.exists():
        raise FileNotFoundError(f'No encontré {path}. Coloca el script en la misma carpeta que los datasets model_ready.')
    print(f'Analizando {antibiotic}...')
    df = pd.read_csv(path, low_memory=False)
    feature_cols = [c for c in df.columns if c.startswith('amr__')]
    X = df[feature_cols].astype(np.int8)
    y = df['y'].astype(int)
    counts = y.value_counts().to_dict()
    all_audit.append({
        'antibiotic': antibiotic,
        'n_samples': len(df),
        'n_features': len(feature_cols),
        'n_susceptible': counts.get(0, 0),
        'n_resistant': counts.get(1, 0),
        'resistant_fraction': counts.get(1, 0) / len(df),
        'duplicated_genome_ids': int(df['genome_id'].duplicated().sum()) if 'genome_id' in df.columns else np.nan
    })

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2, stratify=y, random_state=RANDOM_STATE)
    cv_tune = StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE)

    specs = {
        'logistic_regression_balanced': (
            make_pipeline(StandardScaler(with_mean=False),
                          LogisticRegression(max_iter=3000, class_weight='balanced', solver='liblinear', random_state=RANDOM_STATE)),
            {'logisticregression__C': [0.1, 1, 10], 'logisticregression__penalty': ['l1', 'l2']}
        ),
        'random_forest_balanced': (
            RandomForestClassifier(n_estimators=60, class_weight='balanced', random_state=RANDOM_STATE, n_jobs=-1),
            {'min_samples_leaf': [1, 2], 'max_features': ['sqrt'], 'max_depth': [None]}
        )
    }

    for model_name, (estimator, param_grid) in specs.items():
        print(f'  Modelo {model_name}')
        search = GridSearchCV(estimator, param_grid, scoring='balanced_accuracy', cv=cv_tune, n_jobs=-1, refit=True)
        search.fit(X_train, y_train)
        best = search.best_estimator_
        all_tuning.append({
            'antibiotic': antibiotic,
            'model': model_name,
            'best_inner_cv_balanced_accuracy': search.best_score_,
            'best_params': json.dumps(search.best_params_, ensure_ascii=False)
        })

        pred = best.predict(X_test)
        proba = best.predict_proba(X_test)[:, 1]
        cm = confusion_matrix(y_test, pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        all_holdout.append({
            'antibiotic': antibiotic, 'model': model_name, 'n_train': len(X_train), 'n_test': len(X_test), 'n_features': len(feature_cols),
            'accuracy': accuracy_score(y_test, pred),
            'balanced_accuracy': balanced_accuracy_score(y_test, pred),
            'f1_resistant': f1_score(y_test, pred, pos_label=1),
            'roc_auc': roc_auc_score(y_test, proba),
            'recall_resistant': recall_score(y_test, pred, pos_label=1),
            'precision_resistant': precision_score(y_test, pred, pos_label=1, zero_division=0),
            'specificity_susceptible': specificity_score(y_test, pred),
            'true_susceptible': int(tn), 'false_resistant': int(fp),
            'false_susceptible': int(fn), 'true_resistant': int(tp),
            'confusion_matrix': json.dumps(cm.tolist())
        })
        plot_confusion(antibiotic, model_name, cm)
        plot_roc(antibiotic, model_name, y_test, proba)

        cv_eval = StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE)
        cvres = cross_validate(best, X, y, scoring=scoring, cv=cv_eval, n_jobs=-1)
        cv_row = {'antibiotic': antibiotic, 'model': model_name, 'n_samples': len(df), 'n_features': len(feature_cols)}
        for metric in scoring:
            vals = cvres[f'test_{metric}']
            cv_row[f'{metric}_mean'] = np.nanmean(vals)
            cv_row[f'{metric}_std'] = np.nanstd(vals, ddof=1)
        all_cv.append(cv_row)

        final = clone(best)
        final.fit(X, y)
        imp = extract_importance(final, feature_cols, model_name)
        imp['antibiotic'] = antibiotic
        imp['model'] = model_name
        all_imp.append(imp)

        raw, stability = gene_stability(best, model_name, X, y, feature_cols, antibiotic, top_n=20)
        all_stab_raw.append(raw)
        all_stab.append(stability)
        plot_stable_genes(antibiotic, model_name, stability, n=20)

# Guardar tablas
frames = {
    'data_audit_summary.csv': pd.DataFrame(all_audit),
    'hyperparameter_tuning_summary.csv': pd.DataFrame(all_tuning),
    'holdout_test_metrics_tuned_models.csv': pd.DataFrame(all_holdout),
    'cross_validation_metrics_tuned_models.csv': pd.DataFrame(all_cv),
    'final_feature_importances_all_models.csv': pd.concat(all_imp, ignore_index=True),
    'gene_importance_stability_top20_by_fold.csv': pd.concat(all_stab, ignore_index=True),
    'gene_importance_stability_raw_fold_tops.csv': pd.concat(all_stab_raw, ignore_index=True)
}
for name, table in frames.items():
    out = table.copy()
    for col in out.columns:
        if out[col].dtype.kind in 'fc':
            out[col] = out[col].round(4)
    out.to_csv(OUT / name, index=False)

cvsum = frames['cross_validation_metrics_tuned_models.csv']
holdout = frames['holdout_test_metrics_tuned_models.csv']
stability = frames['gene_importance_stability_top20_by_fold.csv']

selected = []
for antibiotic in DATASETS:
    sub = cvsum[cvsum.antibiotic == antibiotic]
    idx = sub['balanced_accuracy_mean'].idxmax()
    row = sub.loc[idx].to_dict()
    h = holdout[(holdout.antibiotic == antibiotic) & (holdout.model == row['model'])].iloc[0].to_dict()
    row.update({
        'selection_rule': 'highest 3-fold CV balanced accuracy',
        'holdout_balanced_accuracy': h['balanced_accuracy'],
        'holdout_roc_auc': h['roc_auc'],
        'holdout_false_susceptible': h['false_susceptible'],
        'holdout_false_resistant': h['false_resistant']
    })
    selected.append(row)
selected = pd.DataFrame(selected)
for col in selected.columns:
    if selected[col].dtype.kind in 'fc': selected[col] = selected[col].round(4)
selected.to_csv(OUT / 'selected_models_summary.csv', index=False)

stable_selected = []
for _, r in selected.iterrows():
    stable_selected.append(stability[(stability.antibiotic == r.antibiotic) & (stability.model == r.model)].head(25))
stable_selected = pd.concat(stable_selected, ignore_index=True)
stable_selected.to_csv(OUT / 'stable_top_genes_selected_models.csv', index=False)

for antibiotic in DATASETS:
    plot_cv(antibiotic, cvsum)

# Reporte en Markdown
metric_cols = ['antibiotic','model','balanced_accuracy_mean','balanced_accuracy_std','f1_resistant_mean','f1_resistant_std','recall_resistant_mean','recall_resistant_std','roc_auc_mean','roc_auc_std','specificity_susceptible_mean','specificity_susceptible_std']
hold_cols = ['antibiotic','model','accuracy','balanced_accuracy','f1_resistant','roc_auc','recall_resistant','specificity_susceptible','false_susceptible','false_resistant','true_resistant','true_susceptible']
sections = []
for _, r in selected.iterrows():
    genes = stable_selected[(stable_selected.antibiotic == r.antibiotic) & (stable_selected.model == r.model)].head(12)
    sections.append(f"### {r.antibiotic} — {r.model}\n\n" + genes[['gene','selection_frequency','mean_rank','mean_abs_importance']].to_markdown(index=False))

report = f"""
# Resultados finales robustecidos del modelo AMR en *Escherichia coli*

## Auditoría de datos

{frames['data_audit_summary.csv'].to_markdown(index=False)}

## Fortalecimiento del análisis

Se añadieron tres pasos para fortalecer el análisis: ajuste pequeño de hiperparámetros por validación cruzada interna, evaluación por validación cruzada estratificada de tres particiones y análisis de estabilidad de genes importantes. La métrica principal para seleccionar modelos fue balanced accuracy, porque las clases resistente y susceptible están desbalanceadas. También se reportó recall de resistentes, porque un falso susceptible es el error más delicado en resistencia antimicrobiana.

## Validación cruzada de modelos ajustados

{cvsum[metric_cols].to_markdown(index=False)}

## Evaluación en conjunto de prueba retenido

{holdout[hold_cols].to_markdown(index=False)}

## Modelos seleccionados

{selected[['antibiotic','model','balanced_accuracy_mean','balanced_accuracy_std','recall_resistant_mean','roc_auc_mean','holdout_false_susceptible','holdout_false_resistant','selection_rule']].to_markdown(index=False)}

## Genes más estables en los modelos seleccionados

{chr(10).join(sections)}

## Interpretación lista para informe

Los resultados muestran que la presencia/ausencia de genes AMR permite clasificar genomas de *Escherichia coli* como resistentes o susceptibles frente a ciprofloxacin y cefotaxime. La validación cruzada estratificada confirma que el desempeño no depende únicamente de una partición entrenamiento-prueba.

Para ciprofloxacin, los genes importantes se relacionan con mecanismos conocidos de resistencia a fluoroquinolonas, incluyendo genes asociados con topoisomerasas y determinantes plasmídicos. Para cefotaxime, los genes más relevantes pertenecen principalmente a familias de beta-lactamasas, especialmente genes tipo CTX-M y otros determinantes relacionados con resistencia a cefalosporinas.

Estos resultados pueden presentarse como resultados principales del modelo base robustecido. La principal limitación restante es que el análisis usa genes AMR anotados y todavía no incorpora estructura poblacional, linajes bacterianos, mutaciones puntuales no anotadas ni validación externa en una cohorte independiente.
""".strip()
(OUT / 'REPORTE_RESULTADOS_ROBUSTOS.md').write_text(report, encoding='utf-8')

with zipfile.ZipFile('ecoli_amr_robust_final_analysis.zip', 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for f in OUT.rglob('*'):
        z.write(f, f.relative_to(Path('.')))

print('\nAnálisis robustecido completado.')
print('Revisa la carpeta:', OUT)
print('ZIP creado: ecoli_amr_robust_final_analysis.zip')
