from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
    recall_score, precision_score, confusion_matrix
)

RANDOM_STATE = 42
TEST_SIZE = 0.20
MIN_N_GROUP = 100
BASE = Path('/mnt/data/proyecto_genomica_extracted/proyecto genomica')
OUTDIR = Path('/mnt/data/robustez_metodo_laboratorio')
OUTDIR.mkdir(exist_ok=True)

DATASETS = {
    'ciprofloxacin': BASE / 'dataset_ciprofloxacin_model_ready.csv',
    'cefotaxime': BASE / 'dataset_cefotaxime_model_ready.csv',
}


def make_model():
    return make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(
            max_iter=1000,
            class_weight='balanced',
            random_state=RANDOM_STATE,
            solver='liblinear',
            penalty='l1',
            C=0.1,
        ),
    )


def normalize_method(x):
    if pd.isna(x) or str(x).strip() == '':
        return 'No especificado'
    return str(x).strip()


def compute_metrics(y_true, y_pred, y_score=None):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    out = {
        'n': int(len(y_true)),
        'n_susceptible': int((y_true == 0).sum()),
        'n_resistant': int((y_true == 1).sum()),
        'accuracy': accuracy_score(y_true, y_pred),
        'balanced_accuracy': balanced_accuracy_score(y_true, y_pred),
        'recall_resistant': recall_score(y_true, y_pred, zero_division=0),
        'specificity_susceptible': tn / (tn + fp) if (tn + fp) else np.nan,
        'precision_resistant': precision_score(y_true, y_pred, zero_division=0),
        'f1_resistant': f1_score(y_true, y_pred, zero_division=0),
        'true_susceptible': int(tn),
        'false_resistant': int(fp),
        'false_susceptible': int(fn),
        'true_resistant': int(tp),
    }
    if y_score is not None and len(np.unique(y_true)) == 2:
        out['roc_auc'] = roc_auc_score(y_true, y_score)
    else:
        out['roc_auc'] = np.nan
    return out


records_holdout = []
records_leave_method_out = []

for antibiotic, path in DATASETS.items():
    df = pd.read_csv(path, low_memory=False)
    df['method_clean'] = df['laboratory_typing_method'].apply(normalize_method)
    feature_cols = [c for c in df.columns if c.startswith('amr__')]
    X = df[feature_cols].astype(np.int8)
    y = df['y'].astype(int)
    idx = np.arange(len(df))

    # A) Misma partición del análisis principal: entrenar en 80%, evaluar el 20% por subgrupo de método.
    train_idx, test_idx = train_test_split(idx, test_size=TEST_SIZE, stratify=y, random_state=RANDOM_STATE)
    model = make_model()
    model.fit(X.iloc[train_idx], y.iloc[train_idx])
    y_pred = model.predict(X.iloc[test_idx])
    y_score = model.predict_proba(X.iloc[test_idx])[:, 1]
    test_df = df.iloc[test_idx].copy()
    test_df['y_pred'] = y_pred
    test_df['y_score'] = y_score

    overall = compute_metrics(test_df['y'], test_df['y_pred'], test_df['y_score'])
    overall.update({'antibiotic': antibiotic, 'evaluation': 'holdout_subgroup', 'laboratory_typing_method': 'TOTAL_HOLDOUT'})
    records_holdout.append(overall)

    for method, g in test_df.groupby('method_clean'):
        if len(g) < MIN_N_GROUP or g['y'].nunique() < 2:
            continue
        m = compute_metrics(g['y'], g['y_pred'], g['y_score'])
        m.update({'antibiotic': antibiotic, 'evaluation': 'holdout_subgroup', 'laboratory_typing_method': method})
        records_holdout.append(m)

    # B) Validación aproximada externa: entrenar sin un método y probar solo en ese método.
    method_counts = df.groupby('method_clean')['y'].agg(['size', 'sum'])
    candidate_methods = []
    for method, row in method_counts.iterrows():
        n = int(row['size'])
        n_res = int(row['sum'])
        n_sus = n - n_res
        if n >= 500 and n_res >= 20 and n_sus >= 20 and method != 'No especificado':
            candidate_methods.append(method)

    for method in candidate_methods:
        test_mask = df['method_clean'].eq(method)
        train_mask = ~test_mask
        # entrenamiento debe tener ambas clases
        if y[train_mask].nunique() < 2 or y[test_mask].nunique() < 2:
            continue
        model = make_model()
        model.fit(X.loc[train_mask], y.loc[train_mask])
        pred = model.predict(X.loc[test_mask])
        score = model.predict_proba(X.loc[test_mask])[:, 1]
        m = compute_metrics(y.loc[test_mask], pred, score)
        m.update({'antibiotic': antibiotic, 'evaluation': 'leave_method_out', 'laboratory_typing_method': method})
        records_leave_method_out.append(m)

holdout_df = pd.DataFrame(records_holdout)
lomo_df = pd.DataFrame(records_leave_method_out)
all_df = pd.concat([holdout_df, lomo_df], ignore_index=True)

# Orden de columnas
cols = ['evaluation', 'antibiotic', 'laboratory_typing_method', 'n', 'n_susceptible', 'n_resistant',
        'accuracy', 'balanced_accuracy', 'roc_auc', 'recall_resistant', 'specificity_susceptible',
        'precision_resistant', 'f1_resistant', 'false_susceptible', 'false_resistant',
        'true_susceptible', 'true_resistant']
all_df = all_df[cols]
holdout_df = holdout_df[cols]
lomo_df = lomo_df[cols]

all_df.to_csv(OUTDIR / 'robustez_metodo_laboratorio_todos.csv', index=False)
holdout_df.to_csv(OUTDIR / 'robustez_holdout_por_metodo_laboratorio.csv', index=False)
lomo_df.to_csv(OUTDIR / 'robustez_leave_method_out.csv', index=False)

# Figuras simples: balanced accuracy leave-method-out
for evaluation, data, filename, title_prefix in [
    ('holdout_subgroup', holdout_df[holdout_df['laboratory_typing_method'] != 'TOTAL_HOLDOUT'], 'figura_robustez_holdout_por_metodo.png', 'Holdout por método'),
    ('leave_method_out', lomo_df, 'figura_robustez_leave_method_out.png', 'Leave-method-out'),
]:
    if data.empty:
        continue
    plot_df = data.copy().sort_values(['antibiotic', 'balanced_accuracy'])
    labels = plot_df['antibiotic'] + ' — ' + plot_df['laboratory_typing_method']
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.45 * len(plot_df) + 1.5)))
    ax.barh(labels, plot_df['balanced_accuracy'])
    ax.set_xlim(0, 1)
    ax.set_xlabel('Balanced accuracy')
    ax.set_title(f'{title_prefix}: desempeño de regresión logística balanceada')
    for i, val in enumerate(plot_df['balanced_accuracy']):
        ax.text(val + 0.01, i, f'{val:.3f}', va='center', fontsize=9)
    fig.tight_layout()
    fig.savefig(OUTDIR / filename, dpi=250)
    plt.close(fig)

# Reporte markdown breve
md = []
md.append('# Análisis de robustez por método de laboratorio\n')
md.append('Modelo usado: regresión logística balanceada con regularización L1 y `C=0.1`. Para esta comparación rápida se usó `solver=liblinear`, que reproduce casi exactamente el desempeño holdout del modelo seleccionado (`saga`, L1, C=0.1) y permite hacer las particiones por método sin esperar una nueva búsqueda de hiperparámetros.\n')
md.append('Se hicieron dos comparaciones: (1) desempeño por método dentro del holdout original 80/20, y (2) validación leave-method-out, donde se entrena sin un método de laboratorio y se prueba únicamente en ese método.\n')

for antibiotic in DATASETS:
    md.append(f'## {antibiotic}\n')
    md.append('### Holdout por método\n')
    t = holdout_df[holdout_df['antibiotic'].eq(antibiotic)].copy()
    md.append(t[['laboratory_typing_method','n','n_resistant','balanced_accuracy','roc_auc','recall_resistant','specificity_susceptible','false_susceptible','false_resistant']].round(3).to_markdown(index=False))
    md.append('\n\n### Leave-method-out\n')
    t2 = lomo_df[lomo_df['antibiotic'].eq(antibiotic)].copy()
    if t2.empty:
        md.append('No hubo métodos con tamaño suficiente para esta comparación.\n')
    else:
        md.append(t2[['laboratory_typing_method','n','n_resistant','balanced_accuracy','roc_auc','recall_resistant','specificity_susceptible','false_susceptible','false_resistant']].round(3).to_markdown(index=False))
    md.append('\n')

Path(OUTDIR / 'reporte_robustez_metodo_laboratorio.md').write_text('\n'.join(md), encoding='utf-8')

print('Archivos generados en:', OUTDIR)
print(all_df.round(3).to_string(index=False))
