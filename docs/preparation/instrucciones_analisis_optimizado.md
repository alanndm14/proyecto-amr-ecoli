# Script final optimizado para AMR en E. coli

## Archivos necesarios en la misma carpeta

- `dataset_ciprofloxacin_model_ready.csv`
- `dataset_cefotaxime_model_ready.csv`
- `generate_final_optimized_amr_analysis.py`

## Instalar dependencias

En Jupyter:

```python
%pip install pandas numpy matplotlib scikit-learn xgboost tabulate joblib
```

En terminal:

```bash
pip install pandas numpy matplotlib scikit-learn xgboost tabulate joblib
```

## Ejecutar

En Jupyter:

```python
!python generate_final_optimized_amr_analysis.py
```

En terminal:

```bash
python generate_final_optimized_amr_analysis.py
```

## Parámetros por defecto usados

- `RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=42)`
- `RandomizedSearchCV(..., refit="balanced_accuracy", n_jobs=4, pre_dispatch="n_jobs")`
- `n_iter_logreg = 40`
- `n_iter_rf = 50`
- `n_iter_xgb = 60`
- modelos internos con `n_jobs=1` para evitar paralelización anidada

## Salida

El script genera una carpeta:

```text
ecoli_amr_resultados_finales_optimizados
```

y un ZIP:

```text
ecoli_amr_resultados_finales_optimizados.zip
```

Archivos clave:

- `REPORTE_FINAL_OPTIMIZADO.md`
- `modelos_seleccionados_resumen.csv`
- `comparacion_modelos_validacion_cruzada.csv`
- `metricas_holdout_todos_los_modelos.csv`
- `top_genes_estables_modelos_seleccionados.csv`
- figuras PNG en español

## Si tarda demasiado

Puedes hacer una corrida más ligera así:

```bash
python generate_final_optimized_amr_analysis.py --n_iter_logreg 20 --n_iter_rf 25 --n_iter_xgb 30 --n_boot 250
```

Pero para el análisis final usa los valores por defecto.
