# Notas de reproducibilidad

## Entorno versionado

`requirements.txt` fija las versiones instaladas en el entorno local usado para
validar este repositorio.

## Validaciones realizadas

- Los scripts activos compilan correctamente.
- `scripts/prepare_model_data.py` reconstruye exactamente las columnas finales:
  633 genes AMR para ciprofloxacin y 560 para cefotaxime.
- `scripts/run_robustness_by_lab_method.py` reproduce exactamente las filas y
  métricas leave-method-out del resultado final.
- `scripts/run_mlp_multiseed.py` reproduce exactamente cinco de los seis
  ensembles finales en el entorno validado.

## Variación del MLP

El ensemble MLP-32 para ciprofloxacin puede variar ligeramente entre entornos
aunque las semillas sean fijas. Esto ocurre porque el entrenamiento neuronal
puede depender de detalles numéricos de BLAS, CPU y versiones de
scikit-learn/numpy.

Los resultados originales usados en el informe se preservan en `results/mlp/`.
Las conclusiones no cambian: el MLP permanece como cuarto modelo comparativo y
la regresión logística balanceada sigue siendo el modelo seleccionado.

## Corrida optimizada

La corrida principal reportada se generó con:

```bash
python scripts/run_optimized_models.py --n_iter_logreg 5 --n_iter_rf 5 --n_iter_xgb 5 --n_boot 100 --n_jobs_search 1 --verbose 1
```

El uso de `--n_jobs_search 1` evita errores de multiprocessing en algunos
equipos Windows.
