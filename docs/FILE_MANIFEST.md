# Manifiesto de archivos

## Análisis activos

- `scripts/prepare_model_data.py`: reconstruye matrices y datasets model-ready.
- `scripts/run_optimized_models.py`: ejecuta regresión logística, random forest y XGBoost.
- `scripts/run_mlp_multiseed.py`: evalúa el MLP como cuarto modelo.
- `scripts/run_robustness_by_lab_method.py`: ejecuta robustez leave-method-out.
- `scripts/run_all_analyses.py`: ejecuta los tres bloques de modelado.

## Resultados finales

- `results/optimized/`: resultados finales de los tres modelos principales.
- `results/mlp/`: resultados del MLP multisemilla.
- `results/robustness_method/`: robustez por método de laboratorio.
- `results/tablas_finales/`: tablas utilizadas en el informe.
- `figures/`: figuras finales numeradas del 1 al 7.
- `report/informe_final.pdf`: informe final entregado.

## Material histórico

- `scripts/archive/`: scripts originales y versiones previas.
- `results/archive/`: corridas anteriores conservadas para trazabilidad.
- `docs/preparation/`: preparación y planteamiento.
- `docs/review/`: documentos de revisión de resultados.
