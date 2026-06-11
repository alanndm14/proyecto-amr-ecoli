# Análisis robustecido en VSCode/Jupyter

Coloca este script en la misma carpeta donde están:

- `dataset_ciprofloxacin_model_ready.csv`
- `dataset_cefotaxime_model_ready.csv`

## En Jupyter

Ejecuta estas celdas:

```python
%pip install pandas scikit-learn matplotlib tabulate
```

```python
!python generate_robust_final_analysis.py
```

## Qué genera

Una carpeta llamada:

```text
ecoli_amr_robust_final_analysis/
```

Y un ZIP llamado:

```text
ecoli_amr_robust_final_analysis.zip
```

Con:

- validación cruzada estratificada
- ajuste pequeño de hiperparámetros
- métricas en conjunto de prueba retenido
- matrices de confusión
- curvas ROC
- estabilidad de genes importantes
- reporte en Markdown listo para adaptar

## Nota

La corrida puede tardar varios minutos porque entrena modelos varias veces. Si tarda demasiado, cambia en el script `n_estimators=60` por `n_estimators=30` para random forest.
