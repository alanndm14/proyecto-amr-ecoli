# Predicción de resistencia antimicrobiana en *Escherichia coli*

Repositorio reproducible del proyecto de predicción de resistencia a
ciprofloxacin y cefotaxime mediante genes AMR codificados como presencia/ausencia.

## Objetivo

Comparar modelos supervisados capaces de predecir si un genoma de
*Escherichia coli* es resistente o susceptible a cada antibiótico. Se evaluaron:

1. Regresión logística balanceada.
2. Random forest balanceado.
3. XGBoost ponderado.
4. Perceptrón multicapa (MLP) ponderado y evaluado con siete semillas.

La métrica principal fue `balanced_accuracy`. También se revisaron ROC-AUC,
F1 y recall para resistentes, especificidad y falsos susceptibles.

## Origen de los datos

Los genomas, fenotipos de susceptibilidad antimicrobiana y anotaciones de genes
AMR provienen de [BV-BRC](https://www.bv-brc.org/). Se conservaron registros de
*E. coli* con fenotipo de laboratorio Resistant/Susceptible para ciprofloxacin
y cefotaxime.

El archivo crudo de Specialty Genes supera 100 MB y no se versiona en Git
normal. Consulte `data/raw/README.md` y `data/external_links.md`.

## Estructura

```text
proyecto-amr-ecoli/
|-- README.md
|-- requirements.txt
|-- .gitignore
|-- data/
|   |-- raw/
|   |-- processed/
|   `-- external_links.md
|-- scripts/
|-- results/
|-- figures/
|-- report/
`-- docs/
```

- `data/processed/`: datos limpios, matrices de presencia/ausencia y datasets finales.
- `scripts/`: scripts activos reproducibles; `scripts/archive/` conserva originales.
- `results/optimized/`: regresión logística, random forest y XGBoost.
- `results/mlp/`: MLP multisemilla integrado como cuarto modelo.
- `results/robustness_method/`: validación por método y leave-method-out.
- `results/tablas_finales/`: tablas usadas en el informe.
- `figures/`: figuras finales numeradas del 1 al 7.
- `report/`: informe final, borrador y material listo para insertar.
- `docs/`: preparación, revisión y manifiesto de archivos.

## Instalación

Se recomienda Python 3.10 o superior:

```bash
python -m venv .venv
```

En Windows:

```bash
.venv\Scripts\activate
python -m pip install -r requirements.txt
```

En Linux/macOS:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Reproducción

### Ejecutar todos los análisis

Desde la raíz del repositorio:

```bash
python scripts/run_all_analyses.py
```

La secuencia usa la corrida optimizada ligera que generó los resultados
reportados: 5 iteraciones por familia de modelos, 100 bootstraps y un solo
trabajo paralelo. La corrida completa puede tardar alrededor de 1-2 horas.

### Ejecutar bloques individuales

Reconstruir matrices y datasets model-ready desde los datos limpios:

```bash
python scripts/prepare_model_data.py
```

Modelos optimizados principales:

```bash
python scripts/run_optimized_models.py --n_iter_logreg 5 --n_iter_rf 5 --n_iter_xgb 5 --n_boot 100 --n_jobs_search 1 --verbose 1
```

MLP como cuarto modelo:

```bash
python scripts/run_mlp_multiseed.py
```

Robustez leave-method-out por método de laboratorio:

```bash
python scripts/run_robustness_by_lab_method.py
```

## Resultados finales

Los modelos seleccionados fueron regresiones logísticas balanceadas para ambos
antibióticos. En holdout alcanzaron:

| Antibiótico | Balanced accuracy | ROC-AUC |
|---|---:|---:|
| Ciprofloxacin | 0.844 | 0.906 |
| Cefotaxime | 0.827 | 0.893 |

El MLP se conserva como cuarto modelo comparativo. No superó la balanced
accuracy de la regresión logística, pero redujo ligeramente falsos susceptibles
en ambos antibióticos.

Archivos principales:

- `results/optimized/modelos_seleccionados_resumen.csv`
- `results/optimized/metricas_holdout_todos_los_modelos.csv`
- `results/mlp/comparacion_holdout_modelos_previos_mas_mlp.csv`
- `results/robustness_method/robustez_leave_method_out.csv`
- `results/tablas_finales/tabla_2_comparacion_modelos_incluye_mlp.csv`
- `report/informe_final.pdf`

## Figuras y tablas del informe

Las figuras finales son exactamente:

1. `figures/figura_1_flujo_general.png`
2. `figures/figura_2_comparacion_balanced_accuracy_modelos.png`
3. `figures/figura_3_falsos_susceptibles_modelos.png`
4. `figures/figura_4_matrices_confusion_modelos_seleccionados.png`
5. `figures/figura_5_genes_estables_ciprofloxacin.png`
6. `figures/figura_6_genes_estables_cefotaxime.png`
7. `figures/figura_7_robustez_leave_method_out.png`


Las tablas finales están en `results/tablas_finales/`, incluyendo composición de
datasets, comparación de los cuatro modelos, robustez leave-method-out, genes
estables y resumen del MLP.

## Reproducibilidad y trazabilidad

- Semilla principal: `42`.
- Holdout estratificado: 80/20.
- Modelos principales seleccionados mediante validación cruzada repetida.
- MLP evaluado con siete semillas y ensemble de probabilidades.
- Robustez evaluada dejando fuera métodos de laboratorio completos.
- Scripts y resultados anteriores se conservan en carpetas `archive`.

El MLP puede mostrar diferencias numéricas pequeñas entre equipos aún usando
semillas fijas, debido a operaciones de álgebra lineal y versiones de
bibliotecas. Los resultados originales usados en el informe se conservan en
`results/mlp/`; consulte `docs/REPRODUCIBILITY_NOTES.md`.

Consulte `docs/FILE_MANIFEST.md` para el inventario detallado.
Consulte `docs/CODIFICACION.md` para la política de UTF-8 y nombres técnicos.
