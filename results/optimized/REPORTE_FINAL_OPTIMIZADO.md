# Resultados finales optimizados del análisis AMR en *Escherichia coli*

## Enfoque metodológico

Se compararon tres modelos de clasificación supervisada: regresión logística balanceada, bosque aleatorio balanceado y XGBoost ponderado. Para cada antibiótico se separó un conjunto retenido del 20% de los genomas, estratificado por fenotipo. La selección de hiperparámetros se realizó únicamente sobre el conjunto de entrenamiento mediante `RandomizedSearchCV` con validación cruzada repetida y estratificada de 5 particiones × 2 repeticiones.

La métrica principal de selección fue la exactitud balanceada, porque los datos presentan más genomas susceptibles que resistentes. También se guardaron ROC-AUC, F1 y sensibilidad para resistentes, ya que un falso susceptible es el error más delicado en resistencia antimicrobiana.

## Auditoría de datos

| antibiotic    |   n_samples |   n_features |   n_susceptible |   n_resistant |   resistant_fraction |   n_train |   n_test |   train_susceptible |   train_resistant |   scale_pos_weight_xgboost_train |   duplicated_genome_ids |
|:--------------|------------:|-------------:|----------------:|--------------:|---------------------:|----------:|---------:|--------------------:|------------------:|---------------------------------:|------------------------:|
| ciprofloxacin |        7313 |          633 |            5462 |          1851 |             0.253111 |      5850 |     1463 |                4369 |              1481 |                          2.95003 |                       0 |
| cefotaxime    |        6590 |          560 |            5318 |          1272 |             0.19302  |      5272 |     1318 |                4254 |              1018 |                          4.17878 |                       0 |

## Comparación de modelos en validación cruzada

| antibiotic    | model_label                    |   n_iter |   fits_requested |   cv_balanced_accuracy_mean |   cv_balanced_accuracy_std |   cv_f1_resistant_mean |   cv_roc_auc_mean |   cv_recall_resistant_mean | best_params                                                                                                                                                                        |
|:--------------|:-------------------------------|---------:|-----------------:|----------------------------:|---------------------------:|-----------------------:|------------------:|---------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ciprofloxacin | Regresión logística balanceada |        5 |               50 |                    0.845503 |                   0.013063 |               0.765473 |          0.910792 |                   0.776494 | {"logisticregression__solver": "saga", "logisticregression__penalty": "l1", "logisticregression__C": 0.1}                                                                          |
| ciprofloxacin | Bosque aleatorio balanceado    |        5 |               50 |                    0.845084 |                   0.009344 |               0.774077 |          0.928392 |                   0.758261 | {"n_estimators": 500, "min_samples_split": 6, "min_samples_leaf": 1, "max_samples": null, "max_features": "sqrt", "max_depth": 50, "class_weight": "balanced", "bootstrap": true}  |
| ciprofloxacin | XGBoost ponderado              |        5 |               50 |                    0.839836 |                   0.010376 |               0.766324 |          0.911261 |                   0.750167 | {"subsample": 0.85, "reg_lambda": 2, "reg_alpha": 1, "n_estimators": 800, "min_child_weight": 2, "max_depth": 4, "learning_rate": 0.05, "gamma": 0.1, "colsample_bytree": 0.85}    |
| cefotaxime    | Regresión logística balanceada |        5 |               50 |                    0.825456 |                   0.013298 |               0.70849  |          0.887279 |                   0.729895 | {"logisticregression__solver": "saga", "logisticregression__penalty": "l1", "logisticregression__C": 0.1}                                                                          |
| cefotaxime    | Bosque aleatorio balanceado    |        5 |               50 |                    0.794319 |                   0.017673 |               0.641114 |          0.874023 |                   0.708287 | {"n_estimators": 800, "min_samples_split": 2, "min_samples_leaf": 8, "max_samples": 0.85, "max_features": "sqrt", "max_depth": 50, "class_weight": "balanced", "bootstrap": true}  |
| cefotaxime    | XGBoost ponderado              |        5 |               50 |                    0.802768 |                   0.020262 |               0.67781  |          0.885233 |                   0.686166 | {"subsample": 1.0, "reg_lambda": 0.5, "reg_alpha": 1, "n_estimators": 1200, "min_child_weight": 10, "max_depth": 4, "learning_rate": 0.01, "gamma": 0.5, "colsample_bytree": 0.85} |

## Evaluación en conjunto retenido

| antibiotic    | model_label                    | selected_for_final_report   |   accuracy |   balanced_accuracy |   balanced_accuracy_ci95_low |   balanced_accuracy_ci95_high |   f1_resistant |   roc_auc |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |
|:--------------|:-------------------------------|:----------------------------|-----------:|--------------------:|-----------------------------:|------------------------------:|---------------:|----------:|-------------------:|--------------------------:|--------------------:|------------------:|
| ciprofloxacin | Regresión logística balanceada | True                        |   0.874231 |            0.844316 |                     0.82258  |                      0.867551 |       0.759162 |  0.90585  |           0.783784 |                  0.904849 |                  80 |               104 |
| ciprofloxacin | Bosque aleatorio balanceado    | False                       |   0.872864 |            0.837144 |                     0.815206 |                      0.85738  |       0.75266  |  0.922398 |           0.764865 |                  0.909424 |                  87 |                99 |
| ciprofloxacin | XGBoost ponderado              | False                       |   0.874231 |            0.836271 |                     0.813632 |                      0.857633 |       0.753351 |  0.907045 |           0.759459 |                  0.913083 |                  89 |                95 |
| cefotaxime    | Regresión logística balanceada | True                        |   0.890744 |            0.82743  |                     0.792302 |                      0.859126 |       0.71875  |  0.893371 |           0.724409 |                  0.930451 |                  70 |                74 |
| cefotaxime    | Bosque aleatorio balanceado    | False                       |   0.830046 |            0.791335 |                     0.760579 |                      0.819637 |       0.622896 |  0.878893 |           0.728346 |                  0.854323 |                  69 |               155 |
| cefotaxime    | XGBoost ponderado              | False                       |   0.879363 |            0.809891 |                     0.773682 |                      0.837482 |       0.690058 |  0.886171 |           0.69685  |                  0.922932 |                  77 |                82 |

## Modelos seleccionados

| antibiotic    | model_label                    |   cv_balanced_accuracy_mean |   cv_balanced_accuracy_std |   holdout_balanced_accuracy |   holdout_roc_auc |   holdout_recall_resistant |   holdout_specificity_susceptible |   holdout_false_susceptible |   holdout_false_resistant | selection_rule                                                  |
|:--------------|:-------------------------------|----------------------------:|---------------------------:|----------------------------:|------------------:|---------------------------:|----------------------------------:|----------------------------:|--------------------------:|:----------------------------------------------------------------|
| ciprofloxacin | Regresión logística balanceada |                    0.845503 |                   0.013063 |                    0.844316 |          0.90585  |                   0.783784 |                          0.904849 |                          80 |                       104 | Mayor exactitud balanceada media en validación cruzada repetida |
| cefotaxime    | Regresión logística balanceada |                    0.825456 |                   0.013298 |                    0.82743  |          0.893371 |                   0.724409 |                          0.930451 |                          70 |                        74 | Mayor exactitud balanceada media en validación cruzada repetida |

## Interpretación por antibiótico

### ciprofloxacin — Regresión logística balanceada

El modelo seleccionado fue **Regresión logística balanceada**, con exactitud balanceada media de 0.8455 en validación cruzada repetida (DE = 0.013). En el conjunto retenido, obtuvo exactitud balanceada de 0.8443 (IC95%: 0.8226–0.8676), ROC-AUC de 0.9059, sensibilidad para resistentes de 0.7838 y especificidad para susceptibles de 0.9048. Los falsos susceptibles fueron 80.

Genes más estables del modelo seleccionado:

| gene                                                            |   selection_frequency |   selection_fraction |   mean_rank |   mean_abs_importance |
|:----------------------------------------------------------------|----------------------:|---------------------:|------------:|----------------------:|
| escherichia_quinolone_resistant_parc                            |                    15 |             1        |     1.2     |              0.812164 |
| escherichia_coli_parc_conferring_resistance_to_fluoroquinolones |                    15 |             1        |     1.86667 |              0.784133 |
| escherichia_quinolone_resistant_gyra                            |                    15 |             1        |     2.93333 |              0.691956 |
| escherichia_coli_gyra_conferring_resistance_to_fluoroquinolones |                    15 |             1        |     4       |              0.449744 |
| catb3                                                           |                    15 |             1        |     5.33333 |              0.350012 |
| escherichia_quinolone_resistant_pare                            |                    15 |             1        |     5.86667 |              0.317047 |
| extended-spectrum_class_a_beta-lactamase_ctx-m-15               |                    15 |             1        |     9.13333 |              0.254933 |
| blaec                                                           |                    15 |             1        |     9.4     |              0.254733 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-15  |                    14 |             0.933333 |    10.6429  |              0.234347 |
| ctx-m-15                                                        |                    14 |             0.933333 |    10.7143  |              0.231159 |
| escherichia_fosfomycin_resistant_ptsi                           |                    13 |             0.866667 |    12.1538  |              0.21353  |
| vgac                                                            |                    13 |             0.866667 |    12.7692  |              0.213028 |

### cefotaxime — Regresión logística balanceada

El modelo seleccionado fue **Regresión logística balanceada**, con exactitud balanceada media de 0.8255 en validación cruzada repetida (DE = 0.013). En el conjunto retenido, obtuvo exactitud balanceada de 0.8274 (IC95%: 0.7923–0.8591), ROC-AUC de 0.8934, sensibilidad para resistentes de 0.7244 y especificidad para susceptibles de 0.9305. Los falsos susceptibles fueron 70.

Genes más estables del modelo seleccionado:

| gene                                                             |   selection_frequency |   selection_fraction |   mean_rank |   mean_abs_importance |
|:-----------------------------------------------------------------|----------------------:|---------------------:|------------:|----------------------:|
| ctx-m-15                                                         |                    15 |                    1 |     1.4     |              0.652407 |
| ctx-m_family                                                     |                    15 |                    1 |     1.73333 |              0.644212 |
| cmy-2                                                            |                    15 |                    1 |     3.46667 |              0.525094 |
| ctx-m-14                                                         |                    15 |                    1 |     4.66667 |              0.512261 |
| extended-spectrum_class_a_beta-lactamase_ctx-m-15                |                    15 |                    1 |     4.8     |              0.490932 |
| multispecies_class_c_beta-lactamase_cmy-2                        |                    15 |                    1 |     6.13333 |              0.461246 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-15   |                    15 |                    1 |     6.6     |              0.454099 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-14   |                    15 |                    1 |     8.93333 |              0.392454 |
| class_a_beta-lactamase_ec_3.5.2.6_ctx-m_family_extended-spectrum |                    15 |                    1 |     9       |              0.389272 |
| extended-spectrum_class_c_beta-lactamase_cmy-2                   |                    15 |                    1 |     9.13333 |              0.39141  |
| cmycmy-2cfelat_family                                            |                    15 |                    1 |    10.2     |              0.370609 |
| extended-spectrum_class_a_beta-lactamase_ctx-m-14                |                    15 |                    1 |    13.6667  |              0.249763 |


## Conclusión general

La comparación optimizada muestra si la señal de genes AMR permite clasificar de manera consistente genomas resistentes y susceptibles a ciprofloxacin y cefotaxime. El uso de búsqueda aleatoria controlada, validación cruzada repetida, conjunto retenido independiente, intervalos de confianza bootstrap y estabilidad de genes permite presentar estos resultados como análisis final del modelo base optimizado, no como una corrida preliminar.

## Limitaciones

Aunque el análisis es más robusto, todavía debe mencionarse que no incluye una cohorte externa completamente independiente, no incorpora explícitamente mutaciones puntuales no anotadas y no modela estructura poblacional o linajes bacterianos.