# Resultados finales del análisis de resistencia antimicrobiana en *Escherichia coli*

## Enfoque metodológico final

Para evitar presentar resultados meramente preliminares, el análisis final comparó tres modelos: regresión logística balanceada, bosque aleatorio balanceado y XGBoost ponderado. Se separó un conjunto de prueba retenido e independiente del 20% de los genomas. La selección de hiperparámetros se realizó únicamente sobre el 80% restante mediante validación cruzada repetida y estratificada de 5 particiones con 2 repeticiones. El criterio principal para elegir modelo fue la exactitud balanceada, adecuada cuando las clases resistente y susceptible están desbalanceadas. Además, se evaluó la estabilidad de los genes importantes mediante 15 ajustes repetidos del modelo seleccionado.

## Auditoría de datos de entrada

| antibiotic    |   n_samples |   n_features |   n_susceptible |   n_resistant |   resistant_fraction |   scale_pos_weight_used_for_xgboost |   duplicated_genome_ids |
|:--------------|------------:|-------------:|----------------:|--------------:|---------------------:|------------------------------------:|------------------------:|
| ciprofloxacin |        7313 |          633 |            5462 |          1851 |             0.253111 |                             2.95084 |                       0 |
| cefotaxime    |        6590 |          560 |            5318 |          1272 |             0.19302  |                             4.18082 |                       0 |

## Comparación de modelos en validación cruzada interna

| antibiotic    | model_label                    |   cv_balanced_accuracy_mean |   cv_balanced_accuracy_std |   cv_f1_resistant_mean |   cv_roc_auc_mean |   cv_recall_resistant_mean | best_params                                                                                              |
|:--------------|:-------------------------------|----------------------------:|---------------------------:|-----------------------:|------------------:|---------------------------:|:---------------------------------------------------------------------------------------------------------|
| ciprofloxacin | Regresión logística balanceada |                    0.844994 |                  0.0129478 |               0.764675 |          0.91083  |                   0.775818 | {"logisticregression__C": 0.1, "logisticregression__penalty": "l1"}                                      |
| ciprofloxacin | Bosque aleatorio balanceado    |                    0.847151 |                  0.0122584 |               0.761572 |          0.921867 |                   0.791691 | {"max_depth": null, "max_features": "sqrt", "min_samples_leaf": 4}                                       |
| ciprofloxacin | XGBoost ponderado              |                    0.844394 |                  0.0126535 |               0.772208 |          0.916349 |                   0.75894  | {"colsample_bytree": 1.0, "learning_rate": 0.1, "max_depth": 3, "min_child_weight": 1, "subsample": 1.0} |
| cefotaxime    | Regresión logística balanceada |                    0.825103 |                  0.0128722 |               0.707499 |          0.887209 |                   0.729895 | {"logisticregression__C": 0.1, "logisticregression__penalty": "l1"}                                      |
| cefotaxime    | Bosque aleatorio balanceado    |                    0.805941 |                  0.0130149 |               0.674923 |          0.890987 |                   0.702386 | {"max_depth": 25, "max_features": "sqrt", "min_samples_leaf": 4}                                         |
| cefotaxime    | XGBoost ponderado              |                    0.815504 |                  0.013228  |               0.712122 |          0.893686 |                   0.690599 | {"colsample_bytree": 1.0, "learning_rate": 0.1, "max_depth": 5, "min_child_weight": 1, "subsample": 1.0} |

## Resultados en el conjunto retenido

| antibiotic    | model_label                    | selected_for_final_report   |   accuracy |   balanced_accuracy |   balanced_accuracy_ci95_low |   balanced_accuracy_ci95_high |   f1_resistant |   roc_auc |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |
|:--------------|:-------------------------------|:----------------------------|-----------:|--------------------:|-----------------------------:|------------------------------:|---------------:|----------:|-------------------:|--------------------------:|--------------------:|------------------:|
| ciprofloxacin | Regresión logística balanceada | False                       |   0.874231 |            0.84521  |                     0.823376 |                      0.866668 |       0.759791 |  0.905984 |           0.786486 |                  0.903934 |                  79 |               105 |
| ciprofloxacin | Bosque aleatorio balanceado    | True                        |   0.864662 |            0.847745 |                     0.82542  |                      0.867904 |       0.7525   |  0.918328 |           0.813514 |                  0.881976 |                  69 |               129 |
| ciprofloxacin | XGBoost ponderado              | False                       |   0.871497 |            0.835335 |                     0.811571 |                      0.858414 |       0.75     |  0.916374 |           0.762162 |                  0.908509 |                  88 |               100 |
| cefotaxime    | Regresión logística balanceada | True                        |   0.891502 |            0.829399 |                     0.798663 |                      0.857511 |       0.721248 |  0.893964 |           0.728346 |                  0.930451 |                  69 |                74 |
| cefotaxime    | Bosque aleatorio balanceado    | False                       |   0.858118 |            0.807224 |                     0.775099 |                      0.835417 |       0.663063 |  0.894097 |           0.724409 |                  0.890038 |                  70 |               117 |
| cefotaxime    | XGBoost ponderado              | False                       |   0.897572 |            0.818172 |                     0.788485 |                      0.845261 |       0.721649 |  0.892809 |           0.688976 |                  0.947368 |                  79 |                56 |

## Modelos seleccionados para el informe final

| antibiotic    | model_label                    |   cv_balanced_accuracy_mean |   cv_balanced_accuracy_std |   holdout_balanced_accuracy |   holdout_roc_auc |   holdout_recall_resistant |   holdout_specificity_susceptible |   holdout_false_susceptible |   holdout_false_resistant | selection_rule                                                                      |
|:--------------|:-------------------------------|----------------------------:|---------------------------:|----------------------------:|------------------:|---------------------------:|----------------------------------:|----------------------------:|--------------------------:|:------------------------------------------------------------------------------------|
| ciprofloxacin | Bosque aleatorio balanceado    |                    0.847151 |                  0.0122584 |                    0.847745 |          0.918328 |                   0.813514 |                          0.881976 |                          69 |                       129 | Mayor exactitud balanceada en validación cruzada interna (5 folds × 2 repeticiones) |
| cefotaxime    | Regresión logística balanceada |                    0.825103 |                  0.0128722 |                    0.829399 |          0.893964 |                   0.728346 |                          0.930451 |                          69 |                        74 | Mayor exactitud balanceada en validación cruzada interna (5 folds × 2 repeticiones) |

## Interpretación por antibiótico

### ciprofloxacin — Bosque aleatorio balanceado

Modelo seleccionado con exactitud balanceada media en validación cruzada de 0.8472 (DE = 0.012). En el conjunto retenido, obtuvo exactitud balanceada de 0.8477, AUC de 0.9183, sensibilidad para resistentes de 0.8135 y especificidad para susceptibles de 0.882. Los errores más delicados, resistentes predichos como susceptibles, fueron 69.

Genes más estables del modelo seleccionado

| gene                                                            |   selection_frequency |   selection_fraction |   mean_rank |   mean_abs_importance |
|:----------------------------------------------------------------|----------------------:|---------------------:|------------:|----------------------:|
| escherichia_quinolone_resistant_gyra                            |                    15 |                    1 |     1       |             0.0874309 |
| escherichia_quinolone_resistant_parc                            |                    15 |                    1 |     2       |             0.0616019 |
| escherichia_coli_parc_conferring_resistance_to_fluoroquinolones |                    15 |                    1 |     3.2     |             0.0438513 |
| escherichia_coli_gyra_conferring_resistance_to_fluoroquinolones |                    15 |                    1 |     3.8     |             0.0416508 |
| mpha                                                            |                    15 |                    1 |     5       |             0.0310173 |
| mrx                                                             |                    15 |                    1 |     6.06667 |             0.0263894 |
| escherichia_quinolone_resistant_pare                            |                    15 |                    1 |     7.13333 |             0.0214985 |
| ctx-m-15                                                        |                    15 |                    1 |     8.2     |             0.018956  |
| aada5                                                           |                    15 |                    1 |     9.66667 |             0.016835  |
| sul1                                                            |                    15 |                    1 |    10.8667  |             0.0149292 |
| dfra17                                                          |                    15 |                    1 |    11.2     |             0.0147123 |
| teta                                                            |                    15 |                    1 |    11.4667  |             0.0144395 |

### cefotaxime — Regresión logística balanceada

Modelo seleccionado con exactitud balanceada media en validación cruzada de 0.8251 (DE = 0.013). En el conjunto retenido, obtuvo exactitud balanceada de 0.8294, AUC de 0.894, sensibilidad para resistentes de 0.7283 y especificidad para susceptibles de 0.9305. Los errores más delicados, resistentes predichos como susceptibles, fueron 69.

Genes más estables del modelo seleccionado

| gene                                                             |   selection_frequency |   selection_fraction |   mean_rank |   mean_abs_importance |
|:-----------------------------------------------------------------|----------------------:|---------------------:|------------:|----------------------:|
| ctx-m-15                                                         |                    15 |                    1 |     1.4     |              0.657705 |
| ctx-m_family                                                     |                    15 |                    1 |     1.73333 |              0.64599  |
| cmy-2                                                            |                    15 |                    1 |     3.53333 |              0.525622 |
| ctx-m-14                                                         |                    15 |                    1 |     4.66667 |              0.516021 |
| extended-spectrum_class_a_beta-lactamase_ctx-m-15                |                    15 |                    1 |     4.73333 |              0.491617 |
| multispecies_class_c_beta-lactamase_cmy-2                        |                    15 |                    1 |     6.13333 |              0.461162 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-15   |                    15 |                    1 |     6.53333 |              0.454348 |
| class_a_beta-lactamase_ec_3.5.2.6_ctx-m_family_extended-spectrum |                    15 |                    1 |     8.93333 |              0.392413 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-14   |                    15 |                    1 |     9.06667 |              0.39131  |
| extended-spectrum_class_c_beta-lactamase_cmy-2                   |                    15 |                    1 |     9.13333 |              0.391706 |
| cmycmy-2cfelat_family                                            |                    15 |                    1 |    10.2     |              0.371806 |
| extended-spectrum_class_a_beta-lactamase_ctx-m-14                |                    15 |                    1 |    13.7333  |              0.250956 |


## Conclusión general para redactar en resultados

En ambos antibióticos, los modelos basados en presencia y ausencia de genes AMR lograron discriminar entre genomas resistentes y susceptibles con desempeños consistentes tanto en validación cruzada como en evaluación independiente. La incorporación de XGBoost permite contrastar los modelos lineales y de ensamble clásico con un modelo de boosting de gradiente, por lo que la comparación final es más robusta. Además, los genes destacados por los modelos seleccionados mostraron estabilidad a lo largo de múltiples repeticiones, lo que fortalece su interpretación biológica.

## Limitaciones que todavía deben mencionarse

Aunque estos resultados pueden reportarse como resultados finales del modelo base, todavía conviene declarar tres limitaciones. La primera es la ausencia de una cohorte externa completamente independiente. La segunda es que el análisis se basa en genes AMR anotados y no incorpora mutaciones puntuales no anotadas. La tercera es que no se modeló explícitamente la estructura poblacional o los linajes bacterianos.