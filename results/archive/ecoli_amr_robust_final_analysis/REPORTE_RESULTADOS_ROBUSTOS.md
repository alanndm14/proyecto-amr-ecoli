# Resultados finales robustecidos del modelo AMR en *Escherichia coli*

## Auditoría de datos

| antibiotic    |   n_samples |   n_features |   n_susceptible |   n_resistant |   resistant_fraction |   duplicated_genome_ids |
|:--------------|------------:|-------------:|----------------:|--------------:|---------------------:|------------------------:|
| ciprofloxacin |        7313 |          633 |            5462 |          1851 |             0.253111 |                       0 |
| cefotaxime    |        6590 |          560 |            5318 |          1272 |             0.19302  |                       0 |

## Fortalecimiento del análisis

Se añadieron tres pasos para fortalecer el análisis: ajuste pequeño de hiperparámetros por validación cruzada interna, evaluación por validación cruzada estratificada de tres particiones y análisis de estabilidad de genes importantes. La métrica principal para seleccionar modelos fue balanced accuracy, porque las clases resistente y susceptible están desbalanceadas. También se reportó recall de resistentes, porque un falso susceptible es el error más delicado en resistencia antimicrobiana.

## Validación cruzada de modelos ajustados

| antibiotic    | model                        |   balanced_accuracy_mean |   balanced_accuracy_std |   f1_resistant_mean |   f1_resistant_std |   recall_resistant_mean |   recall_resistant_std |   roc_auc_mean |   roc_auc_std |   specificity_susceptible_mean |   specificity_susceptible_std |
|:--------------|:-----------------------------|-------------------------:|------------------------:|--------------------:|-------------------:|------------------------:|-----------------------:|---------------:|--------------:|-------------------------------:|------------------------------:|
| ciprofloxacin | logistic_regression_balanced |                 0.844388 |              0.0096722  |            0.761821 |         0.0157128  |                0.779038 |              0.0117248 |       0.910821 |    0.0124879  |                       0.909738 |                    0.00854074 |
| ciprofloxacin | random_forest_balanced       |                 0.836143 |              0.0152365  |            0.753818 |         0.0211539  |                0.757969 |              0.0269421 |       0.921714 |    0.0168988  |                       0.914317 |                    0.00713899 |
| cefotaxime    | logistic_regression_balanced |                 0.825863 |              0.00520346 |            0.706281 |         0.00813548 |                0.734277 |              0.0134109 |       0.886937 |    0.00970648 |                       0.917449 |                    0.00738775 |
| cefotaxime    | random_forest_balanced       |                 0.76464  |              0.0172529  |            0.64916  |         0.0252992  |                0.577044 |              0.0365629 |       0.891792 |    0.00264759 |                       0.952237 |                    0.00417416 |

## Evaluación en conjunto de prueba retenido

| antibiotic    | model                        |   accuracy |   balanced_accuracy |   f1_resistant |   roc_auc |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |   true_resistant |   true_susceptible |
|:--------------|:-----------------------------|-----------:|--------------------:|---------------:|----------:|-------------------:|--------------------------:|--------------------:|------------------:|-----------------:|-------------------:|
| ciprofloxacin | logistic_regression_balanced |   0.874231 |            0.84521  |       0.759791 |  0.905984 |           0.786486 |                  0.903934 |                  79 |               105 |              291 |                988 |
| ciprofloxacin | random_forest_balanced       |   0.872864 |            0.839826 |       0.754617 |  0.919962 |           0.772973 |                  0.906679 |                  84 |               102 |              286 |                991 |
| cefotaxime    | logistic_regression_balanced |   0.891502 |            0.829399 |       0.721248 |  0.893964 |           0.728346 |                  0.930451 |                  69 |                74 |              185 |                990 |
| cefotaxime    | random_forest_balanced       |   0.881639 |            0.79182  |       0.677686 |  0.899669 |           0.645669 |                  0.93797  |                  90 |                66 |              164 |                998 |

## Modelos seleccionados

| antibiotic    | model                        |   balanced_accuracy_mean |   balanced_accuracy_std |   recall_resistant_mean |   roc_auc_mean |   holdout_false_susceptible |   holdout_false_resistant | selection_rule                      |
|:--------------|:-----------------------------|-------------------------:|------------------------:|------------------------:|---------------:|----------------------------:|--------------------------:|:------------------------------------|
| ciprofloxacin | logistic_regression_balanced |                   0.8444 |                  0.0097 |                  0.779  |         0.9108 |                          79 |                       105 | highest 3-fold CV balanced accuracy |
| cefotaxime    | logistic_regression_balanced |                   0.8259 |                  0.0052 |                  0.7343 |         0.8869 |                          69 |                        74 | highest 3-fold CV balanced accuracy |

## Genes más estables en los modelos seleccionados

### ciprofloxacin — logistic_regression_balanced

| gene                                                            |   selection_frequency |   mean_rank |   mean_abs_importance |
|:----------------------------------------------------------------|----------------------:|------------:|----------------------:|
| escherichia_quinolone_resistant_parc                            |                     3 |     1.33333 |              0.798105 |
| escherichia_coli_parc_conferring_resistance_to_fluoroquinolones |                     3 |     1.66667 |              0.764425 |
| escherichia_quinolone_resistant_gyra                            |                     3 |     3       |              0.700732 |
| escherichia_coli_gyra_conferring_resistance_to_fluoroquinolones |                     3 |     4       |              0.450544 |
| escherichia_quinolone_resistant_pare                            |                     3 |     6       |              0.323412 |
| catb3                                                           |                     3 |     6       |              0.319879 |
| extended-spectrum_class_a_beta-lactamase_ctx-m-15               |                     3 |    10       |              0.237997 |
| ctx-m-15                                                        |                     3 |    12.3333  |              0.226019 |
| escherichia_fosfomycin_resistant_ptsi                           |                     3 |    14       |              0.205083 |
| vgac                                                            |                     3 |    14.6667  |              0.207805 |
| blaec                                                           |                     2 |     6.5     |              0.291245 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-15  |                     2 |     7       |              0.308646 |
### cefotaxime — logistic_regression_balanced

| gene                                                             |   selection_frequency |   mean_rank |   mean_abs_importance |
|:-----------------------------------------------------------------|----------------------:|------------:|----------------------:|
| ctx-m_family                                                     |                     3 |     1.66667 |              0.642848 |
| ctx-m-15                                                         |                     3 |     1.66667 |              0.640607 |
| cmy-2                                                            |                     3 |     3.33333 |              0.51055  |
| ctx-m-14                                                         |                     3 |     5       |              0.475215 |
| extended-spectrum_class_a_beta-lactamase_ctx-m-15                |                     3 |     5       |              0.459928 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-15   |                     3 |     6.66667 |              0.452187 |
| multispecies_class_c_beta-lactamase_cmy-2                        |                     3 |     6.66667 |              0.444285 |
| class_a_beta-lactamase_ec_3.5.2.6_ctx-m_family_extended-spectrum |                     3 |     7.66667 |              0.383784 |
| extended-spectrum_class_c_beta-lactamase_cmy-2                   |                     3 |     8.66667 |              0.378047 |
| multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-14   |                     3 |     9       |              0.365439 |
| cmycmy-2cfelat_family                                            |                     3 |    10.6667  |              0.356064 |
| escherichia_quinolone_resistant_gyra                             |                     3 |    13       |              0.271919 |

## Interpretación lista para informe

Los resultados muestran que la presencia/ausencia de genes AMR permite clasificar genomas de *Escherichia coli* como resistentes o susceptibles frente a ciprofloxacin y cefotaxime. La validación cruzada estratificada confirma que el desempeño no depende únicamente de una partición entrenamiento-prueba.

Para ciprofloxacin, los genes importantes se relacionan con mecanismos conocidos de resistencia a fluoroquinolonas, incluyendo genes asociados con topoisomerasas y determinantes plasmídicos. Para cefotaxime, los genes más relevantes pertenecen principalmente a familias de beta-lactamasas, especialmente genes tipo CTX-M y otros determinantes relacionados con resistencia a cefalosporinas.

Estos resultados pueden presentarse como resultados principales del modelo base robustecido. La principal limitación restante es que el análisis usa genes AMR anotados y todavía no incorpora estructura poblacional, linajes bacterianos, mutaciones puntuales no anotadas ni validación externa en una cohorte independiente.