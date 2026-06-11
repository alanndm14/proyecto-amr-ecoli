## Tabla 1. Composición final de los conjuntos de datos usados para el análisis predictivo
| Antibiótico   | Genomas totales   | Resistentes   | Susceptibles   |   Genes AMR model_ready |
|:--------------|:------------------|:--------------|:---------------|------------------------:|
| Ciprofloxacin | 7,313             | 1,851         | 5,462          |                     633 |
| Cefotaxime    | 6,590             | 1,272         | 5,318          |                     560 |

## Tabla 2. Desempeño de los modelos en el conjunto retenido de prueba
| Antibiótico   | Modelo                              |   Balanced accuracy |   ROC-AUC |   F1 resistente |   Recall resistente |   Especificidad susceptible |   Falsos susceptibles |   Falsos resistentes |
|:--------------|:------------------------------------|--------------------:|----------:|----------------:|--------------------:|----------------------------:|----------------------:|---------------------:|
| Ciprofloxacin | Regresión logística balanceada      |               0.844 |     0.906 |           0.759 |               0.784 |                       0.905 |                    80 |                  104 |
| Ciprofloxacin | Random forest balanceado            |               0.837 |     0.922 |           0.753 |               0.765 |                       0.909 |                    87 |                   99 |
| Ciprofloxacin | XGBoost ponderado                   |               0.836 |     0.907 |           0.753 |               0.759 |                       0.913 |                    89 |                   95 |
| Ciprofloxacin | MLP ponderado (promedio 7 semillas) |               0.835 |     0.903 |           0.734 |               0.797 |                       0.873 |                    75 |                  139 |
| Cefotaxime    | Regresión logística balanceada      |               0.827 |     0.893 |           0.719 |               0.724 |                       0.93  |                    70 |                   74 |
| Cefotaxime    | Random forest balanceado            |               0.791 |     0.879 |           0.623 |               0.728 |                       0.854 |                    69 |                  155 |
| Cefotaxime    | XGBoost ponderado                   |               0.81  |     0.886 |           0.69  |               0.697 |                       0.923 |                    77 |                   82 |
| Cefotaxime    | MLP ponderado (promedio 7 semillas) |               0.818 |     0.893 |           0.676 |               0.744 |                       0.891 |                    65 |                  116 |

## Tabla 3. Validación leave-method-out por método fenotípico de laboratorio
| Antibiótico   | Método dejado fuera para prueba   | n     |   Resistentes |   Balanced accuracy |   ROC-AUC |   Recall resistente |   Falsos susceptibles |   Falsos resistentes |
|:--------------|:----------------------------------|:------|--------------:|--------------------:|----------:|--------------------:|----------------------:|---------------------:|
| Ciprofloxacin | Broth dilution                    | 2,494 |           960 |               0.786 |     0.858 |               0.782 |                   209 |                  322 |
| Ciprofloxacin | Disk diffusion                    | 3,709 |           675 |               0.815 |     0.869 |               0.716 |                   192 |                  257 |
| Ciprofloxacin | MIC                               | 1,074 |           186 |               0.813 |     0.863 |               0.72  |                    52 |                   83 |
| Cefotaxime    | Agar dilution                     | 1,015 |            81 |               0.821 |     0.878 |               0.691 |                    25 |                   47 |
| Cefotaxime    | Broth dilution                    | 2,206 |           829 |               0.768 |     0.829 |               0.651 |                   289 |                  160 |
| Cefotaxime    | Disk diffusion                    | 3,326 |           330 |               0.801 |     0.848 |               0.721 |                    92 |                  357 |

## Tabla 4. Genes AMR más estables del modelo seleccionado
| Antibiótico   | Gen/señal AMR                                                   |   Frecuencia de selección |   Fracción de selección |   Rango medio |   Importancia media absoluta |
|:--------------|:----------------------------------------------------------------|--------------------------:|------------------------:|--------------:|-----------------------------:|
| Cefotaxime    | ctx-m-15                                                        |                        15 |                       1 |         1.4   |                        0.652 |
| Cefotaxime    | ctx-m_family                                                    |                        15 |                       1 |         1.733 |                        0.644 |
| Cefotaxime    | cmy-2                                                           |                        15 |                       1 |         3.467 |                        0.525 |
| Cefotaxime    | ctx-m-14                                                        |                        15 |                       1 |         4.667 |                        0.512 |
| Cefotaxime    | extended-spectrum_class_a_beta-lactamase_ctx-m-15               |                        15 |                       1 |         4.8   |                        0.491 |
| Cefotaxime    | multispecies_class_c_beta-lactamase_cmy-2                       |                        15 |                       1 |         6.133 |                        0.461 |
| Cefotaxime    | multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-15  |                        15 |                       1 |         6.6   |                        0.454 |
| Cefotaxime    | multispecies_class_a_extended-spectrum_beta-lactamase_ctx-m-14  |                        15 |                       1 |         8.933 |                        0.392 |
| Ciprofloxacin | escherichia_quinolone_resistant_parc                            |                        15 |                       1 |         1.2   |                        0.812 |
| Ciprofloxacin | escherichia_coli_parc_conferring_resistance_to_fluoroquinolones |                        15 |                       1 |         1.867 |                        0.784 |
| Ciprofloxacin | escherichia_quinolone_resistant_gyra                            |                        15 |                       1 |         2.933 |                        0.692 |
| Ciprofloxacin | escherichia_coli_gyra_conferring_resistance_to_fluoroquinolones |                        15 |                       1 |         4     |                        0.45  |
| Ciprofloxacin | catb3                                                           |                        15 |                       1 |         5.333 |                        0.35  |
| Ciprofloxacin | escherichia_quinolone_resistant_pare                            |                        15 |                       1 |         5.867 |                        0.317 |
| Ciprofloxacin | extended-spectrum_class_a_beta-lactamase_ctx-m-15               |                        15 |                       1 |         9.133 |                        0.255 |
| Ciprofloxacin | blaec                                                           |                        15 |                       1 |         9.4   |                        0.255 |

## Tabla A1. Resumen del MLP por arquitectura y semillas
| Antibiótico   | Arquitectura MLP   |   Balanced accuracy media |   DE balanced accuracy |   ROC-AUC medio |   F1 resistente medio |   Recall resistente medio |   Falsos susceptibles promedio |   Falsos resistentes promedio |   Semillas |
|:--------------|:-------------------|--------------------------:|-----------------------:|----------------:|----------------------:|--------------------------:|-------------------------------:|------------------------------:|-----------:|
| Cefotaxime    | MLP-32             |                     0.807 |                  0.014 |           0.883 |                 0.659 |                     0.732 |                         68     |                       125     |          7 |
| Cefotaxime    | MLP-16             |                     0.783 |                  0.013 |           0.861 |                 0.608 |                     0.723 |                         70.286 |                       167.857 |          7 |
| Cefotaxime    | MLP-16-8           |                     0.754 |                  0.113 |           0.823 |                 0.593 |                     0.764 |                         60     |                       271.143 |          7 |
| Ciprofloxacin | MLP-32             |                     0.831 |                  0.006 |           0.898 |                 0.726 |                     0.794 |                         76.286 |                       145.143 |          7 |
| Ciprofloxacin | MLP-16             |                     0.825 |                  0.007 |           0.896 |                 0.723 |                     0.775 |                         83.143 |                       136.286 |          7 |
| Ciprofloxacin | MLP-16-8           |                     0.782 |                  0.124 |           0.849 |                 0.681 |                     0.813 |                         69.286 |                       272.714 |          7 |