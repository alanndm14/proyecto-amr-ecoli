# Análisis de robustez por método de laboratorio

Modelo usado: regresión logística balanceada con regularización L1 y `C=0.1`. Para esta comparación rápida se usó `solver=liblinear`, que reproduce casi exactamente el desempeño holdout del modelo seleccionado (`saga`, L1, C=0.1) y permite hacer las particiones por método sin esperar una nueva búsqueda de hiperparámetros.

Se hicieron dos comparaciones: (1) desempeño por método dentro del holdout original 80/20, y (2) validación leave-method-out, donde se entrena sin un método de laboratorio y se prueba únicamente en ese método.

## ciprofloxacin

### Holdout por método

| laboratory_typing_method   |    n |   n_resistant |   balanced_accuracy |   roc_auc |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |
|:---------------------------|-----:|--------------:|--------------------:|----------:|-------------------:|--------------------------:|--------------------:|------------------:|
| TOTAL_HOLDOUT              | 1463 |           370 |               0.845 |     0.906 |              0.786 |                     0.904 |                  79 |               105 |
| Broth dilution             |  541 |           205 |               0.824 |     0.889 |              0.79  |                     0.857 |                  43 |                48 |
| Disk diffusion             |  703 |           122 |               0.864 |     0.925 |              0.803 |                     0.924 |                  24 |                44 |
| MIC                        |  213 |            38 |               0.831 |     0.883 |              0.737 |                     0.926 |                  10 |                13 |


### Leave-method-out

| laboratory_typing_method   |    n |   n_resistant |   balanced_accuracy |   roc_auc |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |
|:---------------------------|-----:|--------------:|--------------------:|----------:|-------------------:|--------------------------:|--------------------:|------------------:|
| Broth dilution             | 2494 |           960 |               0.786 |     0.858 |              0.782 |                     0.79  |                 209 |               322 |
| Disk diffusion             | 3709 |           675 |               0.815 |     0.869 |              0.716 |                     0.915 |                 192 |               257 |
| MIC                        | 1074 |           186 |               0.813 |     0.863 |              0.72  |                     0.907 |                  52 |                83 |


## cefotaxime

### Holdout por método

| laboratory_typing_method   |    n |   n_resistant |   balanced_accuracy |   roc_auc |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |
|:---------------------------|-----:|--------------:|--------------------:|----------:|-------------------:|--------------------------:|--------------------:|------------------:|
| TOTAL_HOLDOUT              | 1318 |           254 |               0.829 |     0.894 |              0.728 |                     0.93  |                  69 |                74 |
| Agar dilution              |  204 |            18 |               0.817 |     0.892 |              0.667 |                     0.968 |                   6 |                 6 |
| Broth dilution             |  451 |           155 |               0.818 |     0.895 |              0.748 |                     0.889 |                  39 |                33 |
| Disk diffusion             |  655 |            75 |               0.817 |     0.867 |              0.693 |                     0.941 |                  23 |                34 |


### Leave-method-out

| laboratory_typing_method   |    n |   n_resistant |   balanced_accuracy |   roc_auc |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |
|:---------------------------|-----:|--------------:|--------------------:|----------:|-------------------:|--------------------------:|--------------------:|------------------:|
| Agar dilution              | 1015 |            81 |               0.821 |     0.878 |              0.691 |                     0.95  |                  25 |                47 |
| Broth dilution             | 2206 |           829 |               0.768 |     0.829 |              0.651 |                     0.884 |                 289 |               160 |
| Disk diffusion             | 3326 |           330 |               0.801 |     0.848 |              0.721 |                     0.881 |                  92 |               357 |

