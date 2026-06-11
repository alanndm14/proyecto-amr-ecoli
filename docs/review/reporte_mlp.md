# Resultados del MLP ponderado

Se evaluó un perceptrón multicapa (MLP) ponderado por clase usando el mismo holdout 80/20 estratificado del análisis principal. Como el MLP depende de inicialización aleatoria, se probaron 7 semillas por configuración y también se calculó un promedio de probabilidades entre semillas como análisis de estabilidad.

## Resumen por configuración y semillas

| antibiotic    | candidate_name   | hidden_layer_sizes   |   balanced_accuracy_mean |   balanced_accuracy_std |   roc_auc_mean |   roc_auc_std |   f1_resistant_mean |   f1_resistant_std |   recall_resistant_mean |   recall_resistant_std |   false_susceptible_mean |   false_resistant_mean |   n_seeds |
|:--------------|:-----------------|:---------------------|-------------------------:|------------------------:|---------------:|--------------:|--------------------:|-------------------:|------------------------:|-----------------------:|-------------------------:|-----------------------:|----------:|
| cefotaxime    | MLP-32           | (32,)                |                 0.807401 |              0.0136303  |       0.882699 |    0.00854697 |            0.659112 |         0.0271949  |                0.732283 |              0.0140119 |                  68      |                125     |         7 |
| cefotaxime    | MLP-16           | (16,)                |                 0.782762 |              0.013174   |       0.860594 |    0.0148359  |            0.607617 |         0.0249514  |                0.723285 |              0.0204212 |                  70.2857 |                167.857 |         7 |
| cefotaxime    | MLP-16-8         | (16, 8)              |                 0.754473 |              0.113098   |       0.822539 |    0.140529   |            0.593268 |         0.121644   |                0.76378  |              0.106324  |                  60      |                271.143 |         7 |
| ciprofloxacin | MLP-32           | (32,)                |                 0.830515 |              0.00603095 |       0.897857 |    0.00418328 |            0.726266 |         0.00851035 |                0.793822 |              0.0111123 |                  76.2857 |                145.143 |         7 |
| ciprofloxacin | MLP-16           | (16,)                |                 0.8253   |              0.00676235 |       0.896222 |    0.0045742  |            0.723329 |         0.00873179 |                0.77529  |              0.0137306 |                  83.1429 |                136.286 |         7 |
| ciprofloxacin | MLP-16-8         | (16, 8)              |                 0.781616 |              0.124275   |       0.848729 |    0.12766    |            0.681072 |         0.122365   |                0.812741 |              0.0847189 |                  69.2857 |                272.714 |         7 |


## MLP ensemble por configuración

| antibiotic    | candidate_name   |   balanced_accuracy |   roc_auc |   f1_resistant |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant | confusion_matrix        |
|:--------------|:-----------------|--------------------:|----------:|---------------:|-------------------:|--------------------------:|--------------------:|------------------:|:------------------------|
| cefotaxime    | MLP-32           |            0.817536 |  0.89275  |       0.676208 |           0.744094 |                  0.890977 |                  65 |               116 | [[948, 116], [65, 189]] |
| cefotaxime    | MLP-16-8         |            0.810842 |  0.884369 |       0.649832 |           0.759843 |                  0.861842 |                  61 |               147 | [[917, 147], [61, 193]] |
| cefotaxime    | MLP-16           |            0.790306 |  0.87275  |       0.622673 |           0.724409 |                  0.856203 |                  70 |               153 | [[911, 153], [70, 184]] |
| ciprofloxacin | MLP-32           |            0.835062 |  0.902502 |       0.733831 |           0.797297 |                  0.872827 |                  75 |               139 | [[954, 139], [75, 295]] |
| ciprofloxacin | MLP-16-8         |            0.831029 |  0.902322 |       0.730238 |           0.786486 |                  0.875572 |                  79 |               136 | [[957, 136], [79, 291]] |
| ciprofloxacin | MLP-16           |            0.829262 |  0.901029 |       0.730038 |           0.778378 |                  0.880146 |                  82 |               131 | [[962, 131], [82, 288]] |


## Comparación holdout contra modelos previos + mejor MLP ensemble

| antibiotic    | model_label                         |   balanced_accuracy |   roc_auc |   f1_resistant |   recall_resistant |   specificity_susceptible |   false_susceptible |   false_resistant |
|:--------------|:------------------------------------|--------------------:|----------:|---------------:|-------------------:|--------------------------:|--------------------:|------------------:|
| cefotaxime    | Regresión logística balanceada      |            0.82743  |  0.893371 |       0.71875  |           0.724409 |                  0.930451 |                  70 |                74 |
| cefotaxime    | MLP ponderado (promedio 7 semillas) |            0.817536 |  0.89275  |       0.676208 |           0.744094 |                  0.890977 |                  65 |               116 |
| cefotaxime    | XGBoost ponderado                   |            0.809891 |  0.886171 |       0.690058 |           0.69685  |                  0.922932 |                  77 |                82 |
| cefotaxime    | Bosque aleatorio balanceado         |            0.791335 |  0.878893 |       0.622896 |           0.728346 |                  0.854323 |                  69 |               155 |
| ciprofloxacin | Regresión logística balanceada      |            0.844316 |  0.90585  |       0.759162 |           0.783784 |                  0.904849 |                  80 |               104 |
| ciprofloxacin | Bosque aleatorio balanceado         |            0.837144 |  0.922398 |       0.75266  |           0.764865 |                  0.909424 |                  87 |                99 |
| ciprofloxacin | XGBoost ponderado                   |            0.836271 |  0.907045 |       0.753351 |           0.759459 |                  0.913083 |                  89 |                95 |
| ciprofloxacin | MLP ponderado (promedio 7 semillas) |            0.835062 |  0.902502 |       0.733831 |           0.797297 |                  0.872827 |                  75 |               139 |


## Interpretación

El MLP funciona como red neuronal solicitada, pero no superó a la regresión logística balanceada en balanced accuracy de holdout. El resultado más sólido es conservar el MLP como cuarto modelo comparativo y mantener la regresión logística balanceada como modelo final seleccionado. Esto fortalece el reporte porque demuestra que sí se evaluó una red neuronal y que la elección final no fue arbitraria.
