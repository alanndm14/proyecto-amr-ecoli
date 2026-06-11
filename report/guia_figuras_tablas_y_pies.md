# Tablas y figuras finales para el reporte AMR

## Orden recomendado en el reporte

**Diseño / Planteamiento:** Figura 1 y Tabla 1.

**Resultados:** Tabla 2, Figura 2, Figura 3, Figura 4, Figuras 5 y 6 y Tabla 3.

**Anexo opcional:** Tabla A1 del MLP por arquitectura/semillas.

## Tabla 1. Composición final de los conjuntos de datos usados para el análisis predictivo

| Antibiótico   | Genomas totales   | Resistentes   | Susceptibles   |   Genes AMR model_ready |
|:--------------|:------------------|:--------------|:---------------|------------------------:|
| Ciprofloxacin | 7,313             | 1,851         | 5,462          |                     633 |
| Cefotaxime    | 6,590             | 1,272         | 5,318          |                     560 |


**Nota.** Los genes AMR corresponden a variables binarias de presencia/ausencia conservadas después del filtrado para generar los datasets `model_ready`.

## Figura 1. Flujo general del análisis genómico predictivo

**Archivo:** `figures/figura_1_flujo_general.png`

**Referencia en texto:** “El flujo general del procesamiento de datos, construcción de matrices, entrenamiento y evaluación de modelos se muestra en la Figura 1.”

**Pie de figura:** El esquema resume el proceso seguido desde la obtención de datos en BV-BRC hasta la interpretación biológica de los modelos. Primero se filtraron los fenotipos AMR de *E. coli* para ciprofloxacin y cefotaxime, conservando únicamente registros Resistant/Susceptible con evidencia de laboratorio. Después se integraron genes AMR de la categoría Antibiotic Resistance, se construyeron matrices binarias de presencia/ausencia y se generaron datasets `model_ready` por antibiótico. Finalmente, se compararon modelos supervisados y se evaluó su desempeño mediante métricas de clasificación e interpretación biológica.

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


**Nota.** La selección principal del modelo se basó en la balanced accuracy de validación cruzada repetida y se confirmó en el conjunto retenido. El MLP corresponde al promedio de probabilidades de siete semillas para la mejor arquitectura evaluada.

## Figura 2. Comparación de modelos por balanced accuracy

**Archivo:** `figures/figura_2_comparacion_balanced_accuracy_modelos.png`

**Referencia en texto:** “La comparación de balanced accuracy entre los cuatro modelos evaluados se resume en la Figura 2.”

**Pie de figura:** La gráfica compara el desempeño balanceado de regresión logística, random forest, XGBoost y MLP en ambos antibióticos. En los dos casos, la regresión logística balanceada se mantuvo como el modelo más consistente según esta métrica, aunque el MLP tuvo valores cercanos y permitió incorporar un enfoque de red neuronal al análisis.

## Figura 3. Falsos susceptibles por modelo

**Archivo:** `figures/figura_3_falsos_susceptibles_modelos.png`

**Referencia en texto:** “Además de la balanced accuracy, se revisaron los falsos susceptibles de cada modelo, como se muestra en la Figura 3.”

**Pie de figura:** Los falsos susceptibles corresponden a genomas resistentes que fueron clasificados como susceptibles. El MLP redujo ligeramente este tipo de error en ambos antibióticos, pero generó más falsos resistentes, por lo que no superó a la regresión logística en el desempeño general.

## Figura 4. Matrices de confusión de los modelos seleccionados

**Archivo:** `figures/figura_4_matrices_confusion_modelos_seleccionados.png`

**Referencia en texto:** “Las matrices de confusión de los modelos seleccionados se muestran en la Figura 4.”

**Pie de figura:** Las matrices muestran los aciertos y errores de la regresión logística balanceada en el conjunto retenido. En ciprofloxacin, el modelo clasificó correctamente 989 susceptibles y 290 resistentes. En cefotaxime, clasificó correctamente 990 susceptibles y 184 resistentes. Los errores restantes corresponden a falsos resistentes y falsos susceptibles.

## Figura 5. Genes AMR más estables en ciprofloxacin

**Archivo:** `figures/figura_5_genes_estables_ciprofloxacin.png`

**Referencia en texto:** “Los genes más estables del modelo seleccionado para ciprofloxacin se presentan en la Figura 5.”

**Pie de figura:** Los genes con mayor frecuencia de selección se relacionaron con resistencia a quinolonas y fluoroquinolonas, especialmente señales asociadas con *parC*, *gyrA* y *parE*. Esta estabilidad indica que esas variables aparecieron de forma repetida entre las más influyentes durante las repeticiones del análisis.

## Figura 6. Genes AMR más estables en cefotaxime

**Archivo:** `figures/figura_6_genes_estables_cefotaxime.png`

**Referencia en texto:** “Para cefotaxime, los genes más estables se muestran en la Figura 6.”

**Pie de figura:** Las variables más estables estuvieron dominadas por genes o familias de β-lactamasas, principalmente CTX-M y CMY. Este patrón coincide con mecanismos conocidos de resistencia a cefalosporinas de tercera generación.

## Tabla 3. Validación leave-method-out por método fenotípico de laboratorio

| Antibiótico   | Método dejado fuera para prueba   | n     |   Resistentes |   Balanced accuracy |   ROC-AUC |   Recall resistente |   Falsos susceptibles |   Falsos resistentes |
|:--------------|:----------------------------------|:------|--------------:|--------------------:|----------:|--------------------:|----------------------:|---------------------:|
| Ciprofloxacin | Broth dilution                    | 2,494 |           960 |               0.786 |     0.858 |               0.782 |                   209 |                  322 |
| Ciprofloxacin | Disk diffusion                    | 3,709 |           675 |               0.815 |     0.869 |               0.716 |                   192 |                  257 |
| Ciprofloxacin | MIC                               | 1,074 |           186 |               0.813 |     0.863 |               0.72  |                    52 |                   83 |
| Cefotaxime    | Agar dilution                     | 1,015 |            81 |               0.821 |     0.878 |               0.691 |                    25 |                   47 |
| Cefotaxime    | Broth dilution                    | 2,206 |           829 |               0.768 |     0.829 |               0.651 |                   289 |                  160 |
| Cefotaxime    | Disk diffusion                    | 3,326 |           330 |               0.801 |     0.848 |               0.721 |                    92 |                  357 |


**Nota.** En esta prueba, cada método de laboratorio se dejó fuera del entrenamiento y después se usó como conjunto de prueba. Por eso representa una evaluación más exigente que una partición aleatoria simple.

## Figura 7. Robustez por método de laboratorio

**Archivo:** `figures/figura_7_robustez_leave_method_out.png`

**Referencia en texto:** “La variación del desempeño bajo la validación leave-method-out se resume en la Figura 7.”

**Pie de figura:** La gráfica muestra que el desempeño se mantuvo razonable al evaluar métodos de laboratorio no vistos durante el entrenamiento, aunque hubo variación entre antibióticos y métodos. El escenario más difícil fue cefotaxime evaluado en broth dilution, con una balanced accuracy menor que las demás condiciones.

## Tabla A1. Resumen del MLP por arquitectura y semillas (anexo opcional)

| Antibiótico   | Arquitectura MLP   |   Balanced accuracy media |   DE balanced accuracy |   ROC-AUC medio |   F1 resistente medio |   Recall resistente medio |   Falsos susceptibles promedio |   Falsos resistentes promedio |   Semillas |
|:--------------|:-------------------|--------------------------:|-----------------------:|----------------:|----------------------:|--------------------------:|-------------------------------:|------------------------------:|-----------:|
| Cefotaxime    | MLP-32             |                     0.807 |                  0.014 |           0.883 |                 0.659 |                     0.732 |                         68     |                       125     |          7 |
| Cefotaxime    | MLP-16             |                     0.783 |                  0.013 |           0.861 |                 0.608 |                     0.723 |                         70.286 |                       167.857 |          7 |
| Cefotaxime    | MLP-16-8           |                     0.754 |                  0.113 |           0.823 |                 0.593 |                     0.764 |                         60     |                       271.143 |          7 |
| Ciprofloxacin | MLP-32             |                     0.831 |                  0.006 |           0.898 |                 0.726 |                     0.794 |                         76.286 |                       145.143 |          7 |
| Ciprofloxacin | MLP-16             |                     0.825 |                  0.007 |           0.896 |                 0.723 |                     0.775 |                         83.143 |                       136.286 |          7 |
| Ciprofloxacin | MLP-16-8           |                     0.782 |                  0.124 |           0.849 |                 0.681 |                     0.813 |                         69.286 |                       272.714 |          7 |


**Nota.** Esta tabla puede ir en anexo si se quiere documentar cómo se eligió el MLP reportado en la Tabla 2.
