# Resultados actuales del modelo base

Estos resultados ya pueden presentarse como **resultados del modelo base**. No los llames necesariamente preliminares si el proyecto es para clase; llámalos “modelo base” o “análisis predictivo inicial”.

## Métricas obtenidas

| Antibiótico | Modelo | Accuracy | Balanced accuracy | F1 Resistant | ROC-AUC |
|---|---:|---:|---:|---:|---:|
| ciprofloxacin | Regresión logística balanceada | 0.8619 | 0.8289 | 0.7363 | 0.8787 |
| ciprofloxacin | Random forest balanceado | 0.8722 | 0.8412 | 0.7549 | 0.9218 |
| cefotaxime | Regresión logística balanceada | 0.8771 | 0.8190 | 0.6943 | 0.8697 |
| cefotaxime | Random forest balanceado | 0.8953 | 0.7928 | 0.6974 | 0.9086 |

## Interpretación breve

Para ciprofloxacin, random forest fue el mejor modelo general, ya que obtuvo la mayor balanced accuracy, F1 para resistentes y ROC-AUC. Para cefotaxime, random forest tuvo mayor accuracy y ROC-AUC, pero la regresión logística presentó mejor balanced accuracy y menor número de falsos susceptibles. En resistencia antimicrobiana, este último punto es importante porque clasificar un genoma resistente como susceptible puede ser el error más delicado.

## Texto para pegar en resultados

Se entrenaron modelos de clasificación supervisada para predecir resistencia antimicrobiana en *Escherichia coli* frente a ciprofloxacin y cefotaxime a partir de una matriz binaria de presencia/ausencia de genes AMR. Para cada antibiótico se evaluaron dos modelos: regresión logística balanceada y random forest balanceado. La variable objetivo se codificó como 1 para genomas resistentes y 0 para genomas susceptibles.

En ciprofloxacin, random forest presentó el mejor desempeño general, con accuracy de 0.8722, balanced accuracy de 0.8412, F1 para resistentes de 0.7549 y ROC-AUC de 0.9218. Estos valores indican que el modelo logró distinguir de manera adecuada entre genomas resistentes y susceptibles.

En cefotaxime, random forest obtuvo la mayor accuracy y ROC-AUC, con valores de 0.8953 y 0.9086, respectivamente. Sin embargo, la regresión logística presentó una balanced accuracy mayor, de 0.8190 frente a 0.7928, y produjo menos falsos susceptibles. Por ello, cuando la prioridad es detectar genomas resistentes, la regresión logística puede considerarse más adecuada para cefotaxime.

En conjunto, los resultados muestran que los genes AMR contienen señal predictiva útil para clasificar genomas de *Escherichia coli* como resistentes o susceptibles a ambos antibióticos.
