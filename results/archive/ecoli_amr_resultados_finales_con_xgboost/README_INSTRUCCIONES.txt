# Instrucciones rápidas

1. Coloca este script en la misma carpeta que:
   - dataset_ciprofloxacin_model_ready.csv
   - dataset_cefotaxime_model_ready.csv

2. Instala dependencias:
   pip install pandas scikit-learn matplotlib tabulate xgboost

   En Jupyter:
   %pip install pandas scikit-learn matplotlib tabulate xgboost

3. Ejecuta:
   python generate_final_nonprelim_analysis_es_xgboost.py

4. El script creará:
   - carpeta ecoli_amr_resultados_finales_con_xgboost
   - archivo ecoli_amr_resultados_finales_con_xgboost.zip

Tiempo esperado:
- Esta versión puede tardar más que la anterior porque XGBoost agrega más combinaciones de hiperparámetros.
- En una laptop moderna puede tardar varios minutos.
- Si tarda demasiado, reduce el grid de XGBoost en el script.