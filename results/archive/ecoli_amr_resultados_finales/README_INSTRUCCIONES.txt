# Instrucciones rápidas

1. Coloca este script en la misma carpeta que:
   - dataset_ciprofloxacin_model_ready.csv
   - dataset_cefotaxime_model_ready.csv
2. Instala dependencias si hace falta:
   pip install pandas scikit-learn matplotlib tabulate
3. Ejecuta:
   python generate_final_nonprelim_analysis_es.py
4. El script creará:
   - carpeta ecoli_amr_resultados_finales
   - archivo ecoli_amr_resultados_finales.zip

Tiempo esperado
- En una laptop moderna, entre ~1 y 6 minutos es normal.
- Si tarda menos, también puede ser normal si tu CPU tiene varios núcleos y los archivos ya están limpios.