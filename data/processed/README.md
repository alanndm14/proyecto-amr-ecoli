# Datos procesados

## Datos finales usados por los modelos

- `dataset_ciprofloxacin_model_ready.csv`: 7,313 genomas y 633 genes AMR.
- `dataset_cefotaxime_model_ready.csv`: 6,590 genomas y 560 genes AMR.

La variable objetivo es `y`:

- `1`: Resistant.
- `0`: Susceptible.

Las variables predictoras comienzan con `amr__` y representan presencia/ausencia.

## Datos intermedios

- `clean_phenotypes_lab_only_cipro_cefotaxime.csv`: fenotipos limpios.
- `clean_specialty_genes_amr_project_genomes.csv`: anotaciones AMR limpias.

`scripts/prepare_model_data.py` reconstruye la matriz binaria, los datasets
completos y los datasets `model_ready` a partir de los dos archivos limpios.
Las matrices completas generadas no se versionan porque son reconstruibles y
duplicarían decenas de megabytes.
