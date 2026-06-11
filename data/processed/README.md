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
- `gene_presence_matrix_all_project_genomes.csv`: matriz binaria completa.
- `dataset_*_gene_presence_full.csv`: cruces completos por antibiótico.
- `specialty_genes_ecoli_amr_filtered_to_project_genomes.csv`: exportación filtrada.

`scripts/prepare_model_data.py` reconstruye la matriz binaria, los datasets
completos y los datasets `model_ready` a partir de los dos archivos limpios.
