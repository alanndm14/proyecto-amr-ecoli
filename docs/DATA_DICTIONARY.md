# Diccionario de datos

## Datasets `model_ready`

Los archivos principales son:

- `data/processed/dataset_ciprofloxacin_model_ready.csv`
- `data/processed/dataset_cefotaxime_model_ready.csv`

Cada fila representa un genoma de *Escherichia coli* con fenotipo AMR.

| Campo | Descripción |
|---|---|
| `genome_id` | Identificador del genoma en BV-BRC |
| `genome_name` | Nombre del genoma |
| `antibiotic` | Antibiótico evaluado |
| `resistant_phenotype` | Fenotipo textual Resistant/Susceptible |
| `y` | Variable objetivo: 1 resistente, 0 susceptible |
| `measurement*` | Medición fenotípica disponible |
| `laboratory_typing_method` | Método de laboratorio usado |
| `testing_standard*` | Estándar de interpretación |
| `evidence` | Evidencia asociada |
| `source` | Fuente del registro |
| `amr__*` | Variable binaria de presencia/ausencia de un gen o señal AMR |

## Dimensiones finales

| Antibiótico | Genomas | Resistentes | Susceptibles | Variables AMR |
|---|---:|---:|---:|---:|
| Ciprofloxacin | 7,313 | 1,851 | 5,462 | 633 |
| Cefotaxime | 6,590 | 1,272 | 5,318 | 560 |

## Datos limpios de origen

- `clean_phenotypes_lab_only_cipro_cefotaxime.csv`: registros fenotípicos
  filtrados a evidencia de laboratorio y clases Resistant/Susceptible.
- `clean_specialty_genes_amr_project_genomes.csv`: pares genoma-predictor AMR
  usados para reconstruir las matrices binarias.
