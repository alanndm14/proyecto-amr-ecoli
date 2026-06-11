# Datos crudos

El archivo crudo principal de Specialty Genes exportado de BV-BRC no se incluye
en el repositorio Git normal porque mide aproximadamente 915 MB:

```text
specialty_genes_ecoli_amr.csv
```

Debe recuperarse desde BV-BRC o almacenarse con Git LFS/enlace externo. El
archivo procesado y filtrado usado por el proyecto si se incluye en
`data/processed/`.

El repositorio conserva datos procesados suficientes para reproducir las
matrices finales y todos los analisis predictivos.
