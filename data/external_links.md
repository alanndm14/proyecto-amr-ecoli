# Fuentes externas y archivos grandes

## BV-BRC

Los datos de genomas, fenotipos de susceptibilidad antimicrobiana y anotaciones
de genes AMR se obtuvieron de BV-BRC:

- Portal: https://www.bv-brc.org/
- Organismo: *Escherichia coli*
- Antibióticos: ciprofloxacin y cefotaxime
- Clases fenotípicas conservadas: Resistant y Susceptible

## Archivo crudo no versionado

| Archivo | Tamaño local aproximado | Manejo recomendado |
|---|---:|---|
| `specialty_genes_ecoli_amr.csv` | 915 MB | Git LFS, almacenamiento institucional o nueva exportación de BV-BRC |

Ejemplo para Git LFS:

```bash
git lfs install
git lfs track "data/raw/*.csv"
```

