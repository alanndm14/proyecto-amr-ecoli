import pandas as pd
import re
from pathlib import Path

# Cambia este nombre por el nombre real de tu archivo grande de Specialty Genes
BIG_FILE = "specialty_genes_ecoli_amr.csv"

# Archivo pequeño con los genome_id que sí nos importan
IDS_FILE = "ecoli_valid_genome_ids.txt"

# Archivo de salida
OUT_FILE = "specialty_genes_ecoli_amr_filtered_to_project_genomes.csv"


def normalize_columns(df):
    df = df.copy()
    df.columns = [
        re.sub(r"_+", "_", re.sub(r"[^0-9a-zA-Z]+", "_", c.strip().lower())).strip("_")
        for c in df.columns
    ]
    return df


valid_ids = set(
    pd.read_csv(IDS_FILE, header=None, dtype=str)[0]
    .dropna()
    .astype(str)
    .str.strip()
)

print(f"Genome IDs válidos cargados: {len(valid_ids):,}")

first = True
total_in = 0
total_out = 0

for chunk in pd.read_csv(BIG_FILE, chunksize=100_000, dtype=str, low_memory=False):
    total_in += len(chunk)
    chunk = normalize_columns(chunk)

    if "genome_id" not in chunk.columns:
        raise ValueError(f"No encontré columna genome_id. Columnas disponibles: {list(chunk.columns)}")

    chunk["genome_id"] = chunk["genome_id"].astype(str).str.strip()

    # 1. Conservar solo genomas que tienen fenotipo limpio
    chunk = chunk[chunk["genome_id"].isin(valid_ids)].copy()

    # 2. Si existe columna property, conservar solo Antibiotic Resistance / AMR
    if "property" in chunk.columns:
        prop = chunk["property"].fillna("").astype(str)
        mask = prop.str.contains("antibiotic resistance|amr", case=False, regex=True)
        if mask.any():
            chunk = chunk[mask].copy()

    total_out += len(chunk)

    chunk.to_csv(
        OUT_FILE,
        index=False,
        mode="w" if first else "a",
        header=first
    )
    first = False

print(f"Filas leídas del archivo grande: {total_in:,}")
print(f"Filas conservadas: {total_out:,}")
print(f"Archivo filtrado creado: {OUT_FILE}")
