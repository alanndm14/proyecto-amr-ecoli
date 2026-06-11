#!/usr/bin/env python3
"""Rebuild gene-presence and model-ready datasets from cleaned project data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ANTIBIOTICS = ("ciprofloxacin", "cefotaxime")
METADATA_COLUMNS = [
    "genome_id",
    "genome_name",
    "antibiotic",
    "resistant_phenotype",
    "y",
    "measurement",
    "measurement_sign",
    "measurement_value",
    "measurement_unit",
    "laboratory_typing_method",
    "testing_standard",
    "testing_standard_year",
    "evidence",
    "source",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data" / "processed")
    parser.add_argument("--min-presence", type=int, default=5)
    parser.add_argument("--min-absence", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    phenotypes = pd.read_csv(
        data_dir / "clean_phenotypes_lab_only_cipro_cefotaxime.csv",
        low_memory=False,
        dtype={"genome_id": str},
    )
    genes = pd.read_csv(
        data_dir / "clean_specialty_genes_amr_project_genomes.csv",
        low_memory=False,
        dtype=str,
    )

    pairs = genes[["genome_id", "predictor"]].dropna().drop_duplicates()
    matrix = pd.crosstab(pairs["genome_id"], pairs["predictor"]).clip(upper=1)
    matrix = matrix.reindex(sorted(matrix.columns), axis=1).reset_index()
    matrix.to_csv(data_dir / "gene_presence_matrix_all_project_genomes.csv", index=False)

    phenotype_summary = []
    prevalence_rows = []
    feature_cols = [c for c in matrix.columns if c.startswith("amr__")]

    for antibiotic in ANTIBIOTICS:
        subset = phenotypes.loc[phenotypes["antibiotic_norm"].eq(antibiotic)].copy()
        full = subset.merge(matrix, on="genome_id", how="inner")
        full = full[METADATA_COLUMNS + feature_cols]
        full_path = data_dir / f"dataset_{antibiotic}_gene_presence_full.csv"
        full.to_csv(full_path, index=False)

        sums = full[feature_cols].sum()
        keep = sums.index[
            (sums >= args.min_presence) & ((len(full) - sums) >= args.min_absence)
        ].tolist()
        ready = full[METADATA_COLUMNS + keep]
        ready.to_csv(data_dir / f"dataset_{antibiotic}_model_ready.csv", index=False)

        phenotype_summary.append(
            {
                "antibiotic": antibiotic,
                "n_genomes": len(ready),
                "n_resistant": int(ready["y"].sum()),
                "n_susceptible": int((ready["y"] == 0).sum()),
                "n_features_model_ready": len(keep),
            }
        )
        for feature in feature_cols:
            prevalence_rows.append(
                {
                    "antibiotic": antibiotic,
                    "feature": feature,
                    "n_present": int(sums[feature]),
                    "prevalence": float(sums[feature] / len(full)),
                    "kept_model_ready": feature in keep,
                }
            )

    pd.DataFrame(phenotype_summary).to_csv(data_dir / "phenotype_balance_summary.csv", index=False)
    pd.DataFrame(prevalence_rows).to_csv(data_dir / "feature_prevalence_summary.csv", index=False)
    print(f"Datos preparados en: {data_dir}")


if __name__ == "__main__":
    main()
