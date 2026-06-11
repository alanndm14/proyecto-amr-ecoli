#!/usr/bin/env python
"""Valida la estructura mínima, codificación y tamaño del repositorio."""

from __future__ import annotations

import ast
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
MAX_FILE_SIZE = 100 * 1024 * 1024

REQUIRED_FILES = [
    "README.md",
    "requirements.txt",
    ".gitignore",
    "CITATION.cff",
    "data/processed/clean_phenotypes_lab_only_cipro_cefotaxime.csv",
    "data/processed/clean_specialty_genes_amr_project_genomes.csv",
    "data/processed/dataset_ciprofloxacin_model_ready.csv",
    "data/processed/dataset_cefotaxime_model_ready.csv",
    "scripts/prepare_model_data.py",
    "scripts/run_all_analyses.py",
    "scripts/run_optimized_models.py",
    "scripts/optimized_models_core.py",
    "scripts/run_mlp_multiseed.py",
    "scripts/run_robustness_by_lab_method.py",
    "results/optimized/metricas_holdout_todos_los_modelos.csv",
    "results/optimized/modelos_seleccionados_resumen.csv",
    "results/mlp/comparacion_holdout_modelos_previos_mas_mlp.csv",
    "results/mlp/mlp_holdout_resumen_semillas.csv",
    "results/robustness_method/robustez_leave_method_out.csv",
    "results/tablas_finales/tabla_1_composicion_datasets.csv",
    "results/tablas_finales/tabla_2_comparacion_modelos_incluye_mlp.csv",
    "results/tablas_finales/tabla_3_robustez_leave_method_out.csv",
    "results/tablas_finales/tabla_4_top_genes_estables.csv",
    "results/tablas_finales/tabla_anexo_mlp_resumen_semillas.csv",
    "report/informe_final.pdf",
]

EXPECTED_FIGURES = {
    "figura_1_flujo_general.png",
    "figura_2_comparacion_balanced_accuracy_modelos.png",
    "figura_3_falsos_susceptibles_modelos.png",
    "figura_4_matrices_confusion_modelos_seleccionados.png",
    "figura_5_genes_estables_ciprofloxacin.png",
    "figura_6_genes_estables_cefotaxime.png",
    "figura_7_robustez_leave_method_out.png",
}

TEXT_SUFFIXES = {".md", ".txt", ".py", ".yml", ".yaml", ".cff"}
MOJIBAKE_MARKERS = (
    "\ufffd",
    "\u00c3\u00a1",
    "\u00c3\u00a9",
    "\u00c3\u00ad",
    "\u00c3\u00b3",
    "\u00c3\u00ba",
    "\u00c3\u00b1",
    "\u00c2",
)


def tracked_candidates() -> list[Path]:
    ignored_parts = {".git", "__pycache__", ".ipynb_checkpoints"}
    return [
        path
        for path in ROOT.rglob("*")
        if path.is_file() and not any(part in ignored_parts for part in path.parts)
    ]


def main() -> int:
    errors: list[str] = []
    files = tracked_candidates()

    for relative_path in REQUIRED_FILES:
        if not (ROOT / relative_path).is_file():
            errors.append(f"Falta archivo requerido: {relative_path}")

    figures_dir = ROOT / "figures"
    actual_figures = (
        {path.name for path in figures_dir.glob("*.png")} if figures_dir.is_dir() else set()
    )
    if actual_figures != EXPECTED_FIGURES:
        missing = sorted(EXPECTED_FIGURES - actual_figures)
        extra = sorted(actual_figures - EXPECTED_FIGURES)
        if missing:
            errors.append(f"Faltan figuras finales: {', '.join(missing)}")
        if extra:
            errors.append(f"Hay figuras PNG no esperadas: {', '.join(extra)}")

    for path in files:
        relative_path = path.relative_to(ROOT)
        if path.stat().st_size > MAX_FILE_SIZE:
            errors.append(f"Archivo mayor a 100 MB: {relative_path}")

        if path.suffix.lower() in TEXT_SUFFIXES:
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                errors.append(f"Archivo de texto no codificado en UTF-8: {relative_path}")
                continue
            if any(marker in text for marker in MOJIBAKE_MARKERS):
                errors.append(f"Texto posiblemente corrupto: {relative_path}")

        if path.suffix.lower() == ".py":
            try:
                ast.parse(path.read_text(encoding="utf-8"), filename=str(relative_path))
            except SyntaxError as exc:
                errors.append(f"Error de sintaxis en {relative_path}: {exc}")

    if errors:
        print("VALIDACION FALLIDA")
        for error in errors:
            print(f"- {error}")
        return 1

    print(
        "VALIDACION CORRECTA: estructura reproducible, 7 figuras finales, "
        "archivos UTF-8 y ningun archivo mayor a 100 MB."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
