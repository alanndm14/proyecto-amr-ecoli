#!/usr/bin/env python3
"""Ejecuta la secuencia reproducible de análisis desde la raíz del repositorio."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent


def run(*args: str) -> None:
    command = [sys.executable, *args]
    print("\n>", " ".join(command))
    subprocess.run(command, check=True)


def main() -> None:
    run(
        str(SCRIPTS / "run_optimized_models.py"),
        "--n_iter_logreg", "5",
        "--n_iter_rf", "5",
        "--n_iter_xgb", "5",
        "--n_boot", "100",
        "--n_jobs_search", "1",
        "--verbose", "1",
    )
    run(str(SCRIPTS / "run_mlp_multiseed.py"))
    run(str(SCRIPTS / "run_robustness_by_lab_method.py"))


if __name__ == "__main__":
    main()
