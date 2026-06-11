#!/usr/bin/env python3
"""Run the exact optimized three-model analysis with repository-relative paths."""

from __future__ import annotations

import importlib.util
import os
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
ORIGINAL_SCRIPT = Path(__file__).with_name("run_optimized_models_original.py")


def load_original():
    spec = importlib.util.spec_from_file_location("optimized_original", ORIGINAL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"No se pudo cargar {ORIGINAL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    module = load_original()
    module.DATASETS = {
        "ciprofloxacin": REPO_ROOT / "data" / "processed" / "dataset_ciprofloxacin_model_ready.csv",
        "cefotaxime": REPO_ROOT / "data" / "processed" / "dataset_cefotaxime_model_ready.csv",
    }
    default_outdir = REPO_ROOT / "results" / "optimized_reproduced"
    if "--outdir" not in sys.argv:
        sys.argv.extend(["--outdir", str(default_outdir)])

    os.chdir(REPO_ROOT)
    args = module.parse_args()
    module.run_analysis(args)

    generated_zip = REPO_ROOT / f"{Path(args.outdir).name}.zip"
    destination_zip = Path(args.outdir).parent / f"{Path(args.outdir).name}.zip"
    if generated_zip.exists() and generated_zip.resolve() != destination_zip.resolve():
        destination_zip.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(generated_zip), str(destination_zip))


if __name__ == "__main__":
    main()
