"""Shared input handling for docking tools."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path
from typing import Sequence


def materialize_ligand_csv(
    inputs_dir: Path,
    csv_path: str | None,
    smiles: Sequence[str] | str | None,
    molecule_names: Sequence[str] | str | None,
) -> tuple[Path, str]:
    """Copy a supplied CSV or persist direct SMILES as a job-local CSV."""
    has_csv = bool(csv_path and str(csv_path).strip())
    has_smiles = smiles is not None
    if has_csv == has_smiles:
        raise ValueError("Supply exactly one of csv_path or smiles")

    if has_csv:
        if molecule_names is not None:
            raise ValueError("molecule_names can only be used with direct smiles input")
        source = Path(str(csv_path)).expanduser().resolve()
        if not source.is_file():
            raise FileNotFoundError(f"input CSV was not found: {source}")
        destination = inputs_dir / source.name
        shutil.copy2(source, destination)
        return destination, "csv"

    values = [smiles] if isinstance(smiles, str) else list(smiles or [])
    normalized_smiles = [str(value).strip() for value in values]
    if not normalized_smiles:
        raise ValueError("smiles must contain at least one molecule")

    if molecule_names is None:
        names = [f"Ligand_{index}" for index in range(1, len(normalized_smiles) + 1)]
    else:
        names = (
            [molecule_names]
            if isinstance(molecule_names, str)
            else [str(value).strip() for value in molecule_names]
        )
        if len(names) != len(normalized_smiles):
            raise ValueError("molecule_names and smiles must contain the same number of items")
        if any(not name for name in names):
            raise ValueError("molecule_names cannot contain empty values")

    destination = inputs_dir / "direct_smiles.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "smiles"])
        writer.writerows(zip(names, normalized_smiles))
    return destination, "direct_smiles"
