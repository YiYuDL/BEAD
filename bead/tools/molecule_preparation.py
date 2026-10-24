"""Small-molecule standardization and 3D preparation for BEAD."""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import traceback
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize


DEFAULT_RESULTS_ROOT = Path.home() / "bead_molecule_preparation_results"
STRUCTURE_COLUMNS = (
    "SMILES",
    "smiles",
    "canonical_smiles",
    "canonical_SMILES",
    "InChI",
    "inchi",
)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _select_column(frame: pd.DataFrame, requested: str | None) -> str:
    if requested:
        if requested not in frame.columns:
            raise ValueError(
                f"Structure column {requested!r} is absent; columns are {list(frame.columns)}"
            )
        return requested
    for candidate in STRUCTURE_COLUMNS:
        if candidate in frame.columns:
            return candidate
    raise ValueError(
        "No SMILES/InChI column was found. Supply structure_column explicitly; "
        f"columns are {list(frame.columns)}"
    )


def _load_inputs(spec: dict[str, Any]) -> list[dict[str, Any]]:
    if spec["input_kind"] == "direct":
        return [{"source_index": 0, "structure": spec["molecular_input"], "name": "molecule_1"}]

    path = Path(spec["molecular_input"])
    suffix = path.suffix.lower()
    rows: list[dict[str, Any]] = []
    if suffix == ".csv":
        frame = pd.read_csv(path)
        structure_column = _select_column(frame, spec.get("structure_column"))
        name_column = spec.get("name_column")
        if name_column and name_column not in frame.columns:
            raise ValueError(
                f"Name column {name_column!r} is absent; columns are {list(frame.columns)}"
            )
        for index, value in frame[structure_column].items():
            name = frame.at[index, name_column] if name_column else f"molecule_{index + 1}"
            rows.append(
                {
                    "source_index": int(index),
                    "structure": "" if pd.isna(value) else str(value).strip(),
                    "name": f"molecule_{index + 1}" if pd.isna(name) else str(name),
                }
            )
    elif suffix in {".smi", ".smiles", ".txt"}:
        for index, line in enumerate(path.read_text(errors="replace").splitlines()):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            rows.append(
                {
                    "source_index": index,
                    "structure": fields[0],
                    "name": " ".join(fields[1:]) if len(fields) > 1 else f"molecule_{index + 1}",
                }
            )
    else:
        raise ValueError("Input files must be CSV, SMI, SMILES, or TXT files")
    maximum = spec.get("max_molecules")
    return rows if maximum is None else rows[: int(maximum)]


def _parse_structure(value: str) -> Chem.Mol | None:
    value = str(value).strip()
    if not value:
        return None
    if value.startswith("InChI="):
        return Chem.MolFromInchi(value)
    return Chem.MolFromSmiles(value)


def _standardize(molecule: Chem.Mol, largest_fragment: bool) -> Chem.Mol:
    molecule = rdMolStandardize.Cleanup(molecule)
    if largest_fragment:
        molecule = rdMolStandardize.FragmentParent(molecule)
    Chem.SanitizeMol(molecule)
    return molecule


def _optimize_conformers(
    molecule: Chem.Mol, conformer_ids: list[int], threads: int, max_iterations: int
) -> tuple[str, list[tuple[int, float]]]:
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        results = AllChem.MMFFOptimizeMoleculeConfs(
            molecule,
            numThreads=threads,
            maxIters=max_iterations,
            mmffVariant="MMFF94",
        )
        return "MMFF94", [(int(status), float(energy)) for status, energy in results]
    if AllChem.UFFHasAllMoleculeParams(molecule):
        results = AllChem.UFFOptimizeMoleculeConfs(
            molecule,
            numThreads=threads,
            maxIters=max_iterations,
        )
        return "UFF", [(int(status), float(energy)) for status, energy in results]
    raise ValueError("Neither MMFF94 nor UFF has parameters for this molecule")


def _assign_gasteiger_charges(molecule: Chem.Mol) -> list[float]:
    AllChem.ComputeGasteigerCharges(molecule)
    charges: list[float] = []
    for atom in molecule.GetAtoms():
        value = float(atom.GetProp("_GasteigerCharge"))
        if not math.isfinite(value):
            raise ValueError("Gasteiger charge assignment produced a non-finite value")
        atom.SetDoubleProp("PartialCharge", value)
        charges.append(value)
    Chem.CreateAtomDoublePropertyList(molecule, "PartialCharge")
    return charges


def _single_conformer(molecule: Chem.Mol, conformer_id: int) -> Chem.Mol:
    result = Chem.Mol(molecule)
    conformer = Chem.Conformer(molecule.GetConformer(conformer_id))
    result.RemoveAllConformers()
    result.AddConformer(conformer, assignId=True)
    return result


def _write_pdbqt(molecule: Chem.Mol, destination: Path) -> None:
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    preparator = MoleculePreparation(charge_model="gasteiger")
    setups = preparator.prepare(molecule)
    if not setups:
        raise RuntimeError("Meeko produced no ligand setup")
    pdbqt, success, message = PDBQTWriterLegacy.write_string(
        setups[0], add_index_map=True, bad_charge_ok=False
    )
    if not success:
        raise RuntimeError(f"Meeko could not write PDBQT: {message}")
    destination.write_text(pdbqt)


def _prepare_one(
    row: dict[str, Any], spec: dict[str, Any], output_dir: Path
) -> tuple[dict[str, Any], Chem.Mol, list[dict[str, Any]]]:
    input_id = f"BEAD_MOL_{row['source_index'] + 1:06d}"
    molecule = _parse_structure(row["structure"])
    if molecule is None:
        raise ValueError("RDKit could not parse the SMILES/InChI input")
    if any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms()):
        raise ValueError("Dummy or attachment-point atoms are not valid prepared molecules")
    if spec["standardize"]:
        molecule = _standardize(molecule, spec["largest_fragment"])
    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    formal_charge = int(Chem.GetFormalCharge(molecule))
    molecule = Chem.AddHs(molecule)
    molecule.SetProp("_Name", input_id)
    molecule.SetProp("Input_Name", str(row["name"]))
    molecule.SetProp("Input_Structure", row["structure"])
    molecule.SetProp("Canonical_SMILES", canonical_smiles)

    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = int(spec["seed"])
    parameters.numThreads = int(spec["rdkit_threads"])
    conformer_ids = list(
        AllChem.EmbedMultipleConfs(
            molecule,
            numConfs=int(spec["num_confs"]),
            params=parameters,
        )
    )
    if not conformer_ids:
        raise RuntimeError("ETKDGv3 could not generate a 3D conformation")
    force_field, optimization = _optimize_conformers(
        molecule,
        conformer_ids,
        int(spec["rdkit_threads"]),
        int(spec["max_iterations"]),
    )
    charges = _assign_gasteiger_charges(molecule)
    best_position = min(range(len(optimization)), key=lambda index: optimization[index][1])
    best_conformer_id = int(conformer_ids[best_position])
    best_status, best_energy = optimization[best_position]

    molecule_dir = output_dir / "molecules" / input_id
    molecule_dir.mkdir(parents=True, exist_ok=True)
    conformers_path = molecule_dir / "conformers.sdf"
    writer = Chem.SDWriter(str(conformers_path))
    try:
        for position, conformer_id in enumerate(conformer_ids):
            status, energy = optimization[position]
            molecule.SetIntProp("Conformer_ID", int(conformer_id))
            molecule.SetDoubleProp("Energy_kcal_mol", energy)
            molecule.SetIntProp("Optimization_Status", status)
            molecule.SetProp("Force_Field", force_field)
            molecule.SetProp("Partial_Charge_Method", "Gasteiger")
            writer.write(molecule, confId=int(conformer_id))
    finally:
        writer.close()

    best = _single_conformer(molecule, best_conformer_id)
    best.SetIntProp("Selected_Conformer_ID", best_conformer_id)
    best.SetDoubleProp("Energy_kcal_mol", best_energy)
    best.SetIntProp("Optimization_Status", best_status)
    best.SetProp("Force_Field", force_field)
    best.SetProp("Partial_Charge_Method", "Gasteiger")
    best_sdf = molecule_dir / "prepared.sdf"
    single_writer = Chem.SDWriter(str(best_sdf))
    single_writer.write(best)
    single_writer.close()
    pdb_path = molecule_dir / "prepared.pdb"
    Chem.MolToPDBFile(best, str(pdb_path))

    warnings: list[str] = []
    pdbqt_path: Path | None = None
    if spec["generate_pdbqt"]:
        try:
            pdbqt_path = molecule_dir / "prepared.pdbqt"
            _write_pdbqt(best, pdbqt_path)
        except Exception as exc:
            pdbqt_path = None
            warnings.append(f"PDBQT export failed: {type(exc).__name__}: {exc}")

    charge_rows = [
        {
            "input_id": input_id,
            "atom_index": atom.GetIdx() + 1,
            "element": atom.GetSymbol(),
            "formal_charge": atom.GetFormalCharge(),
            "gasteiger_partial_charge": charges[atom.GetIdx()],
        }
        for atom in best.GetAtoms()
    ]
    record = {
        "input_id": input_id,
        "source_index": int(row["source_index"]),
        "input_name": str(row["name"]),
        "input_structure": row["structure"],
        "canonical_smiles": canonical_smiles,
        "formal_charge": formal_charge,
        "atom_count_with_hydrogens": best.GetNumAtoms(),
        "conformers_generated": len(conformer_ids),
        "selected_conformer_id": best_conformer_id,
        "force_field": force_field,
        "optimization_status": best_status,
        "energy_kcal_mol": best_energy,
        "partial_charge_method": "Gasteiger",
        "partial_charge_sum": sum(charges),
        "prepared_sdf": str(best_sdf),
        "all_conformers_sdf": str(conformers_path),
        "prepared_pdb": str(pdb_path),
        "prepared_pdbqt": str(pdbqt_path) if pdbqt_path else None,
        "warnings": "; ".join(warnings),
        "status": "prepared",
        "error": None,
    }
    return record, best, charge_rows


def _execute(spec: dict[str, Any]) -> dict[str, Any]:
    job_dir = Path(spec["job_dir"])
    output_dir = job_dir / "output"
    output_dir.mkdir(exist_ok=True)
    rows = _load_inputs(spec)
    if not rows:
        raise ValueError("No molecular structures were found in the input")

    records: list[dict[str, Any]] = []
    charge_rows: list[dict[str, Any]] = []
    prepared_molecules: list[Chem.Mol] = []
    for row in rows:
        input_id = f"BEAD_MOL_{row['source_index'] + 1:06d}"
        try:
            record, molecule, atom_charges = _prepare_one(row, spec, output_dir)
            records.append(record)
            prepared_molecules.append(molecule)
            charge_rows.extend(atom_charges)
        except Exception as exc:
            records.append(
                {
                    "input_id": input_id,
                    "source_index": int(row["source_index"]),
                    "input_name": str(row["name"]),
                    "input_structure": row["structure"],
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
            traceback.print_exc()

    if not prepared_molecules:
        raise RuntimeError("No molecules were prepared successfully")
    combined_sdf = output_dir / "prepared_molecules.sdf"
    writer = Chem.SDWriter(str(combined_sdf))
    try:
        for molecule in prepared_molecules:
            writer.write(molecule)
    finally:
        writer.close()
    results_csv = output_dir / "preparation_results.csv"
    pd.DataFrame(records).to_csv(results_csv, index=False)
    charges_csv = output_dir / "partial_charges.csv"
    pd.DataFrame(charge_rows).to_csv(charges_csv, index=False)
    failed = sum(record["status"] == "failed" for record in records)
    report = {
        "method": {
            "standardization": "RDKit MolStandardize Cleanup and optional FragmentParent",
            "hydrogens": "RDKit explicit hydrogen addition",
            "conformer_generation": "ETKDGv3",
            "random_seed": spec["seed"],
            "requested_conformers_per_molecule": spec["num_confs"],
            "geometry_optimization": "MMFF94 with UFF fallback",
            "maximum_optimization_iterations": spec["max_iterations"],
            "partial_charges": "Gasteiger",
        },
        "input_count": len(rows),
        "prepared_count": len(prepared_molecules),
        "failed_count": failed,
        "outputs": {
            "prepared_molecules_sdf": str(combined_sdf),
            "preparation_results_csv": str(results_csv),
            "partial_charges_csv": str(charges_csv),
            "molecule_directory": str(output_dir / "molecules"),
        },
    }
    report_path = output_dir / "preparation_report.json"
    _write_json(report_path, report)
    return {
        "input_count": len(rows),
        "prepared_count": len(prepared_molecules),
        "failed_count": failed,
        "prepared_molecules_sdf": str(combined_sdf),
        "preparation_results_csv": str(results_csv),
        "partial_charges_csv": str(charges_csv),
        "preparation_report": str(report_path),
        "output_dir": str(output_dir),
    }


def _artifacts(job_dir: Path) -> list[str]:
    keep = {".csv", ".json", ".log", ".txt", ".smi", ".smiles", ".sdf", ".pdb", ".pdbqt"}
    return sorted(
        str(path) for path in job_dir.rglob("*") if path.is_file() and path.suffix.lower() in keep
    )


def _worker(spec_path: Path) -> int:
    spec = _read_json(spec_path)
    job_dir = Path(spec["job_dir"])
    status_path = job_dir / "status.json"
    started = datetime.now().astimezone().isoformat(timespec="seconds")
    _write_json(
        status_path,
        {
            "status": "running",
            "job_type": "Mol_preparation",
            "pid": os.getpid(),
            "started_at": started,
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
        },
    )
    try:
        summary = _execute(spec)
    except Exception as exc:
        _write_json(
            status_path,
            {
                "status": "failed",
                "job_type": "Mol_preparation",
                "pid": os.getpid(),
                "started_at": started,
                "finished_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "job_dir": str(job_dir),
                "log": str(job_dir / "job.log"),
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
                "artifacts": _artifacts(job_dir),
            },
        )
        traceback.print_exc()
        return 1
    _write_json(
        status_path,
        {
            "status": "completed",
            "job_type": "Mol_preparation",
            "pid": os.getpid(),
            "started_at": started,
            "finished_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
            **summary,
            "artifacts": _artifacts(job_dir),
        },
    )
    return 0


def molecule_preparation(
    molecular_input: str,
    structure_column: str | None = None,
    name_column: str | None = None,
    num_confs: int = 50,
    seed: int = 42,
    rdkit_threads: int = 0,
    max_iterations: int = 500,
    standardize: bool = True,
    largest_fragment: bool = True,
    generate_pdbqt: bool = True,
    max_molecules: int | None = None,
    output_root: str | None = None,
    background: bool = True,
) -> str:
    """Prepare direct SMILES/InChI input or a CSV/SMI collection as 3D molecules."""
    if int(num_confs) < 1:
        raise ValueError("num_confs must be at least 1")
    if int(rdkit_threads) < 0:
        raise ValueError("rdkit_threads cannot be negative")
    if int(max_iterations) < 1:
        raise ValueError("max_iterations must be at least 1")
    if max_molecules is not None and int(max_molecules) < 1:
        raise ValueError("max_molecules must be at least 1 when supplied")
    text = str(molecular_input).strip()
    if not text:
        raise ValueError("molecular_input cannot be empty")

    source = Path(text).expanduser()
    try:
        input_kind = "file" if source.is_file() else "direct"
    except OSError:
        input_kind = "direct"
    root = Path(output_root or DEFAULT_RESULTS_ROOT).expanduser().resolve()
    date_dir = root / datetime.now().astimezone().date().isoformat()
    job_id = f"mol-preparation-{datetime.now().strftime('%H%M%S')}-{uuid.uuid4().hex[:8]}"
    job_dir = date_dir / job_id
    inputs_dir = job_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=False)
    if input_kind == "file":
        local_input = inputs_dir / source.name
        shutil.copy2(source.resolve(), local_input)
        stored_input = str(local_input)
    else:
        (inputs_dir / "direct_input.txt").write_text(text + "\n")
        stored_input = text

    spec = {
        "job_id": job_id,
        "job_dir": str(job_dir),
        "molecular_input": stored_input,
        "input_kind": input_kind,
        "structure_column": structure_column,
        "name_column": name_column,
        "num_confs": int(num_confs),
        "seed": int(seed),
        "rdkit_threads": int(rdkit_threads),
        "max_iterations": int(max_iterations),
        "standardize": bool(standardize),
        "largest_fragment": bool(largest_fragment),
        "generate_pdbqt": bool(generate_pdbqt),
        "max_molecules": int(max_molecules) if max_molecules is not None else None,
        "python": sys.executable,
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    spec_path = job_dir / "job.json"
    _write_json(spec_path, spec)
    status_path = job_dir / "status.json"
    _write_json(
        status_path,
        {
            "status": "queued",
            "job_type": "Mol_preparation",
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
        },
    )
    if not background:
        return_code = _worker(spec_path)
        if return_code:
            status = _read_json(status_path)
            raise RuntimeError(status.get("error", "Molecule preparation failed"))
        return json.dumps(_read_json(status_path), indent=2, ensure_ascii=False)

    log_handle = (job_dir / "job.log").open("a")
    try:
        process = subprocess.Popen(
            [sys.executable, "-m", "bead.tools.molecule_preparation", "--worker", str(spec_path)],
            cwd=job_dir,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
        )
    finally:
        log_handle.close()
    response = _read_json(status_path)
    response.update({"pid": process.pid, "status_file": str(status_path)})
    return json.dumps(response, indent=2, ensure_ascii=False)


def molecule_preparation_job_status(job_path: str, log_lines: int = 30) -> str:
    """Return status and recent logs for a Mol_preparation job."""
    path = Path(job_path).expanduser().resolve()
    status_path = path if path.name == "status.json" else path / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"Molecule preparation status file was not found: {status_path}")
    status = _read_json(status_path)
    log_path = Path(status.get("log", status_path.parent / "job.log"))
    if log_path.is_file() and log_lines > 0:
        status["log_tail"] = "\n".join(
            log_path.read_text(errors="replace").splitlines()[-int(log_lines) :]
        )
    status["artifacts"] = _artifacts(status_path.parent)
    return json.dumps(status, indent=2, ensure_ascii=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    args = parser.parse_args()
    return _worker(args.worker.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
