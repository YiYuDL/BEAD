"""Glide molecular docking jobs for BEAD."""

from __future__ import annotations

import argparse
import json
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
from rdkit import Chem, RDLogger

from bead.tools.ligand_input import materialize_ligand_csv


DEFAULT_RESULTS_ROOT = Path.home() / "bead_docking_results"
DEFAULT_SCHRODINGER = Path("/home/kemove/schrodinger2024-1")
REPORT_SCRIPT = Path(__file__).with_name("schrodinger_docking_report.py")
SMILES_COLUMNS = ("SMILES", "smiles", "canonical_smiles", "canonical_SMILES")


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _input_file(value: str | os.PathLike[str], label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} was not found: {path}")
    return path


def _executable(path: Path, label: str) -> Path:
    if not path.is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(f"{label} executable was not found: {path}")
    return path


def _run(command: list[str], cwd: Path, environment: dict[str, str]) -> None:
    printable = " ".join(command)
    print(f"$ {printable}", flush=True)
    subprocess.run(command, cwd=cwd, env=environment, check=True)


def _select_smiles_column(frame: pd.DataFrame, requested: str | None) -> str:
    if requested:
        if requested not in frame.columns:
            raise ValueError(
                f"SMILES column {requested!r} is absent; columns are {list(frame.columns)}"
            )
        return requested
    for candidate in SMILES_COLUMNS:
        if candidate in frame.columns:
            return candidate
    raise ValueError(
        "No SMILES column was found. Supply smiles_column explicitly; "
        f"columns are {list(frame.columns)}"
    )


def _prepare_inputs(spec: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = pd.read_csv(spec["csv_path"])
    smiles_column = _select_smiles_column(frame, spec.get("smiles_column"))
    valid_rows: list[dict[str, Any]] = []
    skipped_rows: list[dict[str, Any]] = []
    maximum = spec.get("max_molecules")
    RDLogger.DisableLog("rdApp.error")
    try:
        for source_index, value in frame[smiles_column].items():
            smiles = "" if pd.isna(value) else str(value).strip()
            reason = None
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            if not smiles:
                reason = "empty SMILES"
            elif molecule is None:
                reason = "RDKit could not parse SMILES"
            elif any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms()):
                reason = "SMILES contains an attachment/dummy atom"
            if reason:
                skipped_rows.append(
                    {
                        "source_index": int(source_index),
                        "source_csv_row": int(source_index) + 2,
                        "smiles": smiles,
                        "reason": reason,
                    }
                )
                continue
            if maximum is not None and len(valid_rows) >= int(maximum):
                continue
            valid_rows.append(
                {
                    "input_id": f"BEAD_{len(valid_rows) + 1:06d}",
                    "source_index": int(source_index),
                    "source_csv_row": int(source_index) + 2,
                    "smiles": smiles,
                }
            )
    finally:
        RDLogger.EnableLog("rdApp.error")
    if not valid_rows:
        raise ValueError("The input CSV contains no dockable SMILES")
    return (
        pd.DataFrame(
            valid_rows,
            columns=["input_id", "source_index", "source_csv_row", "smiles"],
        ),
        pd.DataFrame(
            skipped_rows,
            columns=["source_index", "source_csv_row", "smiles", "reason"],
        ),
    )


def _diagram_files(
    schrodinger: Path,
    output_dir: Path,
    records: list[dict[str, Any]],
    environment: dict[str, str],
) -> tuple[dict[str, str], dict[str, str]]:
    diagram_dir = output_dir / "interaction_diagrams"
    diagram_dir.mkdir(exist_ok=True)
    mapped: dict[str, str] = {}
    warnings: dict[str, str] = {}
    executable = _executable(
        schrodinger / "utilities" / "ligand_interaction_diagram",
        "interaction diagram",
    )
    # Rendering each selected complex independently avoids a Schrodinger 2024-1
    # crash when a rewritten multi-entry pose-viewer is supplied.
    for record in records:
        input_id = record["input_id"]
        final_path = diagram_dir / f"{input_id}.png"
        command = [
            str(executable),
            "-i",
            record["complex_pdb"],
            "-o",
            str(final_path),
            "-showLegend",
            "-showTitle",
        ]
        try:
            _run(command, output_dir, environment)
        except subprocess.CalledProcessError as exc:
            warnings[input_id] = f"interaction diagram failed: {exc}"
            continue
        if final_path.is_file():
            mapped[input_id] = str(final_path)
        else:
            warnings[input_id] = "interaction diagram command produced no PNG"
    return mapped, warnings


def _interaction_table(records: list[dict[str, Any]], output_path: Path) -> None:
    rows: list[dict[str, Any]] = []
    for record in records:
        input_id = record["input_id"]
        interactions = record.get("interactions", {})
        for item in interactions.get("hydrogen_bonds", []):
            rows.append({"input_id": input_id, "interaction_type": "hydrogen_bond", **item})
        for item in interactions.get("salt_bridges", []):
            rows.append({"input_id": input_id, "interaction_type": "salt_bridge", **item})
        for item in interactions.get("close_contact_residues", []):
            rows.append(
                {
                    "input_id": input_id,
                    "interaction_type": "close_contact_residue",
                    "protein_residue": item["residue"],
                    "minimum_distance_angstrom": item["minimum_distance_angstrom"],
                }
            )
    columns = [
        "input_id",
        "interaction_type",
        "protein_atom",
        "protein_residue",
        "ligand_atom",
        "distance_angstrom",
        "minimum_distance_angstrom",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(output_path, index=False)


def _execute(spec: dict[str, Any]) -> dict[str, Any]:
    job_dir = Path(spec["job_dir"])
    work_dir = job_dir / "work"
    output_dir = job_dir / "output"
    work_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    selected, skipped = _prepare_inputs(spec)
    selected.to_csv(output_dir / "selected_inputs.csv", index=False)
    skipped.to_csv(output_dir / "skipped_inputs.csv", index=False)

    ligand_input = selected[["smiles", "input_id"]].rename(
        columns={"input_id": "bead_input_id"}
    )
    ligand_input.to_csv(work_dir / "ligprep_input.csv", index=False)

    schrodinger = Path(spec["schrodinger_root"])
    ligprep = _executable(schrodinger / "ligprep", "LigPrep")
    glide = _executable(schrodinger / "glide", "Glide")
    schrodinger_run = _executable(schrodinger / "run", "Schrodinger run")
    environment = {
        **os.environ,
        "SCHRODINGER": str(schrodinger),
        "PATH": f"{schrodinger}:{schrodinger / 'utilities'}:{os.environ.get('PATH', '')}",
        "QT_QPA_PLATFORM": "offscreen",
    }

    ligprep_command = [
        str(ligprep),
        "-icsv",
        "ligprep_input.csv",
        "-osd",
        "ligprep.sdf",
        "-s",
        str(spec["max_stereoisomers"]),
        "-g",
        "-WAIT",
        "-HOST",
        "localhost:1",
    ]
    if spec["epik"]:
        ligprep_command.extend(["-epik", "-ph", str(spec["ph"]), "-pht", str(spec["ph_tolerance"])])
    else:
        ligprep_command.extend(["-i", "1"])
    _run(ligprep_command, work_dir, environment)
    ligand_sdf = work_dir / "ligprep.sdf"
    if not ligand_sdf.is_file() or ligand_sdf.stat().st_size == 0:
        raise RuntimeError("LigPrep did not produce any prepared ligands")

    glide_input = work_dir / "glide.in"
    glide_input.write_text(
        f"GRIDFILE {spec['grid_path']}\n"
        "LIGANDFILE ligprep.sdf\n"
        f"PRECISION {spec['precision']}\n"
        "POSE_OUTTYPE poseviewer\n"
        f"POSES_PER_LIG {spec['poses_per_ligand']}\n"
        "COMPRESS_POSES False\n"
        "WRITE_CSV True\n"
    )
    _run(
        [str(glide), "-WAIT", "-OVERWRITE", "-HOST", "localhost:1", "glide.in"],
        work_dir,
        environment,
    )
    poseviewer = work_dir / "glide_pv.mae"
    if not poseviewer.is_file():
        raise RuntimeError("Glide did not produce glide_pv.mae; inspect work/glide.log")

    report_json = output_dir / "interaction_details.json"
    _run(
        [
            str(schrodinger_run),
            "python3",
            str(REPORT_SCRIPT),
            "--poseviewer",
            str(poseviewer),
            "--output-dir",
            str(output_dir),
            "--report-json",
            str(report_json),
        ],
        work_dir,
        environment,
    )
    records = _read_json(report_json)
    if not records:
        raise RuntimeError("No docked poses were selected")
    diagrams, diagram_warnings = _diagram_files(
        schrodinger,
        output_dir,
        records,
        environment,
    )
    for record in records:
        record["interaction_diagram"] = diagrams.get(record["input_id"])
        warning = diagram_warnings.get(record["input_id"])
        if warning:
            record["interactions"].setdefault("analysis_warnings", []).append(warning)
    _write_json(report_json, records)
    _interaction_table(records, output_dir / "interaction_details.csv")

    flat_records = [
        {key: value for key, value in record.items() if key != "interactions"}
        for record in records
    ]
    results = selected.merge(pd.DataFrame(flat_records), on="input_id", how="left")
    source = pd.read_csv(spec["csv_path"]).reset_index(names="source_index")
    source = source.rename(
        columns={
            column: f"source_{column}"
            for column in source.columns
            if column != "source_index"
        }
    )
    results = results.merge(source, on="source_index", how="left")
    results["docking_status"] = results["glide_docking_score"].notna().map(
        {True: "docked", False: "no_pose"}
    )
    results = results.sort_values(
        ["glide_docking_score", "source_index"], na_position="last"
    )
    results_path = output_dir / "docking_results.csv"
    results.to_csv(results_path, index=False)
    if (work_dir / "glide.csv").is_file():
        shutil.copy2(work_dir / "glide.csv", output_dir / "glide_raw_results.csv")
    return {
        "input_rows": int(pd.read_csv(spec["csv_path"]).shape[0]),
        "selected_rows": int(len(selected)),
        "skipped_rows": int(len(skipped)),
        "docked_rows": int(results["glide_docking_score"].notna().sum()),
        "results_csv": str(results_path),
        "interaction_details_csv": str(output_dir / "interaction_details.csv"),
        "interaction_details_json": str(report_json),
        "output_dir": str(output_dir),
    }


def _artifacts(job_dir: Path) -> list[str]:
    keep = {".csv", ".json", ".log", ".in", ".pdb", ".sdf", ".mae", ".png"}
    return sorted(
        str(path)
        for path in job_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in keep
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
            "job_type": "mol_docking",
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
                "job_type": "mol_docking",
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
            "job_type": "mol_docking",
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


def mol_docking(
    csv_path: str | None,
    grid_path: str,
    smiles_column: str | None = None,
    smiles: list[str] | None = None,
    molecule_names: list[str] | None = None,
    precision: str = "SP",
    poses_per_ligand: int = 1,
    max_molecules: int | None = None,
    output_root: str | None = None,
    schrodinger_root: str | None = None,
    epik: bool = True,
    ph: float = 7.4,
    ph_tolerance: float = 1.0,
    max_stereoisomers: int = 1,
    background: bool = True,
) -> str:
    """Dock direct SMILES or a SMILES CSV against a prepared Glide grid."""
    grid_file = _input_file(grid_path, "Glide grid")
    precision = str(precision).strip().upper()
    if precision not in {"HTVS", "SP", "XP"}:
        raise ValueError("precision must be HTVS, SP, or XP")
    if int(poses_per_ligand) < 1:
        raise ValueError("poses_per_ligand must be at least 1")
    if max_molecules is not None and int(max_molecules) < 1:
        raise ValueError("max_molecules must be at least 1 when supplied")
    if int(max_stereoisomers) < 1:
        raise ValueError("max_stereoisomers must be at least 1")

    root = Path(output_root or DEFAULT_RESULTS_ROOT).expanduser().resolve()
    date_dir = root / datetime.now().astimezone().date().isoformat()
    job_id = f"mol-docking-{datetime.now().strftime('%H%M%S')}-{uuid.uuid4().hex[:8]}"
    job_dir = date_dir / job_id
    inputs_dir = job_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=False)
    local_csv, input_mode = materialize_ligand_csv(
        inputs_dir,
        csv_path=csv_path,
        smiles=smiles,
        molecule_names=molecule_names,
    )
    local_grid = inputs_dir / grid_file.name
    shutil.copy2(grid_file, local_grid)
    spec = {
        "job_id": job_id,
        "job_dir": str(job_dir),
        "csv_path": str(local_csv),
        "input_mode": input_mode,
        "grid_path": str(local_grid),
        "smiles_column": "smiles" if input_mode == "direct_smiles" else smiles_column,
        "precision": precision,
        "poses_per_ligand": int(poses_per_ligand),
        "max_molecules": int(max_molecules) if max_molecules is not None else None,
        "schrodinger_root": str(
            Path(schrodinger_root or DEFAULT_SCHRODINGER).expanduser().resolve()
        ),
        "epik": bool(epik),
        "ph": float(ph),
        "ph_tolerance": float(ph_tolerance),
        "max_stereoisomers": int(max_stereoisomers),
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    }
    spec_path = job_dir / "job.json"
    _write_json(spec_path, spec)
    status_path = job_dir / "status.json"
    _write_json(
        status_path,
        {
            "status": "queued",
            "job_type": "mol_docking",
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
        },
    )
    if not background:
        return_code = _worker(spec_path)
        if return_code:
            status = _read_json(status_path)
            raise RuntimeError(status.get("error", "Docking failed"))
        return json.dumps(_read_json(status_path), indent=2, ensure_ascii=False)

    log_handle = (job_dir / "job.log").open("a")
    try:
        process = subprocess.Popen(
            [sys.executable, "-m", "bead.tools.docking", "--worker", str(spec_path)],
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


def docking_job_status(job_path: str, log_lines: int = 30) -> str:
    """Return status and recent log lines for one mol_docking job."""
    path = Path(job_path).expanduser().resolve()
    status_path = path if path.name == "status.json" else path / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"Docking status file was not found: {status_path}")
    status = _read_json(status_path)
    log_path = Path(status.get("log", status_path.parent / "job.log"))
    if log_path.is_file() and log_lines > 0:
        status["log_tail"] = "\n".join(
            log_path.read_text(errors="replace").splitlines()[-int(log_lines):]
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
