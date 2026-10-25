"""Target preparation jobs for the BEAD agent.

The preferred backend is Schrödinger Protein Preparation Wizard.  When it is
unavailable, the workflow uses PDBFixer/OpenMM and Meeko so that BEAD can still
produce receptor files suitable for AutoDock Vina.
"""

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


DEFAULT_RESULTS_ROOT = Path.home() / "bead_target_preparation_results"
DEFAULT_SCHRODINGER_ROOT = Path("/home/kemove/schrodinger2024-1")
WATER_RESIDUES = {"HOH", "WAT", "H2O", "DOD"}
AUTO_BOX_EXCLUDED_RESIDUES = {
    # Common crystallization additives, buffers, and solvents. Single-atom ions
    # are rejected separately by the minimum-heavy-atom rule.
    "ACE",
    "ACT",
    "BME",
    "DMS",
    "EDO",
    "EOH",
    "FMT",
    "GOL",
    "IPA",
    "MES",
    "MPD",
    "PEG",
    "PG4",
    "TRS",
}


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _run(command: list[str], cwd: Path, environment: dict[str, str] | None = None) -> None:
    print("$ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, env=environment, check=True)


def _executable(path: Path, label: str) -> Path:
    if path.is_file() and os.access(path, os.X_OK):
        return path.resolve()
    raise FileNotFoundError(f"{label} was not found or is not executable: {path}")


def _environment_executable(name: str) -> Path:
    candidate = Path(sys.executable).resolve().with_name(name)
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate
    found = shutil.which(name)
    if found:
        return Path(found).resolve()
    raise FileNotFoundError(
        f"{name} was not found in the active BEAD environment ({sys.executable})"
    )


def _parse_keep_heterogens(value: str | None) -> set[str]:
    if value is None:
        return set()
    return {item.strip().upper() for item in str(value).split(",") if item.strip()}


def _pdb_element(line: str) -> str:
    """Return a normalized element symbol from a PDB atom record."""
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    atom_name = line[12:16].strip().upper()
    letters = "".join(character for character in atom_name if character.isalpha())
    if not letters:
        return ""
    if letters[:2] in {"BR", "CL"}:
        return letters[:2]
    return letters[0]


def _automatic_ligand_box(
    source: Path,
    ligand_resname: str | None,
    ligand_chain: str | None,
    ligand_resid: int | None,
    box_padding: float,
    minimum_sizes: list[float],
) -> dict[str, Any] | None:
    """Infer a docking box from a co-crystal ligand before PDB cleaning.

    Explicit selectors may identify a ligand with any number of heavy atoms. If
    no selector is supplied, the largest drug-like HETATM residue is selected;
    water, common crystallization additives, and residues with fewer than five
    heavy atoms or two carbon atoms are ignored.
    """
    residues: dict[tuple[str, str, int, str], dict[str, Any]] = {}
    model_seen = False
    first_model_finished = False

    for order, line in enumerate(source.read_text(errors="replace").splitlines()):
        record = line[:6].strip().upper()
        if record == "MODEL":
            if model_seen:
                first_model_finished = True
            model_seen = True
            continue
        if record == "ENDMDL":
            if model_seen:
                first_model_finished = True
            continue
        if first_model_finished or record != "HETATM":
            continue

        residue_name = line[17:20].strip().upper()
        if not residue_name or residue_name in WATER_RESIDUES:
            continue
        chain = line[21:22].strip()
        residue_text = line[22:26].strip()
        insertion_code = line[26:27].strip()
        try:
            residue_id = int(residue_text)
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        alternate_location = line[16:17].strip().upper()
        if alternate_location not in {"", "A", "1"}:
            continue
        element = _pdb_element(line)
        if element in {"H", "D"}:
            continue

        key = (residue_name, chain, residue_id, insertion_code)
        residue = residues.setdefault(
            key,
            {
                "resname": residue_name,
                "chain": chain,
                "resid": residue_id,
                "insertion_code": insertion_code,
                "atoms": {},
                "first_order": order,
            },
        )
        atom_name = line[12:16].strip().upper()
        residue["atoms"].setdefault(atom_name, (x, y, z, element))

    requested_resname = ligand_resname.strip().upper() if ligand_resname else None
    requested_chain = ligand_chain.strip() if ligand_chain else None
    selector_supplied = any(
        value is not None for value in (requested_resname, requested_chain, ligand_resid)
    )
    candidates: list[dict[str, Any]] = []
    for residue in residues.values():
        if requested_resname is not None and residue["resname"] != requested_resname:
            continue
        if requested_chain is not None and residue["chain"] != requested_chain:
            continue
        if ligand_resid is not None and residue["resid"] != ligand_resid:
            continue
        atoms = list(residue["atoms"].values())
        carbon_count = sum(atom[3] == "C" for atom in atoms)
        residue["heavy_atom_count"] = len(atoms)
        residue["carbon_atom_count"] = carbon_count
        if not selector_supplied:
            if residue["resname"] in AUTO_BOX_EXCLUDED_RESIDUES:
                continue
            if len(atoms) < 5 or carbon_count < 2:
                continue
        if atoms:
            candidates.append(residue)

    if not candidates:
        if selector_supplied:
            selector = ", ".join(
                part
                for part in (
                    f"resname={requested_resname}" if requested_resname else "",
                    f"chain={requested_chain}" if requested_chain is not None else "",
                    f"resid={ligand_resid}" if ligand_resid is not None else "",
                )
                if part
            )
            raise ValueError(f"No HETATM ligand matched the automatic-box selector ({selector})")
        return None

    selected = max(
        candidates,
        key=lambda residue: (
            residue["heavy_atom_count"],
            residue["carbon_atom_count"],
            -residue["first_order"],
        ),
    )
    coordinates = [atom[:3] for atom in selected["atoms"].values()]
    axes = list(zip(*coordinates))
    center = [sum(axis) / len(axis) for axis in axes]
    extents = [max(axis) - min(axis) for axis in axes]
    sizes = [
        max(float(minimum), float(extent) + 2.0 * float(box_padding))
        for minimum, extent in zip(minimum_sizes, extents)
    ]
    return {
        "selection_method": (
            "explicit_selector" if selector_supplied else "largest_drug_like_heterogen"
        ),
        "ligand": {
            "resname": selected["resname"],
            "chain": selected["chain"],
            "resid": selected["resid"],
            "insertion_code": selected["insertion_code"],
            "heavy_atom_count": selected["heavy_atom_count"],
            "carbon_atom_count": selected["carbon_atom_count"],
        },
        "center_method": "heavy_atom_centroid",
        "box_center": center,
        "ligand_extents": extents,
        "box_padding": float(box_padding),
        "minimum_box_size": [float(value) for value in minimum_sizes],
        "box_size": sizes,
        "eligible_candidate_count": len(candidates),
    }


def _clean_pdb(
    source: Path,
    destination: Path,
    remove_waters: bool,
    remove_heterogens: bool,
    keep_heterogens: set[str],
) -> dict[str, Any]:
    """Apply explicit, auditable solvent/heterogen filtering to a PDB file."""
    output: list[str] = []
    input_atoms = 0
    output_atoms = 0
    removed_water_atoms = 0
    removed_heterogen_atoms = 0
    removed_water_residues: set[tuple[str, str, str]] = set()
    removed_heterogen_residues: set[tuple[str, str, str]] = set()
    retained_heterogens: set[str] = set()
    model_seen = False
    first_model_finished = False

    for line in source.read_text(errors="replace").splitlines():
        record = line[:6].strip().upper()
        if record == "MODEL":
            if model_seen:
                first_model_finished = True
            model_seen = True
            if not first_model_finished:
                output.append(line)
            continue
        if record == "ENDMDL":
            if model_seen and not first_model_finished:
                output.append(line)
                first_model_finished = True
            continue
        if first_model_finished and record in {"ATOM", "HETATM", "ANISOU", "TER"}:
            continue
        if record == "ANISOU":
            continue
        if record == "CONECT":
            # Connectivity can refer to atoms removed by the filter. Both downstream
            # preparation backends reconstruct the relevant chemistry.
            continue
        if record not in {"ATOM", "HETATM"}:
            if record != "END":
                output.append(line)
            continue

        input_atoms += 1
        residue_name = line[17:20].strip().upper()
        residue_key = (line[21:22].strip(), line[22:27].strip(), residue_name)
        is_water = residue_name in WATER_RESIDUES
        if remove_waters and is_water:
            removed_water_atoms += 1
            removed_water_residues.add(residue_key)
            continue
        if (
            record == "HETATM"
            and not is_water
            and remove_heterogens
            and residue_name not in keep_heterogens
        ):
            removed_heterogen_atoms += 1
            removed_heterogen_residues.add(residue_key)
            continue
        if record == "HETATM":
            retained_heterogens.add(residue_name)
        output.append(line)
        output_atoms += 1

    if output_atoms == 0:
        raise ValueError("No receptor atoms remain after water/heterogen filtering")
    output.append("END")
    destination.write_text("\n".join(output) + "\n")
    return {
        "input_atom_count": input_atoms,
        "cleaned_atom_count": output_atoms,
        "removed_water_atom_count": removed_water_atoms,
        "removed_water_residue_count": len(removed_water_residues),
        "removed_heterogen_atom_count": removed_heterogen_atoms,
        "removed_heterogen_residue_count": len(removed_heterogen_residues),
        "retained_heterogen_residue_names": sorted(retained_heterogens),
    }


def _schrodinger_root(spec: dict[str, Any]) -> Path | None:
    prepwizard_on_path = shutil.which("prepwizard")
    path_root = (
        Path(prepwizard_on_path).resolve().parent.parent if prepwizard_on_path else None
    )
    choices = [
        spec.get("schrodinger_root"),
        os.environ.get("SCHRODINGER"),
        str(path_root) if path_root else None,
        str(DEFAULT_SCHRODINGER_ROOT),
    ]
    for value in choices:
        if not value:
            continue
        root = Path(value).expanduser().resolve()
        if (root / "utilities" / "prepwizard").is_file():
            return root
    return None


def _prepare_with_schrodinger(
    spec: dict[str, Any], cleaned: Path, work_dir: Path, output_dir: Path, root: Path
) -> tuple[Path, Path]:
    prepwizard = _executable(root / "utilities" / "prepwizard", "Protein Preparation Wizard")
    structconvert = _executable(root / "utilities" / "structconvert", "structconvert")
    maestro_work = work_dir / "prepared_target.maegz"
    maestro_path = output_dir / "prepared_target.maegz"
    pdb_path = output_dir / "prepared_target.pdb"
    command = [
        str(prepwizard),
        "-WAIT",
        "-rehtreat",
        "-disulfides",
        "-epik_pH",
        str(spec["ph"]),
        "-epik_pHt",
        str(spec["ph_tolerance"]),
        "-propka_pH",
        str(spec["ph"]),
    ]
    if spec["fill_missing_sidechains"]:
        command.append("-fillsidechains")
    if spec["fill_missing_loops"]:
        command.append("-fillloops")
    if spec["minimize"]:
        command.extend(["-rmsd", str(spec["minimize_rmsd"])])
    else:
        command.append("-noimpref")
    command.extend([str(cleaned), str(maestro_work)])
    environment = {
        **os.environ,
        "SCHRODINGER": str(root),
        "PATH": f"{root}:{root / 'utilities'}:{os.environ.get('PATH', '')}",
    }
    _run(command, work_dir, environment)
    if not maestro_work.is_file() or maestro_work.stat().st_size == 0:
        raise RuntimeError("Protein Preparation Wizard did not produce prepared_target.maegz")
    shutil.copy2(maestro_work, maestro_path)
    _run([str(structconvert), str(maestro_path), str(pdb_path)], work_dir, environment)
    if not pdb_path.is_file() or pdb_path.stat().st_size == 0:
        raise RuntimeError("structconvert did not produce prepared_target.pdb")
    return pdb_path, maestro_path


def _prepare_with_pdbfixer(spec: dict[str, Any], cleaned: Path, output_dir: Path) -> Path:
    try:
        from openmm import app
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise RuntimeError(
            "The AutoDock fallback requires pdbfixer and openmm in the BEAD environment"
        ) from exc

    fixer = PDBFixer(filename=str(cleaned))
    fixer.findMissingResidues()
    if not spec["fill_missing_loops"]:
        fixer.missingResidues = {}
    if spec["fill_missing_sidechains"] or spec["fill_missing_loops"]:
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
    fixer.addMissingHydrogens(float(spec["ph"]))
    pdb_path = output_dir / "prepared_target.pdb"
    with pdb_path.open("w") as handle:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    if not pdb_path.is_file() or pdb_path.stat().st_size == 0:
        raise RuntimeError("PDBFixer did not produce prepared_target.pdb")
    return pdb_path


def _prepare_pdbqt(prepared_pdb: Path, work_dir: Path, output_dir: Path) -> Path:
    prepare = _environment_executable("mk_prepare_receptor.py")
    basename = work_dir / "prepared_target"
    _run(
        [str(prepare), "--read_pdb", str(prepared_pdb), "-o", str(basename), "-p", "-x"],
        work_dir,
    )
    generated = basename.with_suffix(".pdbqt")
    if not generated.is_file() or generated.stat().st_size == 0:
        raise RuntimeError("Meeko did not produce prepared_target.pdbqt")
    destination = output_dir / "prepared_target.pdbqt"
    shutil.copy2(generated, destination)
    return destination


def _write_vina_box(spec: dict[str, Any], output_dir: Path) -> Path | None:
    center = spec.get("box_center")
    if center is None:
        return None
    size = spec["box_size"]
    path = output_dir / "vina_box.txt"
    path.write_text(
        "\n".join(
            [
                f"center_x = {center[0]}",
                f"center_y = {center[1]}",
                f"center_z = {center[2]}",
                f"size_x = {size[0]}",
                f"size_y = {size[1]}",
                f"size_z = {size[2]}",
            ]
        )
        + "\n"
    )
    return path


def _prepare_glide_grid(
    spec: dict[str, Any], maestro_path: Path, work_dir: Path, output_dir: Path, root: Path
) -> Path | None:
    center = spec.get("box_center")
    if center is None:
        return None
    glide = _executable(root / "glide", "Glide")
    size = spec["box_size"]
    input_path = work_dir / "glide-grid.in"
    input_path.write_text(
        "\n".join(
            [
                f"RECEP_FILE {maestro_path}",
                f"GRID_CENTER {center[0]}, {center[1]}, {center[2]}",
                "INNERBOX 10, 10, 10",
                f"OUTERBOX {size[0]}, {size[1]}, {size[2]}",
            ]
        )
        + "\n"
    )
    environment = {
        **os.environ,
        "SCHRODINGER": str(root),
        "PATH": f"{root}:{root / 'utilities'}:{os.environ.get('PATH', '')}",
    }
    _run([str(glide), "-WAIT", "-OVERWRITE", str(input_path)], work_dir, environment)
    generated = work_dir / "glide-grid.zip"
    if not generated.is_file() or generated.stat().st_size == 0:
        candidates = sorted(work_dir.glob("*.zip"), key=lambda path: path.stat().st_mtime)
        if not candidates:
            raise RuntimeError("Glide did not produce a grid archive")
        generated = candidates[-1]
    destination = output_dir / "glide-grid.zip"
    shutil.copy2(generated, destination)
    shutil.copy2(input_path, output_dir / input_path.name)
    return destination


def _count_pdb_atoms(path: Path) -> dict[str, int]:
    atom_count = 0
    hydrogen_count = 0
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom_count += 1
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            atom_name = line[12:16].strip().upper()
            if element == "H" or (not element and atom_name.startswith("H")):
                hydrogen_count += 1
    return {"prepared_atom_count": atom_count, "prepared_hydrogen_count": hydrogen_count}


def _execute(spec: dict[str, Any]) -> dict[str, Any]:
    job_dir = Path(spec["job_dir"])
    work_dir = job_dir / "work"
    output_dir = job_dir / "output"
    work_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    warnings: list[str] = []
    automatic_box: dict[str, Any] | None = None
    box_source = "manual" if spec.get("box_center") is not None else "none"
    if spec.get("auto_box", False) and spec.get("box_center") is None:
        automatic_box = _automatic_ligand_box(
            Path(spec["pdb_path"]),
            spec.get("ligand_resname"),
            spec.get("ligand_chain"),
            spec.get("ligand_resid"),
            float(spec.get("box_padding", 4.0)),
            [float(value) for value in spec["box_size"]],
        )
        if automatic_box is None:
            warnings.append(
                "Automatic docking-box detection found no drug-like HETATM residue; "
                "the receptor was prepared without a docking box. Supply a ligand selector "
                "or manual center coordinates before docking."
            )
        else:
            spec["box_center"] = automatic_box["box_center"]
            spec["box_size"] = automatic_box["box_size"]
            box_source = "co_crystal_ligand"
            ligand = automatic_box["ligand"]
            print(
                "Automatic docking box: "
                f"{ligand['resname']} {ligand['chain'] or '-'}{ligand['resid']} "
                f"center={spec['box_center']} size={spec['box_size']}",
                flush=True,
            )
    cleaned = output_dir / "cleaned_input.pdb"
    cleaning = _clean_pdb(
        Path(spec["pdb_path"]),
        cleaned,
        spec["remove_waters"],
        spec["remove_heterogens"],
        set(spec["keep_heterogens"]),
    )

    requested_backend = spec["backend"]
    root = _schrodinger_root(spec)
    maestro_path: Path | None = None
    if requested_backend == "schrodinger" and root is None:
        raise FileNotFoundError("Schrödinger Protein Preparation Wizard was not found")

    if requested_backend in {"auto", "schrodinger"} and root is not None:
        try:
            prepared_pdb, maestro_path = _prepare_with_schrodinger(
                spec, cleaned, work_dir, output_dir, root
            )
            backend_used = "schrodinger"
        except Exception as exc:
            if requested_backend == "schrodinger" or not spec["allow_fallback"]:
                raise
            warnings.append(
                f"Schrodinger preparation failed; used AutoDock fallback: {type(exc).__name__}: {exc}"
            )
            prepared_pdb = _prepare_with_pdbfixer(spec, cleaned, output_dir)
            backend_used = "autodock"
            maestro_path = None
    else:
        prepared_pdb = _prepare_with_pdbfixer(spec, cleaned, output_dir)
        backend_used = "autodock"

    if backend_used == "autodock" and spec["minimize"]:
        warnings.append(
            "The open-source fallback repairs atoms and protonates the receptor but does not "
            "perform Schrodinger restrained minimization."
        )

    pdbqt_path = _prepare_pdbqt(prepared_pdb, work_dir, output_dir)
    vina_box = _write_vina_box(spec, output_dir)
    glide_grid = None
    if backend_used == "schrodinger" and maestro_path is not None and root is not None:
        glide_grid = _prepare_glide_grid(spec, maestro_path, work_dir, output_dir, root)

    report = {
        "requested_backend": requested_backend,
        "backend_used": backend_used,
        "input_pdb": spec["pdb_path"],
        "ph": spec["ph"],
        "remove_waters": spec["remove_waters"],
        "remove_heterogens": spec["remove_heterogens"],
        "keep_heterogens": spec["keep_heterogens"],
        "fill_missing_sidechains": spec["fill_missing_sidechains"],
        "fill_missing_loops": spec["fill_missing_loops"],
        "minimize": spec["minimize"] if backend_used == "schrodinger" else False,
        "auto_box_requested": bool(spec.get("auto_box", False)),
        "box_source": box_source,
        "automatic_box": automatic_box,
        "box_center": spec.get("box_center"),
        "box_size": spec.get("box_size") if spec.get("box_center") else None,
        **cleaning,
        **_count_pdb_atoms(prepared_pdb),
        "warnings": warnings,
        "outputs": {
            "cleaned_pdb": str(cleaned),
            "prepared_pdb": str(prepared_pdb),
            "prepared_pdbqt": str(pdbqt_path),
            "prepared_maegz": str(maestro_path) if maestro_path else None,
            "vina_box": str(vina_box) if vina_box else None,
            "glide_grid": str(glide_grid) if glide_grid else None,
        },
    }
    report_path = output_dir / "preparation_report.json"
    _write_json(report_path, report)
    return {
        "backend_used": backend_used,
        "box_source": box_source,
        "box_center": spec.get("box_center"),
        "box_size": spec.get("box_size") if spec.get("box_center") else None,
        "box_ligand": automatic_box["ligand"] if automatic_box else None,
        "prepared_pdb": str(prepared_pdb),
        "prepared_pdbqt": str(pdbqt_path),
        "prepared_maegz": str(maestro_path) if maestro_path else None,
        "vina_box": str(vina_box) if vina_box else None,
        "glide_grid": str(glide_grid) if glide_grid else None,
        "preparation_report": str(report_path),
        "output_dir": str(output_dir),
        "warnings": warnings,
    }


def _artifacts(job_dir: Path) -> list[str]:
    keep = {".json", ".log", ".txt", ".pdb", ".pdbqt", ".mae", ".maegz", ".zip", ".in"}
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
            "job_type": "Target_preparation",
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
                "job_type": "Target_preparation",
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
            "job_type": "Target_preparation",
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


def target_preparation(
    pdb_path: str,
    backend: str = "auto",
    ph: float = 7.4,
    ph_tolerance: float = 1.0,
    remove_waters: bool = True,
    remove_heterogens: bool = True,
    keep_heterogens: str | None = None,
    fill_missing_sidechains: bool = True,
    fill_missing_loops: bool = False,
    minimize: bool = True,
    minimize_rmsd: float = 0.3,
    center_x: float | None = None,
    center_y: float | None = None,
    center_z: float | None = None,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    auto_box: bool = True,
    ligand_resname: str | None = None,
    ligand_chain: str | None = None,
    ligand_resid: int | None = None,
    box_padding: float = 4.0,
    output_root: str | None = None,
    schrodinger_root: str | None = None,
    allow_fallback: bool = True,
    background: bool = True,
) -> str:
    """Clean and prepare a raw PDB receptor for Glide and AutoDock Vina."""
    source = Path(pdb_path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Input PDB was not found: {source}")
    if source.suffix.lower() not in {".pdb", ".ent"}:
        raise ValueError("pdb_path must point to a PDB or ENT file")
    backend = str(backend).strip().lower()
    if backend not in {"auto", "schrodinger", "autodock"}:
        raise ValueError("backend must be 'auto', 'schrodinger', or 'autodock'")
    if not 0.0 <= float(ph) <= 14.0:
        raise ValueError("ph must be between 0 and 14")
    if float(ph_tolerance) < 0 or float(minimize_rmsd) < 0:
        raise ValueError("ph_tolerance and minimize_rmsd cannot be negative")
    if float(box_padding) < 0:
        raise ValueError("box_padding cannot be negative")
    center_values = [center_x, center_y, center_z]
    if any(value is not None for value in center_values) and not all(
        value is not None for value in center_values
    ):
        raise ValueError("center_x, center_y, and center_z must be supplied together")
    sizes = [float(size_x), float(size_y), float(size_z)]
    if any(value <= 0 for value in sizes):
        raise ValueError("All box dimensions must be greater than zero")
    if not auto_box and any(
        value is not None for value in (ligand_resname, ligand_chain, ligand_resid)
    ):
        raise ValueError("Ligand selectors require auto_box=True")

    root = Path(output_root or DEFAULT_RESULTS_ROOT).expanduser().resolve()
    date_dir = root / datetime.now().astimezone().date().isoformat()
    job_id = f"target-preparation-{datetime.now().strftime('%H%M%S')}-{uuid.uuid4().hex[:8]}"
    job_dir = date_dir / job_id
    inputs_dir = job_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=False)
    local_pdb = inputs_dir / source.name
    shutil.copy2(source, local_pdb)
    spec = {
        "job_id": job_id,
        "job_dir": str(job_dir),
        "pdb_path": str(local_pdb),
        "backend": backend,
        "ph": float(ph),
        "ph_tolerance": float(ph_tolerance),
        "remove_waters": bool(remove_waters),
        "remove_heterogens": bool(remove_heterogens),
        "keep_heterogens": sorted(_parse_keep_heterogens(keep_heterogens)),
        "fill_missing_sidechains": bool(fill_missing_sidechains),
        "fill_missing_loops": bool(fill_missing_loops),
        "minimize": bool(minimize),
        "minimize_rmsd": float(minimize_rmsd),
        "box_center": [float(value) for value in center_values] if center_x is not None else None,
        "box_size": sizes,
        "auto_box": bool(auto_box),
        "ligand_resname": ligand_resname.strip().upper() if ligand_resname else None,
        "ligand_chain": ligand_chain.strip() if ligand_chain else None,
        "ligand_resid": int(ligand_resid) if ligand_resid is not None else None,
        "box_padding": float(box_padding),
        "schrodinger_root": schrodinger_root,
        "allow_fallback": bool(allow_fallback),
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
            "job_type": "Target_preparation",
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
        },
    )
    if not background:
        return_code = _worker(spec_path)
        if return_code:
            status = _read_json(status_path)
            raise RuntimeError(status.get("error", "Target preparation failed"))
        return json.dumps(_read_json(status_path), indent=2, ensure_ascii=False)

    log_handle = (job_dir / "job.log").open("a")
    try:
        process = subprocess.Popen(
            [sys.executable, "-m", "bead.tools.target_preparation", "--worker", str(spec_path)],
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


def target_preparation_job_status(job_path: str, log_lines: int = 30) -> str:
    """Return status and recent logs for a Target_preparation job."""
    path = Path(job_path).expanduser().resolve()
    status_path = path if path.name == "status.json" else path / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"Target preparation status file was not found: {status_path}")
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
