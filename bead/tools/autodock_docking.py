"""AutoDock Vina molecular docking jobs for BEAD.

The workflow mirrors :mod:`bead.tools.docking` while remaining independent of
Schrodinger.  It prepares ligands with Molscrub and Meeko, docks them with the
Vina Python API, and uses ProLIF for protein-ligand interaction annotations.
"""

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
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw

from bead.tools.ligand_input import materialize_ligand_csv


DEFAULT_RESULTS_ROOT = Path.home() / "bead_autodock_results"
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


def _run(command: list[str], cwd: Path) -> None:
    print("$ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


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
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            reason = None
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


def _pdbqt_to_pdb(pdbqt_path: Path, pdb_path: Path) -> None:
    lines: list[str] = []
    for line in pdbqt_path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atom_type = line.split()[-1]
            element = "".join(char for char in atom_type if char.isalpha())[:2].title()
            lines.append(f"{line[:76]:<76}{element:>2}")
    if not lines:
        raise ValueError(f"The receptor PDBQT contains no atoms: {pdbqt_path}")
    pdb_path.write_text("\n".join(lines) + "\nEND\n")


def _prepare_receptor(spec: dict[str, Any], work_dir: Path, output_dir: Path) -> tuple[Path, Path]:
    source = Path(spec["receptor_path"])
    receptor_pdbqt = work_dir / "receptor.pdbqt"
    receptor_pdb = output_dir / "receptor.pdb"
    if source.suffix.lower() == ".pdbqt":
        shutil.copy2(source, receptor_pdbqt)
        _pdbqt_to_pdb(receptor_pdbqt, receptor_pdb)
    elif source.suffix.lower() in {".pdb", ".ent"}:
        prepare = _environment_executable("mk_prepare_receptor.py")
        basename = work_dir / "receptor"
        _run(
            [
                str(prepare),
                "--read_pdb",
                str(source),
                "-o",
                str(basename),
                "-p",
                "--write_pdb",
                str(receptor_pdb),
                "-x",
                "--box_center",
                *(str(value) for value in spec["box_center"]),
                "--box_size",
                *(str(value) for value in spec["box_size"]),
                "-v",
            ],
            work_dir,
        )
        generated = basename.with_suffix(".pdbqt")
        if generated != receptor_pdbqt and generated.is_file():
            shutil.move(generated, receptor_pdbqt)
    else:
        raise ValueError("receptor_path must be a prepared PDBQT or a PDB/ENT file")
    if not receptor_pdbqt.is_file() or receptor_pdbqt.stat().st_size == 0:
        raise RuntimeError("Meeko did not produce a receptor PDBQT file")
    if not receptor_pdb.is_file():
        _pdbqt_to_pdb(receptor_pdbqt, receptor_pdb)
    config = work_dir / "vina_box.txt"
    config.write_text(
        "\n".join(
            [
                f"center_x = {spec['box_center'][0]}",
                f"center_y = {spec['box_center'][1]}",
                f"center_z = {spec['box_center'][2]}",
                f"size_x = {spec['box_size'][0]}",
                f"size_y = {spec['box_size'][1]}",
                f"size_z = {spec['box_size'][2]}",
                f"exhaustiveness = {spec['exhaustiveness']}",
                f"num_modes = {spec['num_modes']}",
                f"energy_range = {spec['energy_range']}",
                f"seed = {spec['seed']}",
            ]
        )
        + "\n"
    )
    shutil.copy2(config, output_dir / config.name)
    return receptor_pdbqt, receptor_pdb


def _limit_sdf_states(source: Path, destination: Path, maximum: int, input_id: str) -> int:
    molecules = [mol for mol in Chem.SDMolSupplier(str(source), removeHs=False) if mol]
    molecules = molecules[:maximum]
    if not molecules:
        return 0
    writer = Chem.SDWriter(str(destination))
    try:
        for number, molecule in enumerate(molecules, start=1):
            molecule.SetProp("_Name", f"{input_id}_state{number:02d}")
            writer.write(molecule)
    finally:
        writer.close()
    return len(molecules)


def _prepare_ligand_states(
    row: Any,
    spec: dict[str, Any],
    ligand_work_dir: Path,
) -> list[Path]:
    ligand_work_dir.mkdir(parents=True, exist_ok=True)
    scrub = _environment_executable("scrub.py")
    meeko = _environment_executable("mk_prepare_ligand.py")
    scrubbed = ligand_work_dir / "scrubbed.sdf"
    limited = ligand_work_dir / "prepared.sdf"
    _run(
        [
            str(scrub),
            row.smiles,
            "-o",
            str(scrubbed),
            "--ph_low",
            str(spec["ph"] - spec["ph_tolerance"]),
            "--ph_high",
            str(spec["ph"] + spec["ph_tolerance"]),
            "--cpu",
            "1",
            "--etkdg_rng_seed",
            str(spec["seed"] + int(row.source_index)),
            "--num_internal_confs",
            "1",
            "--min_output_confs",
            "1",
        ],
        ligand_work_dir,
    )
    count = _limit_sdf_states(
        scrubbed,
        limited,
        int(spec["max_states"]),
        row.input_id,
    )
    if count == 0:
        raise RuntimeError("Molscrub did not produce a prepared ligand")
    pdbqt_dir = ligand_work_dir / "pdbqt"
    pdbqt_dir.mkdir(exist_ok=True)
    _run(
        [
            str(meeko),
            "-i",
            str(limited),
            "--multimol_outdir",
            str(pdbqt_dir),
            "--multimol_prefix",
            row.input_id,
            "--add_index_map",
        ],
        ligand_work_dir,
    )
    states = sorted(pdbqt_dir.glob("*.pdbqt"))
    if not states:
        raise RuntimeError("Meeko did not produce a ligand PDBQT file")
    return states


def _pdbqt_poses_to_sdf(pdbqt_path: Path, sdf_path: Path) -> list[Chem.Mol]:
    from meeko import PDBQTMolecule, RDKitMolCreate

    docked = PDBQTMolecule.from_file(str(pdbqt_path))
    sdf_text, failures = RDKitMolCreate.write_sd_string(docked)
    if failures:
        print(f"Meeko export warnings for {pdbqt_path.name}: {failures}", flush=True)
    sdf_path.write_text(sdf_text)
    return [mol for mol in Chem.SDMolSupplier(str(sdf_path), removeHs=False) if mol]


def _element_from_pdb_line(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    name = line[12:16].strip()
    return "".join(char for char in name if char.isalpha())[:1].upper()


def _nearby_residues(
    receptor_pdb: Path,
    ligand: Chem.Mol,
    pocket_cutoff: float = 6.0,
    contact_cutoff: float = 4.0,
) -> tuple[list[str], list[dict[str, Any]]]:
    conformer = ligand.GetConformer()
    ligand_xyz = np.asarray(
        [
            list(conformer.GetAtomPosition(index))
            for index, atom in enumerate(ligand.GetAtoms())
            if atom.GetAtomicNum() != 1
        ],
        dtype=float,
    )
    residues: dict[tuple[str, str, str, str], list[tuple[str, np.ndarray, str]]] = {}
    for line in receptor_pdb.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            xyz = np.asarray(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
        except ValueError:
            continue
        key = (line[21:22].strip(), line[17:20].strip(), line[22:26].strip(), line[26:27].strip())
        residues.setdefault(key, []).append((line, xyz, _element_from_pdb_line(line)))

    pocket_lines: list[str] = []
    contacts: list[dict[str, Any]] = []
    for (chain, name, number, insertion), atoms in residues.items():
        minimum = float("inf")
        for _, xyz, element in atoms:
            if element == "H":
                continue
            minimum = min(minimum, float(np.linalg.norm(ligand_xyz - xyz, axis=1).min()))
        if minimum <= pocket_cutoff:
            pocket_lines.extend(line for line, _, _ in atoms)
        if minimum <= contact_cutoff:
            residue = f"{chain or '_'}:{name}{number}{insertion}"
            contacts.append(
                {
                    "residue": residue,
                    "minimum_distance_angstrom": round(minimum, 3),
                }
            )
    contacts.sort(key=lambda item: item["minimum_distance_angstrom"])
    return pocket_lines, contacts


def _atom_label(molecule: Any, indices: Any) -> str | None:
    if indices is None:
        return None
    if isinstance(indices, int):
        indices = (indices,)
    labels: list[str] = []
    for index in indices:
        atom = molecule.GetAtomWithIdx(int(index))
        info = atom.GetPDBResidueInfo()
        name = info.GetName().strip() if info else f"{atom.GetSymbol()}{int(index) + 1}"
        labels.append(name)
    return ",".join(labels)


def _interaction_analysis(
    receptor_pdb: Path,
    ligand: Chem.Mol,
    analysis_dir: Path,
) -> dict[str, Any]:
    details: dict[str, Any] = {
        "hydrogen_bonds": [],
        "salt_bridges": [],
        "hydrophobic_contacts": [],
        "pi_pi_interactions": [],
        "pi_cation_interactions": [],
        "close_contact_residues": [],
        "analysis_warnings": [],
    }
    pocket_lines, contacts = _nearby_residues(receptor_pdb, ligand)
    details["close_contact_residues"] = contacts
    if not pocket_lines:
        details["analysis_warnings"].append("No receptor residues were found within 6 A")
        return details
    pocket_path = analysis_dir / "interaction_pocket.pdb"
    pocket_path.write_text("\n".join(pocket_lines) + "\nEND\n")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import MDAnalysis as mda
            import prolif as plf

            protein = plf.Molecule.from_mda(mda.Universe(str(pocket_path)))
            prolif_ligand = plf.Molecule.from_rdkit(ligand)
            fingerprint = plf.Fingerprint()
            fingerprint.run_from_iterable(
                [prolif_ligand],
                protein,
                progress=False,
                n_jobs=1,
            )
        mapping = {
            "HBDonor": "hydrogen_bonds",
            "HBAcceptor": "hydrogen_bonds",
            "Anionic": "salt_bridges",
            "Cationic": "salt_bridges",
            "Hydrophobic": "hydrophobic_contacts",
            "PiStacking": "pi_pi_interactions",
            "PiCation": "pi_cation_interactions",
            "CationPi": "pi_cation_interactions",
        }
        for (ligand_residue, protein_residue), interactions in fingerprint.ifp.get(0, {}).items():
            residue = f"{protein_residue.chain or '_'}:{protein_residue.name}{protein_residue.number}"
            for prolif_type, occurrences in interactions.items():
                destination = mapping.get(prolif_type)
                if destination is None:
                    continue
                for occurrence in occurrences:
                    indices = occurrence.get("indices", {})
                    item: dict[str, Any] = {
                        "protein_residue": residue,
                        "protein_atom": _atom_label(protein, indices.get("protein")),
                        "ligand_atom": _atom_label(prolif_ligand, indices.get("ligand")),
                        "prolif_type": prolif_type,
                    }
                    if occurrence.get("distance") is not None:
                        item["distance_angstrom"] = round(float(occurrence["distance"]), 3)
                    if occurrence.get("DHA_angle") is not None:
                        item["angle_degree"] = round(float(occurrence["DHA_angle"]), 3)
                    details[destination].append(item)
    except Exception as exc:
        details["analysis_warnings"].append(f"ProLIF analysis failed: {type(exc).__name__}: {exc}")
    return details


def _complex_pdb(receptor_pdb: Path, ligand: Chem.Mol, output_path: Path) -> None:
    receptor_lines = receptor_pdb.read_text(errors="replace").splitlines()
    retained = [line for line in receptor_lines if not line.startswith(("END", "CONECT"))]
    maximum_serial = 0
    for line in retained:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                maximum_serial = max(maximum_serial, int(line[6:11]))
            except ValueError:
                pass
    ligand_lines = Chem.MolToPDBBlock(ligand).splitlines()
    serial_map: dict[int, int] = {}
    added: list[str] = []
    for line in ligand_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        old_serial = int(line[6:11])
        new_serial = maximum_serial + len(serial_map) + 1
        serial_map[old_serial] = new_serial
        rewritten = (
            f"HETATM{new_serial:5d}{line[11:17]}LIG Z{999:4d}{line[26:]}"
        )
        added.append(rewritten)
    for line in ligand_lines:
        if not line.startswith("CONECT"):
            continue
        fields = line.split()[1:]
        mapped = [serial_map.get(int(value)) for value in fields]
        mapped = [value for value in mapped if value is not None]
        if len(mapped) > 1:
            added.append("CONECT" + "".join(f"{value:5d}" for value in mapped))
    output_path.write_text("\n".join(retained + ["TER"] + added + ["END"]) + "\n")


def _interaction_diagram(
    ligand: Chem.Mol,
    details: dict[str, Any],
    score: float,
    input_id: str,
    output_path: Path,
) -> None:
    """Render a LigPlot-style 2D ligand-residue interaction diagram."""
    drawing_molecule = Chem.RemoveHs(Chem.Mol(ligand))
    Chem.rdDepictor.Compute2DCoords(drawing_molecule)
    molecule_image = Draw.MolToImage(drawing_molecule, size=(620, 470))
    canvas = Image.new("RGB", (1200, 900), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        title_font = ImageFont.truetype("DejaVuSans-Bold.ttf", 28)
        node_font = ImageFont.truetype("DejaVuSans.ttf", 17)
        legend_font = ImageFont.truetype("DejaVuSans.ttf", 16)
    except OSError:
        title_font = node_font = legend_font = ImageFont.load_default()

    draw.text(
        (24, 18),
        f"{input_id} | AutoDock Vina affinity: {score:.3f} kcal/mol",
        fill="#17202a",
        font=title_font,
    )
    center = (600, 460)
    residue_types: dict[str, set[str]] = {}
    categories = [
        ("H-bond", "hydrogen_bonds", "#2980b9"),
        ("Salt bridge", "salt_bridges", "#c0392b"),
        ("Hydrophobic", "hydrophobic_contacts", "#27ae60"),
        ("Pi-stacking", "pi_pi_interactions", "#8e44ad"),
        ("Pi-cation", "pi_cation_interactions", "#d35400"),
    ]
    colors = {label: color for label, _, color in categories}
    for label, key, _ in categories:
        for item in details.get(key, []):
            residue_types.setdefault(item["protein_residue"], set()).add(label)
    for contact in details.get("close_contact_residues", []):
        residue_types.setdefault(contact["residue"], set()).add("Contact")
    colors["Contact"] = "#7f8c8d"

    residues = sorted(
        residue_types,
        key=lambda residue: (
            residue_types[residue] == {"Contact"},
            residue,
        ),
    )[:18]
    positions: dict[str, tuple[int, int]] = {}
    for index, residue in enumerate(residues):
        angle = -math.pi / 2 + 2 * math.pi * index / max(len(residues), 1)
        positions[residue] = (
            int(center[0] + 500 * math.cos(angle)),
            int(center[1] + 340 * math.sin(angle)),
        )

    # Draw connectors first so that the ligand depiction remains unobscured.
    for residue, position in positions.items():
        interaction_types = residue_types[residue]
        primary = next(
            (label for label, _, _ in categories if label in interaction_types),
            "Contact",
        )
        width = 5 if primary != "Contact" else 2
        draw.line([center, position], fill=colors[primary], width=width)

    canvas.paste(molecule_image, (290, 225))
    for residue, (x, y) in positions.items():
        interaction_types = residue_types[residue]
        primary = next(
            (label for label, _, _ in categories if label in interaction_types),
            "Contact",
        )
        label = residue + "\n" + "/".join(sorted(interaction_types))
        box = (x - 84, y - 28, x + 84, y + 28)
        draw.rounded_rectangle(box, radius=10, fill="white", outline=colors[primary], width=4)
        bbox = draw.multiline_textbbox((0, 0), label, font=node_font, align="center", spacing=2)
        width = bbox[2] - bbox[0]
        height = bbox[3] - bbox[1]
        draw.multiline_text(
            (x - width / 2, y - height / 2),
            label,
            fill="#17202a",
            font=node_font,
            align="center",
            spacing=2,
        )

    legend_x, legend_y = 22, 820
    for label, _, color in categories + [("Contact", "", colors["Contact"])]:
        draw.line([(legend_x, legend_y + 8), (legend_x + 34, legend_y + 8)], fill=color, width=5)
        draw.text((legend_x + 42, legend_y), label, fill="#17202a", font=legend_font)
        legend_x += 178
    canvas.save(output_path)


def _binding_3d_diagram(
    receptor_pdb: Path,
    ligand: Chem.Mol,
    details: dict[str, Any],
    score: float,
    input_id: str,
    output_path: Path,
) -> None:
    """Render a static 3D pocket view with ligand and interaction guides."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    conformer = ligand.GetConformer()
    ligand_xyz = np.asarray(
        [list(conformer.GetAtomPosition(index)) for index in range(ligand.GetNumAtoms())],
        dtype=float,
    )
    ligand_heavy = np.asarray(
        [
            ligand_xyz[index]
            for index, atom in enumerate(ligand.GetAtoms())
            if atom.GetAtomicNum() != 1
        ]
    )
    receptor_residues: dict[str, list[tuple[str, str, np.ndarray]]] = {}
    for line in receptor_pdb.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = _element_from_pdb_line(line)
        if element == "H":
            continue
        try:
            xyz = np.asarray(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
        except ValueError:
            continue
        chain = line[21:22].strip() or "_"
        residue = f"{chain}:{line[17:20].strip()}{line[22:26].strip()}{line[26:27].strip()}"
        if float(np.linalg.norm(ligand_heavy - xyz, axis=1).min()) <= 5.0:
            receptor_residues.setdefault(residue, []).append((line[12:16].strip(), element, xyz))

    element_colors = {
        "C": "#b0b0b0",
        "N": "#2878d0",
        "O": "#d62728",
        "S": "#e6b800",
        "P": "#e67e22",
        "F": "#2ca02c",
        "CL": "#2ca02c",
        "BR": "#8c564b",
        "I": "#9467bd",
    }
    interaction_colors = {
        "hydrogen_bonds": "#00a6d6",
        "salt_bridges": "#e74c3c",
        "hydrophobic_contacts": "#27ae60",
        "pi_pi_interactions": "#8e44ad",
        "pi_cation_interactions": "#f39c12",
    }
    interacting_residues = {
        item["protein_residue"]
        for key in interaction_colors
        for item in details.get(key, [])
    }

    figure = plt.figure(figsize=(10, 8), dpi=180)
    axis = figure.add_subplot(111, projection="3d")
    axis.set_facecolor("white")

    for residue, atoms in receptor_residues.items():
        coordinates = np.asarray([atom[2] for atom in atoms])
        for first in range(len(atoms)):
            for second in range(first + 1, len(atoms)):
                distance = float(np.linalg.norm(coordinates[first] - coordinates[second]))
                if 0.5 < distance <= 1.95:
                    axis.plot(
                        [coordinates[first, 0], coordinates[second, 0]],
                        [coordinates[first, 1], coordinates[second, 1]],
                        [coordinates[first, 2], coordinates[second, 2]],
                        color="#9aa0a6",
                        linewidth=1.2,
                        alpha=0.65,
                    )
        for _, element, xyz in atoms:
            axis.scatter(
                xyz[0],
                xyz[1],
                xyz[2],
                color=element_colors.get(element, "#9e9e9e"),
                s=18,
                alpha=0.7,
                depthshade=True,
            )
        nearest_index = int(
            np.argmin([float(np.linalg.norm(ligand_heavy - xyz, axis=1).min()) for _, _, xyz in atoms])
        )
        label_xyz = atoms[nearest_index][2]
        axis.text(
            label_xyz[0],
            label_xyz[1],
            label_xyz[2],
            residue,
            fontsize=7,
            color="#641e16" if residue in interacting_residues else "#34495e",
        )

    for bond in ligand.GetBonds():
        first = ligand_xyz[bond.GetBeginAtomIdx()]
        second = ligand_xyz[bond.GetEndAtomIdx()]
        axis.plot(
            [first[0], second[0]],
            [first[1], second[1]],
            [first[2], second[2]],
            color="#202020",
            linewidth=3.0,
            zorder=10,
        )
    for index, atom in enumerate(ligand.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue
        xyz = ligand_xyz[index]
        axis.scatter(
            xyz[0],
            xyz[1],
            xyz[2],
            color=element_colors.get(atom.GetSymbol().upper(), "#333333"),
            edgecolor="black",
            linewidth=0.5,
            s=85,
            depthshade=True,
            zorder=12,
        )

    for key, color in interaction_colors.items():
        for item in details.get(key, []):
            residue = item["protein_residue"]
            atoms = receptor_residues.get(residue)
            if not atoms:
                continue
            receptor_xyz = np.asarray([atom[2] for atom in atoms])
            distances = np.linalg.norm(receptor_xyz[:, None, :] - ligand_heavy[None, :, :], axis=2)
            receptor_index, ligand_index = np.unravel_index(np.argmin(distances), distances.shape)
            start = receptor_xyz[receptor_index]
            end = ligand_heavy[ligand_index]
            axis.plot(
                [start[0], end[0]],
                [start[1], end[1]],
                [start[2], end[2]],
                color=color,
                linestyle="--",
                linewidth=2.4,
                alpha=0.95,
                zorder=11,
            )

    all_xyz = np.vstack(
        [ligand_heavy]
        + [np.asarray([atom[2] for atom in atoms]) for atoms in receptor_residues.values()]
    )
    center = all_xyz.mean(axis=0)
    radius = max(float(np.ptp(all_xyz[:, dimension])) / 2 for dimension in range(3)) + 1.5
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1))
    axis.view_init(elev=22, azim=-58)
    axis.set_axis_off()
    axis.set_title(
        f"{input_id} — AutoDock Vina binding pose ({score:.3f} kcal/mol)",
        fontsize=13,
        pad=16,
    )
    legend = [
        Line2D([0], [0], marker="o", color="black", label="Docked ligand", markersize=7),
        Line2D([0], [0], color="#9aa0a6", label="Binding-site residues", linewidth=2),
    ]
    legend.extend(
        Line2D([0], [0], color=color, linestyle="--", label=label, linewidth=2)
        for label, color in [
            ("Hydrogen bond", interaction_colors["hydrogen_bonds"]),
            ("Salt bridge", interaction_colors["salt_bridges"]),
            ("Hydrophobic", interaction_colors["hydrophobic_contacts"]),
            ("Pi-stacking", interaction_colors["pi_pi_interactions"]),
            ("Pi-cation", interaction_colors["pi_cation_interactions"]),
        ]
    )
    axis.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.0, 0.98), fontsize=8)
    figure.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _interaction_table(records: list[dict[str, Any]], output_path: Path) -> None:
    rows: list[dict[str, Any]] = []
    groups = [
        ("hydrogen_bond", "hydrogen_bonds"),
        ("salt_bridge", "salt_bridges"),
        ("hydrophobic_contact", "hydrophobic_contacts"),
        ("pi_pi", "pi_pi_interactions"),
        ("pi_cation", "pi_cation_interactions"),
    ]
    for record in records:
        interactions = record["interactions"]
        for interaction_type, key in groups:
            for item in interactions.get(key, []):
                rows.append(
                    {
                        "input_id": record["input_id"],
                        "interaction_type": interaction_type,
                        **item,
                    }
                )
        for item in interactions.get("close_contact_residues", []):
            rows.append(
                {
                    "input_id": record["input_id"],
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
        "angle_degree",
        "prolif_type",
    ]
    pd.DataFrame(rows).reindex(columns=columns).to_csv(output_path, index=False)


def _execute(spec: dict[str, Any]) -> dict[str, Any]:
    from vina import Vina

    job_dir = Path(spec["job_dir"])
    work_dir = job_dir / "work"
    output_dir = job_dir / "output"
    work_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    for name in (
        "docked_ligands",
        "raw_poses",
        "complexes",
        "interaction_diagrams",
        "binding_3d",
    ):
        (output_dir / name).mkdir(exist_ok=True)

    selected, skipped = _prepare_inputs(spec)
    selected.to_csv(output_dir / "selected_inputs.csv", index=False)
    skipped.to_csv(output_dir / "skipped_inputs.csv", index=False)
    receptor_pdbqt, receptor_pdb = _prepare_receptor(spec, work_dir, output_dir)

    vina_engine = Vina(
        sf_name=spec["scoring"],
        cpu=int(spec["cpu"]),
        seed=int(spec["seed"]),
        verbosity=0,
    )
    vina_engine.set_receptor(str(receptor_pdbqt))
    vina_engine.compute_vina_maps(
        center=[float(value) for value in spec["box_center"]],
        box_size=[float(value) for value in spec["box_size"]],
    )

    records: list[dict[str, Any]] = []
    raw_rows: list[dict[str, Any]] = []
    failures: dict[str, str] = {}
    ligand_work_root = work_dir / "ligands"
    for row_number, row in enumerate(selected.itertuples(index=False), start=1):
        print(f"[{row_number}/{len(selected)}] Preparing and docking {row.input_id}", flush=True)
        try:
            states = _prepare_ligand_states(
                row,
                spec,
                ligand_work_root / row.input_id,
            )
            best: tuple[float, int, Path, Chem.Mol, list[float]] | None = None
            for state_number, state_path in enumerate(states, start=1):
                vina_engine.set_ligand_from_file(str(state_path))
                vina_engine.dock(
                    exhaustiveness=int(spec["exhaustiveness"]),
                    n_poses=int(spec["num_modes"]),
                )
                energies = vina_engine.energies(n_poses=int(spec["num_modes"]))
                if energies.size == 0:
                    continue
                pose_path = output_dir / "raw_poses" / f"{row.input_id}_state{state_number:02d}.pdbqt"
                vina_engine.write_poses(
                    str(pose_path),
                    n_poses=int(spec["num_modes"]),
                    energy_range=float(spec["energy_range"]),
                    overwrite=True,
                )
                state_sdf = ligand_work_root / row.input_id / f"state{state_number:02d}_poses.sdf"
                pose_molecules = _pdbqt_poses_to_sdf(pose_path, state_sdf)
                if not pose_molecules:
                    continue
                for mode, energy in enumerate(energies.tolist(), start=1):
                    raw_rows.append(
                        {
                            "input_id": row.input_id,
                            "state": state_number,
                            "mode": mode,
                            "vina_affinity_kcal_mol": float(energy[0]),
                            "energy_components": json.dumps(energy),
                            "pose_pdbqt": str(pose_path),
                        }
                    )
                score = float(energies[0][0])
                if best is None or score < best[0]:
                    best = (score, state_number, pose_path, pose_molecules[0], energies[0].tolist())
            if best is None:
                failures[row.input_id] = "Vina produced no docked pose"
                continue

            score, state_number, pose_path, best_molecule, components = best
            best_molecule.SetProp("_Name", row.input_id)
            best_molecule.SetProp("vina_affinity_kcal_mol", f"{score:.6f}")
            best_molecule.SetProp("source_smiles", row.smiles)
            best_molecule.SetIntProp("prepared_state", state_number)
            ligand_path = output_dir / "docked_ligands" / f"{row.input_id}.sdf"
            writer = Chem.SDWriter(str(ligand_path))
            writer.write(best_molecule)
            writer.close()

            analysis_dir = ligand_work_root / row.input_id / "analysis"
            analysis_dir.mkdir(exist_ok=True)
            interactions = _interaction_analysis(receptor_pdb, best_molecule, analysis_dir)
            complex_path = output_dir / "complexes" / f"Inter_{row_number}.pdb"
            _complex_pdb(receptor_pdb, best_molecule, complex_path)
            diagram_path = output_dir / "interaction_diagrams" / f"{row.input_id}.png"
            _interaction_diagram(best_molecule, interactions, score, row.input_id, diagram_path)
            binding_3d_path = output_dir / "binding_3d" / f"{row.input_id}.png"
            _binding_3d_diagram(
                receptor_pdb,
                best_molecule,
                interactions,
                score,
                row.input_id,
                binding_3d_path,
            )
            records.append(
                {
                    "input_id": row.input_id,
                    "vina_affinity_kcal_mol": score,
                    "scoring_function": spec["scoring"],
                    "prepared_state": state_number,
                    "energy_components": components,
                    "complex_pdb": str(complex_path),
                    "docked_ligand_sdf": str(ligand_path),
                    "raw_pose_pdbqt": str(pose_path),
                    "interaction_diagram": str(diagram_path),
                    "interaction_diagram_2d": str(diagram_path),
                    "binding_3d_image": str(binding_3d_path),
                    "hydrogen_bond_count": len(interactions["hydrogen_bonds"]),
                    "salt_bridge_count": len(interactions["salt_bridges"]),
                    "hydrophobic_contact_count": len(interactions["hydrophobic_contacts"]),
                    "pi_pi_count": len(interactions["pi_pi_interactions"]),
                    "pi_cation_count": len(interactions["pi_cation_interactions"]),
                    "contact_residue_count": len(interactions["close_contact_residues"]),
                    "interactions": interactions,
                }
            )
        except Exception as exc:
            failures[row.input_id] = f"{type(exc).__name__}: {exc}"
            traceback.print_exc()

    if not records:
        raise RuntimeError(f"No molecules were docked successfully: {failures}")

    report_json = output_dir / "interaction_details.json"
    _write_json(report_json, records)
    interaction_csv = output_dir / "interaction_details.csv"
    _interaction_table(records, interaction_csv)
    raw_csv = output_dir / "vina_raw_results.csv"
    pd.DataFrame(raw_rows).to_csv(raw_csv, index=False)

    flat_records = [
        {key: value for key, value in record.items() if key not in {"interactions", "energy_components"}}
        for record in records
    ]
    results = selected.merge(pd.DataFrame(flat_records), on="input_id", how="left")
    results["docking_status"] = results["vina_affinity_kcal_mol"].notna().map(
        {True: "docked", False: "failed"}
    )
    results["failure_reason"] = results["input_id"].map(failures)
    source = pd.read_csv(spec["csv_path"]).reset_index(names="source_index")
    source = source.rename(
        columns={
            column: f"source_{column}"
            for column in source.columns
            if column != "source_index"
        }
    )
    results = results.merge(source, on="source_index", how="left")
    results = results.sort_values(
        ["vina_affinity_kcal_mol", "source_index"],
        na_position="last",
    )
    results_path = output_dir / "docking_results.csv"
    results.to_csv(results_path, index=False)
    return {
        "input_rows": int(pd.read_csv(spec["csv_path"]).shape[0]),
        "selected_rows": int(len(selected)),
        "skipped_rows": int(len(skipped)),
        "docked_rows": int(len(records)),
        "failed_rows": int(len(failures)),
        "results_csv": str(results_path),
        "interaction_details_csv": str(interaction_csv),
        "interaction_details_json": str(report_json),
        "receptor_pdbqt": str(receptor_pdbqt),
        "output_dir": str(output_dir),
    }


def _artifacts(job_dir: Path) -> list[str]:
    keep = {".csv", ".json", ".log", ".txt", ".pdb", ".pdbqt", ".sdf", ".png"}
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
            "job_type": "Docking_autodock",
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
                "job_type": "Docking_autodock",
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
            "job_type": "Docking_autodock",
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


def autodock_docking(
    csv_path: str | None,
    receptor_path: str,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    smiles_column: str | None = None,
    smiles: list[str] | None = None,
    molecule_names: list[str] | None = None,
    scoring: str = "vina",
    exhaustiveness: int = 16,
    num_modes: int = 9,
    energy_range: float = 3.0,
    seed: int = 20260928,
    cpu: int = 1,
    ph: float = 7.4,
    ph_tolerance: float = 1.0,
    max_states: int = 8,
    max_molecules: int | None = None,
    output_root: str | None = None,
    background: bool = True,
) -> str:
    """Dock direct SMILES or a SMILES CSV against a receptor with AutoDock Vina."""
    receptor_file = _input_file(receptor_path, "receptor")
    scoring = str(scoring).strip().lower()
    if scoring not in {"vina", "vinardo"}:
        raise ValueError("scoring must be 'vina' or 'vinardo'")
    box_size = [float(size_x), float(size_y), float(size_z)]
    if any(value <= 0 for value in box_size):
        raise ValueError("All search-box dimensions must be greater than zero")
    if int(exhaustiveness) < 1 or int(num_modes) < 1 or int(cpu) < 1:
        raise ValueError("exhaustiveness, num_modes, and cpu must be at least 1")
    if float(energy_range) < 0 or float(ph_tolerance) < 0:
        raise ValueError("energy_range and ph_tolerance cannot be negative")
    if int(max_states) < 1:
        raise ValueError("max_states must be at least 1")
    if max_molecules is not None and int(max_molecules) < 1:
        raise ValueError("max_molecules must be at least 1 when supplied")

    root = Path(output_root or DEFAULT_RESULTS_ROOT).expanduser().resolve()
    date_dir = root / datetime.now().astimezone().date().isoformat()
    job_id = f"autodock-docking-{datetime.now().strftime('%H%M%S')}-{uuid.uuid4().hex[:8]}"
    job_dir = date_dir / job_id
    inputs_dir = job_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=False)
    local_csv, input_mode = materialize_ligand_csv(
        inputs_dir,
        csv_path=csv_path,
        smiles=smiles,
        molecule_names=molecule_names,
    )
    local_receptor = inputs_dir / receptor_file.name
    shutil.copy2(receptor_file, local_receptor)
    spec = {
        "job_id": job_id,
        "job_dir": str(job_dir),
        "csv_path": str(local_csv),
        "input_mode": input_mode,
        "receptor_path": str(local_receptor),
        "smiles_column": "smiles" if input_mode == "direct_smiles" else smiles_column,
        "box_center": [float(center_x), float(center_y), float(center_z)],
        "box_size": box_size,
        "scoring": scoring,
        "exhaustiveness": int(exhaustiveness),
        "num_modes": int(num_modes),
        "energy_range": float(energy_range),
        "seed": int(seed),
        "cpu": int(cpu),
        "ph": float(ph),
        "ph_tolerance": float(ph_tolerance),
        "max_states": int(max_states),
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
            "job_type": "Docking_autodock",
            "job_dir": str(job_dir),
            "log": str(job_dir / "job.log"),
        },
    )
    if not background:
        return_code = _worker(spec_path)
        if return_code:
            status = _read_json(status_path)
            raise RuntimeError(status.get("error", "AutoDock docking failed"))
        return json.dumps(_read_json(status_path), indent=2, ensure_ascii=False)

    log_handle = (job_dir / "job.log").open("a")
    try:
        process = subprocess.Popen(
            [sys.executable, "-m", "bead.tools.autodock_docking", "--worker", str(spec_path)],
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


def autodock_job_status(job_path: str, log_lines: int = 30) -> str:
    """Return status and recent log lines for one AutoDock docking job."""
    path = Path(job_path).expanduser().resolve()
    status_path = path if path.name == "status.json" else path / "status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"AutoDock status file was not found: {status_path}")
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
