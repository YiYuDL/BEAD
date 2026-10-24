"""Create BEAD docking artifacts with the Schrodinger Python runtime.

This module intentionally supports Python 3.8 because it is launched through
``$SCHRODINGER/run`` rather than through the BEAD Conda environment.
"""

from __future__ import print_function

import argparse
import json
import math
from pathlib import Path

from schrodinger import structure
from schrodinger.structutils.interactions import hbond, pi, salt_bridge


ID_PROPERTY_SUFFIX = "bead\\_input\\_id"


def _atom_label(atom):
    residue = (getattr(atom, "pdbres", "") or "UNK").strip()
    chain = (getattr(atom, "chain", "") or "").strip()
    number = getattr(atom, "resnum", 0)
    insertion = (getattr(atom, "inscode", "") or "").strip()
    atom_name = (getattr(atom, "pdbname", "") or "").strip()
    atom_name = atom_name or (getattr(atom, "atom_name", "") or "").strip()
    atom_name = atom_name or "%s%s" % (getattr(atom, "element", "X"), atom.index)
    residue_id = "%s:%s%s%s" % (chain or "_", residue, number, insertion)
    return "%s:%s" % (residue_id, atom_name), residue_id


def _distance(atom1, atom2):
    return math.sqrt(sum((float(a) - float(b)) ** 2 for a, b in zip(atom1.xyz, atom2.xyz)))


def _protein_ligand_atoms(atom1, atom2, receptor):
    """Return a pair in protein, ligand order (H-bond output is donor-first)."""
    if atom1._ct == receptor:
        return atom1, atom2
    return atom2, atom1


def _input_id(ligand, fallback):
    for key, value in ligand.property.items():
        if key.endswith(ID_PROPERTY_SUFFIX) or key.replace("\\_", "_").endswith("bead_input_id"):
            return str(value)
    return fallback


def _score(ligand):
    for key in ("r_i_docking_score", "r_i_glide_gscore"):
        if key in ligand.property:
            return float(ligand.property[key])
    return None


def _close_contact_residues(receptor, ligand, cutoff=4.0):
    contacts = {}
    receptor_atoms = [atom for atom in receptor.atom if atom.element != "H"]
    ligand_atoms = [atom for atom in ligand.atom if atom.element != "H"]
    cutoff_squared = cutoff * cutoff
    for rec_atom in receptor_atoms:
        rec_xyz = rec_atom.xyz
        best = None
        for lig_atom in ligand_atoms:
            lig_xyz = lig_atom.xyz
            distance_squared = sum(
                (float(rec_xyz[index]) - float(lig_xyz[index])) ** 2
                for index in range(3)
            )
            if distance_squared <= cutoff_squared and (best is None or distance_squared < best):
                best = distance_squared
        if best is not None:
            _, residue_id = _atom_label(rec_atom)
            distance = math.sqrt(best)
            contacts[residue_id] = min(distance, contacts.get(residue_id, distance))
    return [
        {"residue": residue, "minimum_distance_angstrom": round(distance, 3)}
        for residue, distance in sorted(contacts.items(), key=lambda item: item[1])
    ]


def _interactions(receptor, ligand):
    details = {
        "hydrogen_bonds": [],
        "salt_bridges": [],
        "pi_pi_count": 0,
        "pi_cation_count": 0,
        "close_contact_residues": [],
        "analysis_warnings": [],
    }
    try:
        for atom1, atom2 in hbond.get_hydrogen_bonds(receptor, st2=ligand):
            rec_atom, lig_atom = _protein_ligand_atoms(atom1, atom2, receptor)
            rec_label, rec_residue = _atom_label(rec_atom)
            lig_label, _ = _atom_label(lig_atom)
            details["hydrogen_bonds"].append(
                {
                    "protein_atom": rec_label,
                    "protein_residue": rec_residue,
                    "ligand_atom": lig_label,
                    "distance_angstrom": round(_distance(rec_atom, lig_atom), 3),
                }
            )
    except Exception as exc:
        details["analysis_warnings"].append("hydrogen bonds: %s" % exc)

    try:
        bridges = salt_bridge.get_salt_bridges(
            receptor,
            struc2=ligand,
            order_by=salt_bridge.OrderBy.InputOrder,
        )
        for atom1, atom2 in bridges:
            rec_atom, lig_atom = _protein_ligand_atoms(atom1, atom2, receptor)
            rec_label, rec_residue = _atom_label(rec_atom)
            lig_label, _ = _atom_label(lig_atom)
            details["salt_bridges"].append(
                {
                    "protein_atom": rec_label,
                    "protein_residue": rec_residue,
                    "ligand_atom": lig_label,
                    "distance_angstrom": round(_distance(rec_atom, lig_atom), 3),
                }
            )
    except Exception as exc:
        details["analysis_warnings"].append("salt bridges: %s" % exc)

    try:
        details["pi_pi_count"] = len(pi.find_pi_pi_interactions(receptor, struct2=ligand))
    except Exception as exc:
        details["analysis_warnings"].append("pi-pi interactions: %s" % exc)

    try:
        interactions = pi.find_pi_cation_interactions(receptor, struct2=ligand)
        details["pi_cation_count"] = len(interactions)
    except Exception as exc:
        details["analysis_warnings"].append("pi-cation interactions: %s" % exc)

    try:
        details["close_contact_residues"] = _close_contact_residues(receptor, ligand)
    except Exception as exc:
        details["analysis_warnings"].append("close contacts: %s" % exc)
    return details


def build_report(poseviewer, output_dir, report_json):
    output_dir.mkdir(parents=True, exist_ok=True)
    complexes_dir = output_dir / "complexes"
    ligands_dir = output_dir / "docked_ligands"
    complexes_dir.mkdir(exist_ok=True)
    ligands_dir.mkdir(exist_ok=True)

    with structure.StructureReader(str(poseviewer)) as reader:
        structures = list(reader)
    if len(structures) < 2:
        raise RuntimeError("The Glide pose-viewer contains no docked ligand poses")

    receptor = structures[0]
    receptor.write(str(output_dir / "receptor.pdb"))
    best = {}
    for number, ligand in enumerate(structures[1:], start=1):
        input_id = _input_id(ligand, "BEAD_%06d" % number)
        score = _score(ligand)
        current = best.get(input_id)
        if current is None or (
            score is not None
            and (current[0] is None or score < current[0])
        ):
            best[input_id] = (score, ligand)

    def order_key(item):
        input_id = item[0]
        try:
            return int(input_id.rsplit("_", 1)[1])
        except Exception:
            return input_id

    selected = sorted(best.items(), key=order_key)
    poseviewer_best = output_dir / "poseviewer_best.mae"
    rows = []
    with structure.StructureWriter(str(poseviewer_best)) as writer:
        writer.append(receptor)
        for output_index, (input_id, (score, ligand)) in enumerate(selected, start=1):
            ligand.title = input_id
            writer.append(ligand)
            ligand_path = ligands_dir / (input_id + ".sdf")
            ligand.write(str(ligand_path))
            complex_path = complexes_dir / ("Inter_%d.pdb" % output_index)
            complex_structure = receptor.merge(ligand, copy_props=True)
            complex_structure.title = input_id
            complex_structure.write(str(complex_path))
            details = _interactions(receptor, ligand)
            rows.append(
                {
                    "input_id": input_id,
                    "glide_docking_score": score,
                    "glide_gscore": ligand.property.get("r_i_glide_gscore"),
                    "glide_emodel": ligand.property.get("r_i_glide_emodel"),
                    "glide_energy": ligand.property.get("r_i_glide_energy"),
                    "complex_pdb": str(complex_path),
                    "docked_ligand_sdf": str(ligand_path),
                    "hydrogen_bond_count": len(details["hydrogen_bonds"]),
                    "salt_bridge_count": len(details["salt_bridges"]),
                    "pi_pi_count": details["pi_pi_count"],
                    "pi_cation_count": details["pi_cation_count"],
                    "contact_residue_count": len(details["close_contact_residues"]),
                    "interactions": details,
                }
            )

    report_json.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--poseviewer", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report-json", type=Path, required=True)
    args = parser.parse_args()
    build_report(args.poseviewer.resolve(), args.output_dir.resolve(), args.report_json.resolve())


if __name__ == "__main__":
    main()
