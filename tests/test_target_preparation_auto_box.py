from pathlib import Path

import pytest

from bead.tools.target_preparation import _automatic_ligand_box


def _hetatm(
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"HETATM{serial:5d} {atom_name:<4} {resname:>3} {chain:1}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2}"
    )


def _write_fixture(path: Path) -> None:
    lines = [
        _hetatm(1, "ZN", "ZN", "A", 1, 50.0, 50.0, 50.0, "ZN"),
        _hetatm(2, "C1", "GOL", "A", 2, 20.0, 0.0, 0.0, "C"),
        _hetatm(3, "C2", "GOL", "A", 2, 21.0, 0.0, 0.0, "C"),
        _hetatm(4, "C3", "GOL", "A", 2, 22.0, 0.0, 0.0, "C"),
        _hetatm(5, "O1", "GOL", "A", 2, 20.0, 1.0, 0.0, "O"),
        _hetatm(6, "O2", "GOL", "A", 2, 21.0, 1.0, 0.0, "O"),
        _hetatm(7, "O3", "GOL", "A", 2, 22.0, 1.0, 0.0, "O"),
        _hetatm(8, "C1", "LIG", "B", 401, 0.0, 0.0, 0.0, "C"),
        _hetatm(9, "C2", "LIG", "B", 401, 2.0, 0.0, 0.0, "C"),
        _hetatm(10, "C3", "LIG", "B", 401, 4.0, 0.0, 0.0, "C"),
        _hetatm(11, "C4", "LIG", "B", 401, 6.0, 0.0, 0.0, "C"),
        _hetatm(12, "N1", "LIG", "B", 401, 2.0, 2.0, 0.0, "N"),
        _hetatm(13, "O1", "LIG", "B", 401, 4.0, 2.0, 0.0, "O"),
        "END",
    ]
    path.write_text("\n".join(lines) + "\n")


def test_auto_box_selects_largest_drug_like_heterogen(tmp_path: Path) -> None:
    pdb_path = tmp_path / "target.pdb"
    _write_fixture(pdb_path)

    result = _automatic_ligand_box(
        pdb_path,
        ligand_resname=None,
        ligand_chain=None,
        ligand_resid=None,
        box_padding=2.0,
        minimum_sizes=[10.0, 10.0, 10.0],
    )

    assert result is not None
    assert result["ligand"] == {
        "resname": "LIG",
        "chain": "B",
        "resid": 401,
        "insertion_code": "",
        "heavy_atom_count": 6,
        "carbon_atom_count": 4,
    }
    assert result["box_center"] == pytest.approx([3.0, 2.0 / 3.0, 0.0])
    assert result["ligand_extents"] == pytest.approx([6.0, 2.0, 0.0])
    assert result["box_size"] == pytest.approx([10.0, 10.0, 10.0])
    assert result["eligible_candidate_count"] == 1


def test_explicit_selector_can_choose_an_excluded_additive(tmp_path: Path) -> None:
    pdb_path = tmp_path / "target.pdb"
    _write_fixture(pdb_path)

    result = _automatic_ligand_box(
        pdb_path,
        ligand_resname="GOL",
        ligand_chain="A",
        ligand_resid=2,
        box_padding=1.0,
        minimum_sizes=[4.0, 4.0, 4.0],
    )

    assert result is not None
    assert result["selection_method"] == "explicit_selector"
    assert result["ligand"]["resname"] == "GOL"
    assert result["box_center"] == pytest.approx([21.0, 0.5, 0.0])


def test_explicit_selector_reports_a_missing_ligand(tmp_path: Path) -> None:
    pdb_path = tmp_path / "target.pdb"
    _write_fixture(pdb_path)

    with pytest.raises(ValueError, match="No HETATM ligand matched"):
        _automatic_ligand_box(
            pdb_path,
            ligand_resname="BAD",
            ligand_chain="A",
            ligand_resid=999,
            box_padding=4.0,
            minimum_sizes=[20.0, 20.0, 20.0],
        )
