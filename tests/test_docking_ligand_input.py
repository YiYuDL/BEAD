import csv
from pathlib import Path

import pytest

from bead.tools.ligand_input import materialize_ligand_csv


def test_materializes_named_direct_smiles(tmp_path: Path) -> None:
    output, mode = materialize_ligand_csv(
        tmp_path,
        csv_path=None,
        smiles=["CCO", "CCN"],
        molecule_names=["Ref", "New"],
    )

    assert mode == "direct_smiles"
    with output.open(newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["name", "smiles"],
            ["Ref", "CCO"],
            ["New", "CCN"],
        ]


def test_copies_an_existing_csv(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    input_dir = tmp_path / "inputs"
    source_dir.mkdir()
    input_dir.mkdir()
    source = source_dir / "molecules.csv"
    source.write_text("name,smiles\nRef,CCO\n")

    output, mode = materialize_ligand_csv(
        input_dir,
        csv_path=str(source),
        smiles=None,
        molecule_names=None,
    )

    assert mode == "csv"
    assert output.read_text() == source.read_text()


@pytest.mark.parametrize(
    ("csv_path", "smiles"),
    [(None, None), ("molecules.csv", ["CCO"])],
)
def test_requires_exactly_one_input_mode(
    tmp_path: Path, csv_path: str | None, smiles: list[str] | None
) -> None:
    with pytest.raises(ValueError, match="exactly one"):
        materialize_ligand_csv(
            tmp_path,
            csv_path=csv_path,
            smiles=smiles,
            molecule_names=None,
        )


def test_requires_one_name_per_smiles(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="same number"):
        materialize_ligand_csv(
            tmp_path,
            csv_path=None,
            smiles=["CCO", "CCN"],
            molecule_names=["Ref"],
        )
