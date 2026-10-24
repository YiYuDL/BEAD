"""Lightweight AutoDock Vina score node for REINVENT/LibINVENT feedback."""

from pathlib import Path
from types import SimpleNamespace
from typing import Annotated, Literal

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem

from maize.core.interface import FileParameter, Input, Output, Parameter, Suffix
from maize.core.node import Node

from bead.tools.autodock_docking import _prepare_ligand_states


class AutoDockVinaScore(Node):
    """Return one Vina affinity per input SMILES in the original order.

    The node uses the same Molscrub/Meeko ligand preparation and Vina Python
    API as :func:`bead.tools.autodock_docking.autodock_docking`, while omitting
    pose rendering and interaction analysis that would be too costly inside an
    RL epoch.
    """

    tags = {"chemistry", "docking", "scorer", "reinforcement-learning"}

    inp: Input[list[str]] = Input()
    """SMILES sampled by REINVENT or LibINVENT."""

    out: Output[NDArray[np.float32]] = Output()
    """Best Vina affinity for every input SMILES."""

    receptor: FileParameter[Annotated[Path, Suffix("pdbqt")]] = FileParameter()
    """Prepared rigid receptor in PDBQT format."""

    search_center: Parameter[tuple[float, float, float]] = Parameter()
    """Docking-box center in Angstrom."""

    search_range: Parameter[tuple[float, float, float]] = Parameter(
        default=(20.0, 20.0, 20.0)
    )
    """Docking-box dimensions in Angstrom."""

    scoring: Parameter[Literal["vina", "vinardo"]] = Parameter(default="vina")
    exhaustiveness: Parameter[int] = Parameter(default=8)
    num_modes: Parameter[int] = Parameter(default=1)
    seed: Parameter[int] = Parameter(default=20260928)
    cpu: Parameter[int] = Parameter(default=1)
    ph: Parameter[float] = Parameter(default=7.4)
    ph_tolerance: Parameter[float] = Parameter(default=1.0)
    max_states: Parameter[int] = Parameter(default=1)
    failure_score: Parameter[float] = Parameter(default=0.0)

    def run(self) -> None:
        from vina import Vina

        smiles_batch = self.inp.receive()
        center = tuple(float(value) for value in self.search_center.value)
        box_size = tuple(float(value) for value in self.search_range.value)
        engine = Vina(
            sf_name=self.scoring.value,
            cpu=int(self.cpu.value),
            seed=int(self.seed.value),
            verbosity=0,
        )
        engine.set_receptor(str(self.receptor.filepath))
        engine.compute_vina_maps(center=list(center), box_size=list(box_size))

        ligand_root = Path("ligands").resolve()
        ligand_root.mkdir(exist_ok=True)
        preparation_spec = {
            "ph": float(self.ph.value),
            "ph_tolerance": float(self.ph_tolerance.value),
            "max_states": int(self.max_states.value),
            "seed": int(self.seed.value),
        }
        # Match Maize ExtractScores: the port is typed as float32 for workflow
        # compatibility, while float64 scalars remain JSON-serializable in the
        # REINVENT interceptor.
        scores = np.full(
            len(smiles_batch), float(self.failure_score.value), dtype=np.float64
        )
        for index, raw_smiles in enumerate(smiles_batch):
            smiles = str(raw_smiles).strip()
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            if molecule is None or any(
                atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms()
            ):
                self.logger.warning("Invalid generated SMILES at index %s: %s", index, smiles)
                continue
            input_id = f"RL_{index + 1:06d}"
            row = SimpleNamespace(
                smiles=smiles,
                input_id=input_id,
                source_index=index,
            )
            try:
                states = _prepare_ligand_states(
                    row,
                    preparation_spec,
                    ligand_root / input_id,
                )
                best_score: float | None = None
                for state_path in states:
                    engine.set_ligand_from_file(str(state_path))
                    engine.dock(
                        exhaustiveness=int(self.exhaustiveness.value),
                        n_poses=int(self.num_modes.value),
                    )
                    energies = engine.energies(n_poses=int(self.num_modes.value))
                    if energies.size:
                        score = float(energies[0][0])
                        if best_score is None or score < best_score:
                            best_score = score
                if best_score is not None:
                    scores[index] = best_score
                    self.logger.info(
                        "Vina score %s/%s: %.3f kcal/mol",
                        index + 1,
                        len(smiles_batch),
                        best_score,
                    )
                else:
                    self.logger.warning("Vina returned no pose for SMILES index %s", index)
            except Exception as exc:  # keep one failed molecule from aborting the epoch
                self.logger.warning(
                    "Vina scoring failed for SMILES index %s: %s: %s",
                    index,
                    type(exc).__name__,
                    exc,
                    exc_info=True,
                )
        self.out.send(scores)
