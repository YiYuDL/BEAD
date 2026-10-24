# Small-molecule preparation tool

`Mol_preparation` standardizes 1D molecular representations and creates
charged three-dimensional structures for subsequent modeling tasks. The input
may be a direct SMILES/InChI string or a CSV, SMI, SMILES, or TXT file.

The workflow uses RDKit MolStandardize, adds explicit hydrogens, generates an
ensemble with ETKDGv3, and optimizes every conformer with MMFF94 using the same
default random seed (`42`) and iteration limit (`500`) as BEAD's
`Similarity_prediction` calculation. UFF is used only when MMFF94 parameters
are unavailable. The lowest-energy conformer is selected, and Gasteiger
partial charges are assigned.

Outputs include:

- `prepared_molecules.sdf`, containing the selected structures;
- individual optimized conformer ensembles and selected SDF structures;
- individual PDB and AutoDock-compatible PDBQT files;
- `partial_charges.csv`, containing per-atom partial charges;
- `preparation_results.csv` and `preparation_report.json`.

Example for one molecule:

```python
from bead.agents.make_tools import Mol_preparation

result = Mol_preparation.invoke({
    "molecular_input": "CC(=O)Oc1ccccc1C(=O)O",
    "background": False,
})
```

Example for a CSV collection:

```python
result = Mol_preparation.invoke({
    "molecular_input": "/path/to/molecules.csv",
    "structure_column": "smiles",
    "name_column": "name",
})
```

Background jobs can be inspected with
`molecule_preparation_job_status_tool`.
