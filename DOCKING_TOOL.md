# BEAD `Docking_schrodinger` tool

`Docking_schrodinger` docks molecules from direct SMILES or a CSV file against
an existing Schrodinger Glide grid. It runs entirely from the `bead-py311-test`
environment while invoking the local Schrodinger 2024-1 executables for
LigPrep, Glide, complex export, and interaction diagrams.

Required inputs:

- exactly one of `csv_path` (a CSV containing `SMILES`, `smiles`, or
  `canonical_smiles`) or `smiles` (a direct list of SMILES)
- `molecule_names`: optional labels corresponding to direct `smiles`
- `grid_path`: prepared Glide `.zip` grid containing its receptor

Main optional inputs:

- `smiles_column`: use a non-standard SMILES column
- `precision`: `HTVS`, `SP` (default), or `XP`
- `poses_per_ligand`: Glide poses retained per prepared ligand
- `max_molecules`: cap the number of valid inputs, useful for a quick check
- `background`: defaults to `True`; query with `docking_job_status_tool`

Example Python call:

```python
from bead.tools.docking import mol_docking

result = mol_docking(
    csv_path="/path/to/molecules.csv",
    grid_path="/path/to/glide-grid.zip",
    precision="SP",
)
print(result)
```

Each job is saved below `~/bead_docking_results/YYYY-MM-DD/`. Its `output`
directory contains:

- `docking_results.csv`: ranked Glide scores, interaction counts, source CSV
  fields, and artifact paths
- `complexes/Inter_1.pdb`, etc.: protein-ligand complexes
- `docked_ligands/*.sdf`: selected best docked poses
- `interaction_diagrams/*.png`: Schrodinger 2D interaction diagrams
- `interaction_details.csv` and `.json`: hydrogen bonds, salt bridges, pi
  interaction counts, and protein residues within 4 Å
- `poseviewer_best.mae` and `receptor.pdb`
- `skipped_inputs.csv`: empty, invalid, or attachment-point SMILES

The tool expects a prepared Glide grid. It does not prepare a protein from a
raw PDB file.
