# BEAD `Docking_autodock` tool

`Docking_autodock` docks molecules from direct SMILES or a SMILES CSV with
AutoDock Vina without using Schrodinger. It accepts a curated receptor in PDB
or prepared PDBQT format and requires the center and dimensions of the
binding-site search box in Angstrom.

The workflow validates SMILES with RDKit, enumerates protonation and tautomeric
states at pH 7.4 +/- 1.0 with Molscrub, generates PDBQT files with Meeko, and
docks every prepared state with the Vina or Vinardo scoring function. The
lowest-affinity pose across the prepared states is retained for each input
molecule. ProLIF annotates hydrogen bonds, ionic, hydrophobic, pi-stacking, and
pi-cation interactions; residues with a heavy atom within 4 Angstrom are also
reported.

Required inputs:

- exactly one of `csv_path` (a CSV containing `SMILES`, `smiles`, or
  `canonical_smiles`) or `smiles` (a direct list of SMILES)
- `molecule_names`: optional labels corresponding to direct `smiles`
- `receptor_path`: curated receptor `.pdb`, `.ent`, or prepared `.pdbqt`
- `center_x`, `center_y`, `center_z`: search-box center in Angstrom

Important optional inputs:

- `size_x`, `size_y`, `size_z`: search-box dimensions; each defaults to 20 A
- `scoring`: `vina` (default) or `vinardo`
- `exhaustiveness`: search effort, default 16
- `num_modes`: maximum poses generated per prepared state, default 9
- `seed`: deterministic Vina and conformer-generation seed
- `cpu`: CPU threads assigned to the Vina calculation
- `max_molecules`: input limit for validation runs
- `background`: start a background job when true

Example:

```python
from bead.tools.autodock_docking import autodock_docking

result = autodock_docking(
    csv_path="/path/to/molecules.csv",
    receptor_path="/path/to/receptor.pdb",
    center_x=21.705,
    center_y=-46.419,
    center_z=-9.836,
    size_x=20,
    size_y=24,
    size_z=20,
)
print(result)
```

Jobs are stored under
`~/bead_autodock_results/YYYY-MM-DD/<job-id>/` by default. The output directory
contains `docking_results.csv`, `vina_raw_results.csv`, the prepared receptor,
docked SDF and PDBQT poses, protein-ligand complexes, interaction CSV/JSON
reports, LigPlot-style 2D interaction PNGs under `interaction_diagrams/`, and
static 3D pocket-pose PNGs under `binding_3d/`. Use
`autodock_job_status_tool` to query a background job.

The automatic PDB-to-PDBQT conversion cannot decide biological choices such as
missing-loop reconstruction, relevant waters, cofactors, metals, alternate
locations, or residue protonation. Those choices must be resolved before the
tool is called. Vina affinities are empirical ranking scores and should not be
compared numerically with GlideScore.
