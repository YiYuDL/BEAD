# Target preparation tool

`Target_preparation` converts a raw PDB target into receptor files that can be
used by BEAD's two docking tools. In its default `auto` mode, it first uses the
Schrödinger Protein Preparation Wizard. If Schrödinger is absent or preparation
fails, it uses PDBFixer/OpenMM and Meeko.

Before removing crystallographic waters and non-protein `HETATM` records, the
workflow can calculate a docking box from a co-crystal ligand. It then repairs
missing protein atoms, adds hydrogens for the requested pH, and writes:

- `prepared_target.pdb` for structure inspection and downstream conversion;
- `prepared_target.pdbqt` for `Docking_autodock`;
- `prepared_target.maegz` when Schrödinger is used;
- `vina_box.txt` when a binding-site center is provided or detected;
- `glide-grid.zip` when a binding-site center is provided or detected and
  Schrödinger is used;
- `preparation_report.json`, recording the backend, settings, removals, warnings,
  and output paths.

Required cofactors, prosthetic groups, and metals must be named explicitly with
`keep_heterogens`, for example `"HEM,ZN"`. The program cannot infer whether every
heterogen is biologically required.

Automatic box detection is enabled by default when manual center coordinates
are omitted. Without an explicit selector, the tool selects the largest
drug-like `HETATM` residue in the first PDB model while ignoring water, common
crystallization additives, and small ions. For targets with multiple ligands or
large cofactors, select the intended ligand with `ligand_resname`,
`ligand_chain`, and `ligand_resid`. The center is the ligand heavy-atom centroid.
Each box dimension is the greater of the corresponding `size_x/y/z` value
(20 Å by default) and the ligand span plus `box_padding` on both sides.

Example:

```python
from bead.agents.make_tools import Target_preparation

result = Target_preparation.invoke({
    "pdb_path": "/path/to/raw_target.pdb",
    "keep_heterogens": "HEM,ZN",
    "ligand_resname": "LIG",
    "ligand_chain": "A",
    "ligand_resid": 401,
    "box_padding": 4.0,
    "background": False,
})
```

Manual `center_x`, `center_y`, and `center_z` values take precedence over
automatic detection. Set `auto_box=False` to prepare an apo receptor without a
docking box.

For a background job, pass `background=True` and inspect the returned job
directory with `target_preparation_job_status_tool`.
