# Target preparation tool

`Target_preparation` converts a raw PDB target into receptor files that can be
used by BEAD's two docking tools. In its default `auto` mode, it first uses the
Schrödinger Protein Preparation Wizard. If Schrödinger is absent or preparation
fails, it uses PDBFixer/OpenMM and Meeko.

The workflow removes crystallographic waters and non-protein `HETATM` records,
repairs missing protein atoms, adds hydrogens for the requested pH, and writes:

- `prepared_target.pdb` for structure inspection and downstream conversion;
- `prepared_target.pdbqt` for `Docking_autodock`;
- `prepared_target.maegz` when Schrödinger is used;
- `vina_box.txt` when a binding-site center is provided;
- `glide-grid.zip` when a binding-site center is provided and Schrödinger is used;
- `preparation_report.json`, recording the backend, settings, removals, warnings,
  and output paths.

Required cofactors, prosthetic groups, and metals must be named explicitly with
`keep_heterogens`, for example `"HEM,ZN"`. The program cannot infer whether every
heterogen is biologically required.

Example:

```python
from bead.agents.make_tools import Target_preparation

result = Target_preparation.invoke({
    "pdb_path": "/path/to/raw_target.pdb",
    "keep_heterogens": "HEM,ZN",
    "center_x": 12.3,
    "center_y": -4.5,
    "center_z": 8.1,
    "background": False,
})
```

For a background job, pass `background=True` and inspect the returned job
directory with `target_preparation_job_status_tool`.
