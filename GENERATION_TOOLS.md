# BEAD REINVENT and LibINVENT tools

BEAD exposes three LangChain tools:

- `Denovo_molgen`: starts a de novo REINVENT job with either LigPrep/Glide or
  open-source Molscrub/Meeko/AutoDock Vina scoring.
- `Constrained_molgen`: starts a scaffold-conditioned LibINVENT job with either
  LigPrep/Glide or open-source Molscrub/Meeko/AutoDock Vina scoring.
- `generation_job_status_tool`: reports job state, recent log output, and output artifacts.

Generation jobs are background processes. Each job records `job.json`, `status.json`, `job.log`,
the generated runtime configuration, and result files under
`~/bead_generation_results/YYYY-MM-DD/<job-id>/` by default.

Both tools use the current interpreter (`sys.executable`). The default physical devices are GPU 2
for REINVENT and GPU 7 for LibINVENT. Inside each isolated process the selected GPU is exposed as
`cuda:0`. `batch_size` defaults to 32 and can be lowered for short validation runs.

LibINVENT requires a scaffold file containing valid REINVENT attachment points.
The Schrödinger backend requires a Glide grid. The AutoDock backend requires a
prepared receptor PDBQT plus the binding-site box center and dimensions. Both
generation methods use the bundled `reinvent.prior` and `libinvent.prior` files.

The Python stack runs entirely in `bead-py311-test`. Glide and LigPrep remain external Schrodinger
executables and require a valid Schrodinger license. Their default root is
`/home/kemove/schrodinger2024-1`.

Both `Denovo_molgen` and `Constrained_molgen` default to `backend="auto"`. At job creation each verifies
the Glide and LigPrep executables and checks out `GLIDE_MAIN` and `LIGPREP_MAIN`
with `lictest`. If that preflight fails, it selects AutoDock Vina. The selected
backend is written to `job.json` and `status.json` and remains fixed throughout
the run, so Glide and Vina scores are never mixed within one RL trajectory.

The AutoDock feedback node uses the same ligand-preparation and scoring path as
`Docking_autodock`: Molscrub generates protonation/tautomer states, Meeko writes
PDBQT, and the Vina Python API returns the best affinity in kcal/mol. Rendering
and interaction analysis are omitted inside RL epochs. Failed structures receive
`failure_score` (default `0.0`) without changing the order of the score vector.
The raw feedback is recorded by REINVENT in the `maize (raw)` output column.

Example de novo AutoDock run:

```python
from bead.agents.make_tools import Denovo_molgen

job = Denovo_molgen.invoke({
    "backend": "autodock",
    "receptor_path": "/path/to/prepared_target.pdbqt",
    "center_x": 12.3,
    "center_y": -4.5,
    "center_z": 8.1,
    "size_x": 20.0,
    "size_y": 20.0,
    "size_z": 20.0,
    "max_epoch": 5,
})
```

Example scaffold-constrained AutoDock run:

```python
from bead.agents.make_tools import Constrained_molgen

job = Constrained_molgen.invoke({
    "scaffold_path": "/path/to/scaffolds.smi",
    "backend": "autodock",
    "receptor_path": "/path/to/prepared_target.pdbqt",
    "center_x": 12.3,
    "center_y": -4.5,
    "center_z": 8.1,
    "size_x": 20.0,
    "size_y": 20.0,
    "size_z": 20.0,
    "max_epoch": 5,
})
```

Vina and Glide score distributions are target-dependent. Calibrate `score_low`
and `score_high` for the selected backend; lower docking energies are rewarded
when `reverse_score=True`.

To reconstruct the generation dependencies on this server, run:

```bash
bash /home/kemove/BEAD/scripts/install_generation_dependencies.sh
```
