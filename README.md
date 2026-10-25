# BEAD
Implementation of the Paper "[BEAD: a Human-in-the-Loop Agentic framework for explainable and prospective drug discovery](xxx)" (Unpublished) by Yi Yu et al.. BEAD is a novel LLM-based agent that Bridges the gap between Experiment- and Artificial intelligence-based Drug discovery.

<img src="example/Figure 1.png" width="100%" height="100%">

## Installation (Python 3.11)

BEAD and its local ADMET-AI predictor now run in one Python 3.11 environment.
Create the environment and install this repository in editable mode:

```bash
cd BEAD
conda env create -f environment.yml
conda activate bead
python -m pip install -e .
```

The supplied environment uses the CUDA 12.6 builds of PyTorch and torchvision.
It contains the BEAD agent dependencies, RDKit, and ADMET-AI 2.0.1; a second
ADMET environment is not required.

Confirm the unified environment after installation:

```bash
python -c "import sys, admet_ai, rdkit; print(sys.version); print(admet_ai.__version__); print(rdkit.__version__)"
```

### Local ADMET-AI resources

BEAD includes the ADMET-AI 2.0.1 model checkpoints and DrugBank reference
data under `bead/resources/admet_ai`. No remote ADMET prediction API is used.
`bead.tools.tools.ADMET_prediction` uses the official ADMET-AI Python API and
passes the repository-local model and reference-data paths directly. The CSV
CLI remains available as an optional fallback. The output contains SMILES,
human-scaled hepatocyte intrinsic clearance as
`Clearance_Hepatocyte_AZ_L_h_kg`, half-life in hours as `Half_Life_Obach`,
`Bioavailability_Ma`, and `hERG`. The default clearance conversion assumes
120 million hepatocytes/g liver, 21.4 g liver/kg, and an incubation unbound
fraction of 1; it is not an observed systemic-clearance measurement.

### Target preparation

The `Target_preparation` agent tool cleans raw PDB receptors, removes water and
unwanted heterogens, repairs missing atoms, adds hydrogens, and emits prepared
PDB and AutoDock PDBQT files. It prefers Schrödinger Protein Preparation Wizard
and falls back to PDBFixer/OpenMM plus Meeko when Schrödinger is unavailable.
Required cofactors and metals must be supplied through `keep_heterogens`. When
no binding-site center is provided, it can identify a co-crystal ligand and
calculate the box before removing that ligand. It writes Vina box settings and,
with Schrödinger, a ready-to-use Glide grid. See
[`TARGET_PREPARATION_TOOL.md`](TARGET_PREPARATION_TOOL.md) for usage.

### Small-molecule preparation

The `Mol_preparation` agent tool accepts a direct SMILES/InChI string or a
CSV/SMI collection. It standardizes each structure, applies the same ETKDGv3
conformer-generation strategy used by `Similarity_prediction`, optimizes the
conformers with MMFF94 (UFF fallback), adds explicit hydrogens and Gasteiger
partial charges, and exports SDF, PDB, and AutoDock PDBQT files. See
[`MOL_PREPARATION_TOOL.md`](MOL_PREPARATION_TOOL.md) for usage.

## Getting start
BEAD defaults to the OpenAI API model `gpt-4.1`. Export the API key before
starting the agent; the key is read from the environment and is never stored in
the repository:

```bash
export OPENAI_API_KEY="your-api-key"
```

The default OpenAI API endpoint is selected by `langchain-openai`, so no custom
`base_url` is required.

## Scaffold hopping by human-in-the-loop dialogues
1st query:
```python
from bead.app import BEAD
'''Substructure Search'''
BEAD("Subsearch the structures based on Scaffold_1.mol and mol_list.csv.", "thread_003")
```
1st output is:
<img src="example/Figure 2-1.png" width="100%" height="100%">

```bash
Answer_1: The substructure search has been completed, and the results are saved in the file: `/home/user/save/search_result_0.csv`.
```

2nd query
```python
from bead.app import BEAD
'''Sim calculation and Rev_subsearch'''
BEAD("Filter out the structures conforming to the initial scaffold. Then calculate the similarity within the target molecule against 226192.", "thread_003")
```
2nd output is:
<img src="example/Figure 2-2.png" width="100%" height="100%">

```bash
Answer_2: The filtering and similarity calculation for the molecule with ID 226192 has been completed. The results are saved in the file: `/home/user/save/sim_result_0.csv`.
```

3rd query
```python
from bead.app import BEAD

'''Scaffold hopping validation'''
BEAD(
    "Use the two ligands Ref (SMILES: CN(C)c1nc(=O)n(-c2ccccc2)c2cc(Cl)ccc12) "
    "and New (SMILES: Cc1ncccc1C(=O)n1nc(N)c2ccc(Cl)nc21), together with the "
    "protein file /path/to/mat2a.pdb, to simulate their "
    "drug-target interactions. Return the ranked docking results and the Ref "
    "and New protein-ligand PDB complexes.",
    "thread_003",
)
```
3rd output is:
<img src="example/Figure 2-3.png" width="100%" height="100%">

```bash
Answer_3: The drug-target interaction simulation has completed. The response reports the selected docking backend and the generated `docking_results.csv`. The Ref and New complexes are saved as `output/complexes/Inter_1.pdb` and `output/complexes/Inter_2.pdb`, respectively.
```
