

from langchain_core.tools import tool

# Note: Import tools.py from the root-level tools package
from bead.tools import tools as core_tools

def _normalize_scaffold(value: str) -> str:
    scaffold_aliases = {
        "initial scaffold",
        "initial_scaffold",
        "core scaffold",
        "core_scaffold",
        "scaffold",
        "scaffold_1",
    }
    text = str(value).strip()
    if text.lower() in scaffold_aliases:
        return "Scaffold_1.mol"
    return text

@tool
def subsearch_tool(query_smiles: str, dataset_path: str) -> str:
    """Use this tool to identify molecules containing a specific chemical core within a dataset.
    Input: query_smiles (the core scaffold SMILES), dataset_path (path to the database)."""
    return core_tools.Subsearch(_normalize_scaffold(query_smiles), dataset_path)

@tool
def rev_subsearch_tool(toxic_smiles: str, dataset_path: str) -> str:
    """Crucial Negative Filter Tool. Use this tool to perform a reverse substructure search 
    to proactively eliminate molecules containing user-specified motifs."""
    return core_tools.Rev_subsearch(_normalize_scaffold(toxic_smiles), dataset_path)

@tool
def similarity_prediction_tool(reference_id: int) -> str:
    """Use this tool to calculate 3D structural similarities against a reference molecule.
    Input: reference_id, the molecule ID/index to compare against, such as 226192."""
    return core_tools.Similarity_prediction(reference_id)

@tool
def filter_then_similarity_tool(
    reference_id: int,
    scaffold_file: str = "Scaffold_1.mol",
    min_similarity: float = 0.7,
    max_similarity: float = 0.95,
    num_confs: int = 50,
    parallel: bool = True,
    workers: int | None = None,
    rdkit_threads: int = 1
) -> str:
    """Use this tool for the query pattern: filter out structures sharing the initial scaffold,
    then calculate similarity against a reference molecule.
    Workflow: find molecules in mollist.csv with the same component_synonym as reference_id,
    remove molecules containing scaffold_file with reverse substructure search, then calculate
    similarities on that filtered CSV and save rows where min_similarity < similarity < max_similarity.
    num_confs controls the number of generated conformations per candidate; default is 50.
    parallel enables multi-CPU execution across molecules. When parallel is true, use workers to set
    the number of CPU worker processes and keep rdkit_threads at 1 to avoid CPU oversubscription.
    Use this directly for queries like: "Filter out the structures conforming to the initial scaffold.
    Then calculate the similarity within the target molecule against 226192."."""
    return core_tools.Filter_then_similarity(
        reference_id=reference_id,
        scaffold_file=_normalize_scaffold(scaffold_file),
        min_similarity=min_similarity,
        max_similarity=max_similarity,
        num_confs=num_confs,
        parallel=parallel,
        workers=workers,
        rdkit_threads=rdkit_threads
    )

@tool
def docking_autodock_tool(ligand_pdbqt: str, receptor_pdbqt: str) -> str:
    """Use this tool to evaluate protein-ligand binding affinity using AutoDock Vina."""
    return core_tools.Docking_autodock(ligand_pdbqt, receptor_pdbqt)

@tool
def target_preparation_tool(pdb_file: str) -> str:
    """Use this tool to clean and prepare raw PDB files before performing docking simulations."""
    return core_tools.Target_preparation(pdb_file)

@tool
def mol_gen_tool(scaffold_smiles: str, pocket_info: str) -> str:
    """Use this tool for scaffold-constrained structure-based molecular generation."""
    return core_tools.Mol_gen(scaffold_smiles, pocket_info)

# Aggregate all tools for agent.py to invoke
bead_toolset = [
    subsearch_tool,
    rev_subsearch_tool,
    similarity_prediction_tool,
    filter_then_similarity_tool,
    docking_autodock_tool,
    target_preparation_tool,
    mol_gen_tool
]
