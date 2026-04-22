


from langchain_core.tools import tool

# Note: Import tools.py from the root-level tools package
from tools import tools as core_tools

@tool
def subsearch_tool(query_smiles: str, dataset_path: str) -> str:
    """Use this tool to identify molecules containing a specific chemical core within a dataset.
    Input: query_smiles (the core scaffold SMILES), dataset_path (path to the database)."""
    return core_tools.Subsearch(query_smiles, dataset_path)

@tool
def rev_subsearch_tool(toxic_smiles: str, dataset_path: str) -> str:
    """Crucial Negative Filter Tool. Use this tool to perform a reverse substructure search 
    to proactively eliminate molecules containing toxic motifs (e.g., hERG liabilities)."""
    return core_tools.Rev_subsearch(toxic_smiles, dataset_path)

@tool
def similarity_prediction_tool(target_smiles: str, reference_smiles: str) -> str:
    """Use this tool to calculate the 3D structural similarity between a newly designed molecule 
    and a reference compound. Crucial for explainable scaffold hopping."""
    return core_tools.Similarity_prediction(target_smiles, reference_smiles)

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
    docking_autodock_tool,
    target_preparation_tool,
    mol_gen_tool
]
