
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdShapeAlign

def Sim3D(smiles_1: str, smiles_2: str, num_confs: int = 50, rdkit_threads: int = 0) -> float:
    """
    Calculate the highest 3D shape similarity between two molecules (input via SMILES).
    
    Args:
        smiles_1 (str): SMILES of the reference molecule
        smiles_2 (str): SMILES of the probe molecule
        
    Returns:
        float: The highest 3D shape similarity (between 0.0 and 1.0)
    """
    num_confs = int(num_confs)
    rdkit_threads = int(rdkit_threads)
    if num_confs < 1:
        raise ValueError("num_confs must be at least 1")
    if rdkit_threads < 0:
        raise ValueError("rdkit_threads cannot be negative")

    # 1. Parse SMILES before adding hydrogens so invalid input raises a clear error.
    mol_1 = Chem.MolFromSmiles(smiles_1) if isinstance(smiles_1, str) else None
    mol_2 = Chem.MolFromSmiles(smiles_2) if isinstance(smiles_2, str) else None
    if mol_1 is None or mol_2 is None:
        raise ValueError("SMILES parsing failed. Please check the format of smiles_1 or smiles_2.")
    mol_1 = Chem.AddHs(mol_1)
    mol_2 = Chem.AddHs(mol_2)

    # 2. Configure 3D conformation generation parameters (ETKDGv3 algorithm)
    embed_params = AllChem.ETKDGv3()
    embed_params.randomSeed = 42
    embed_params.numThreads = int(rdkit_threads)

    # 3. Generate a single static conformation for the reference molecule and optimize with force field
    # Note: If EmbedMolecule fails, it returns -1
    if AllChem.EmbedMolecule(mol_1, embed_params) == -1:
        return 0.0
    AllChem.MMFFOptimizeMolecule(mol_1, maxIters=500)

    # 4. Generate conformations for the probe molecule and perform parallel force field optimization
    cids = AllChem.EmbedMultipleConfs(mol_2, num_confs, embed_params)
    if not cids:
        return 0.0
    AllChem.MMFFOptimizeMoleculeConfs(mol_2, numThreads=rdkit_threads, maxIters=500)

    best_shape_sim = 0.0

    # 5. Iterate through all conformations of the probe molecule to calculate Shape-it similarity
    for cid in cids:
        # Ignore color (pharmacophore) scoring to accelerate computation.
        shape_sim, _ = rdShapeAlign.AlignMol(
            mol_1,
            mol_2,
            refConfId=0,
            probeConfId=cid,
            useColors=False
        )
        
        # Record the highest similarity
        if shape_sim > best_shape_sim:
            best_shape_sim = shape_sim

    return best_shape_sim
