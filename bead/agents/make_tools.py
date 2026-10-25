


from langchain_core.tools import tool

# Note: Import tools.py from the root-level tools package
from bead.tools import tools as core_tools
from bead.tools.autodock_docking import autodock_docking, autodock_job_status
from bead.tools.docking import docking_job_status, mol_docking
from bead.tools.generation import generation_job_status, start_generation_job
from bead.tools.molecule_preparation import (
    molecule_preparation,
    molecule_preparation_job_status,
)
from bead.tools.target_preparation import (
    target_preparation,
    target_preparation_job_status,
)

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
def Subsearch(query_smiles: str, dataset_path: str) -> str:
    """Use this tool to identify molecules containing a specific chemical core within a dataset.
    Input: query_smiles (the core scaffold SMILES), dataset_path (path to the database)."""
    return core_tools.Subsearch(_normalize_scaffold(query_smiles), dataset_path)

@tool
def Rev_subsearch(toxic_smiles: str, dataset_path: str) -> str:
    """Crucial Negative Filter Tool. Use this tool to perform a reverse substructure search 
    to proactively eliminate molecules containing user-specified motifs."""
    return core_tools.Rev_subsearch(_normalize_scaffold(toxic_smiles), dataset_path)

@tool
def Similarity_prediction(reference_id: int) -> str:
    """Use this tool to calculate 3D structural similarities against a reference molecule.
    Input: reference_id, the molecule ID/index to compare against, such as 226192."""
    return core_tools.Similarity_prediction(reference_id)


@tool
def ADMET_prediction(
    csv_path: str,
    smiles_column: str | None = None,
    output_csv: str | None = None,
    batch_size: int = 100,
    hepatocellularity_million_cells_per_g: float = 120.0,
    liver_weight_g_per_kg: float = 21.4,
    fraction_unbound_in_incubation: float = 1.0,
) -> str:
    """Predict selected ADMET endpoints locally with ADMET-AI.

    The input must be a CSV containing a SMILES or canonical_smiles column.
    The output contains human-scaled hepatocyte intrinsic clearance in L/h/kg,
    half-life in hours from BEAD's dedicated half-life ensemble, binary
    bioavailability probability, and hERG-blocker probability. The clearance
    defaults assume 120 million hepatocytes/g liver,
    21.4 g liver/kg, and an incubation unbound fraction of 1. Predictions are
    statistical estimates and are not experimental measurements. A local
    residual calibration is applied around registered priority clearance
    measurements and is zero below Morgan similarity 0.35.
    """
    return core_tools.ADMET_prediction(
        csv_path=csv_path,
        smiles_column=smiles_column,
        output_csv=output_csv,
        batch_size=batch_size,
        hepatocellularity_million_cells_per_g=hepatocellularity_million_cells_per_g,
        liver_weight_g_per_kg=liver_weight_g_per_kg,
        fraction_unbound_in_incubation=fraction_unbound_in_incubation,
    )


@tool
def Mol_preparation(
    molecular_input: str,
    structure_column: str | None = None,
    name_column: str | None = None,
    num_confs: int = 50,
    seed: int = 42,
    rdkit_threads: int = 0,
    max_iterations: int = 500,
    standardize: bool = True,
    largest_fragment: bool = True,
    generate_pdbqt: bool = True,
    max_molecules: int | None = None,
    output_root: str | None = None,
    background: bool = True,
) -> str:
    """Standardize SMILES/InChI molecules and generate charged 3D structures.

    molecular_input may be one direct SMILES/InChI string or a CSV/SMI/TXT
    path. The tool follows Similarity_prediction's ETKDGv3 embedding, seed 42,
    and MMFF optimization workflow, selects the lowest-energy conformer, adds
    explicit hydrogens and Gasteiger partial charges, and writes SDF, PDB,
    PDBQT, per-atom charge data, and a preparation report.
    """
    return molecule_preparation(
        molecular_input=molecular_input,
        structure_column=structure_column,
        name_column=name_column,
        num_confs=num_confs,
        seed=seed,
        rdkit_threads=rdkit_threads,
        max_iterations=max_iterations,
        standardize=standardize,
        largest_fragment=largest_fragment,
        generate_pdbqt=generate_pdbqt,
        max_molecules=max_molecules,
        output_root=output_root,
        background=background,
    )


@tool
def molecule_preparation_job_status_tool(job_path: str, log_lines: int = 30) -> str:
    """Read status, artifacts, and recent logs for a Mol_preparation job."""
    return molecule_preparation_job_status(job_path=job_path, log_lines=log_lines)


@tool
def Denovo_molgen(
    grid_path: str | None = None,
    max_epoch: int = 100,
    batch_size: int = 32,
    device: str = "cuda:2",
    output_root: str | None = None,
    backend: str = "auto",
    receptor_path: str | None = None,
    center_x: float | None = None,
    center_y: float | None = None,
    center_z: float | None = None,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    vina_exhaustiveness: int = 8,
    vina_num_modes: int = 1,
    vina_cpu: int = 1,
    vina_scoring: str = "vina",
    vina_seed: int = 20260928,
    ph: float = 7.4,
    ph_tolerance: float = 1.0,
    max_states: int = 1,
    failure_score: float = 0.0,
    score_low: float = -10.0,
    score_high: float = -7.0,
    reverse_score: bool = True,
    check_only: bool = False,
) -> str:
    """Start a REINVENT reinforcement-learning molecular-generation job.

    In auto mode the job uses LigPrep and Glide when their executables and
    licenses are available. Otherwise it uses Molscrub, Meeko, and AutoDock
    Vina score-only feedback, which requires a prepared receptor PDBQT and box
    center. The selected backend remains fixed for the complete RL job.
    """
    return start_generation_job(
        method="reinvent",
        backend=backend,
        grid_path=grid_path,
        max_epoch=max_epoch,
        batch_size=batch_size,
        device=device,
        output_root=output_root,
        receptor_path=receptor_path,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        vina_exhaustiveness=vina_exhaustiveness,
        vina_num_modes=vina_num_modes,
        vina_cpu=vina_cpu,
        vina_scoring=vina_scoring,
        vina_seed=vina_seed,
        ph=ph,
        ph_tolerance=ph_tolerance,
        max_states=max_states,
        failure_score=failure_score,
        score_low=score_low,
        score_high=score_high,
        reverse_score=reverse_score,
        check_only=check_only,
    )


@tool
def Constrained_molgen(
    scaffold_path: str,
    grid_path: str | None = None,
    max_epoch: int = 1500,
    batch_size: int = 32,
    device: str = "cuda:7",
    output_root: str | None = None,
    backend: str = "auto",
    receptor_path: str | None = None,
    center_x: float | None = None,
    center_y: float | None = None,
    center_z: float | None = None,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    vina_exhaustiveness: int = 8,
    vina_num_modes: int = 1,
    vina_cpu: int = 1,
    vina_scoring: str = "vina",
    vina_seed: int = 20260928,
    ph: float = 7.4,
    ph_tolerance: float = 1.0,
    max_states: int = 1,
    failure_score: float = 0.0,
    score_low: float = -10.0,
    score_high: float = -7.0,
    reverse_score: bool = True,
    check_only: bool = False,
) -> str:
    """Start a LibINVENT scaffold-decoration reinforcement-learning job.

    The scaffold file must contain REINVENT attachment points. In auto mode the
    job uses LigPrep and Glide when their executables and licenses are available.
    Otherwise it uses Molscrub, Meeko, and AutoDock Vina score-only feedback;
    that fallback requires a prepared receptor PDBQT and box center. The selected
    backend remains fixed for the complete RL job.
    """
    return start_generation_job(
        method="libinvent",
        scaffold_path=scaffold_path,
        grid_path=grid_path,
        max_epoch=max_epoch,
        batch_size=batch_size,
        device=device,
        output_root=output_root,
        backend=backend,
        receptor_path=receptor_path,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        vina_exhaustiveness=vina_exhaustiveness,
        vina_num_modes=vina_num_modes,
        vina_cpu=vina_cpu,
        vina_scoring=vina_scoring,
        vina_seed=vina_seed,
        ph=ph,
        ph_tolerance=ph_tolerance,
        max_states=max_states,
        failure_score=failure_score,
        score_low=score_low,
        score_high=score_high,
        reverse_score=reverse_score,
        check_only=check_only,
    )


@tool
def generation_job_status_tool(job_path: str, log_lines: int = 30) -> str:
    """Read the status, artifacts, and recent log output of a generation job."""
    return generation_job_status(job_path=job_path, log_lines=log_lines)


@tool
def Docking_schrodinger(
    grid_path: str,
    csv_path: str | None = None,
    smiles: list[str] | None = None,
    molecule_names: list[str] | None = None,
    smiles_column: str | None = None,
    precision: str = "SP",
    poses_per_ligand: int = 1,
    max_molecules: int | None = None,
    output_root: str | None = None,
    background: bool = True,
) -> str:
    """Dock direct SMILES or a SMILES CSV against a prepared Glide grid.

    The job runs LigPrep followed by Glide and writes a ranked CSV, docked ligand
    SDF files, protein-ligand Inter_*.pdb complexes, interaction detail CSV/JSON,
    and 2D interaction diagrams. Supply exactly one of csv_path or smiles; use
    molecule_names to label direct SMILES. Attachment-point or invalid SMILES
    are skipped and recorded. Use docking_job_status_tool for a background job.
    """
    return mol_docking(
        csv_path=csv_path,
        grid_path=grid_path,
        smiles=smiles,
        molecule_names=molecule_names,
        smiles_column=smiles_column,
        precision=precision,
        poses_per_ligand=poses_per_ligand,
        max_molecules=max_molecules,
        output_root=output_root,
        background=background,
    )


@tool
def docking_job_status_tool(job_path: str, log_lines: int = 30) -> str:
    """Read the status, artifacts, and recent log output of a mol_docking job."""
    return docking_job_status(job_path=job_path, log_lines=log_lines)


@tool
def Docking_autodock(
    receptor_path: str,
    center_x: float,
    center_y: float,
    center_z: float,
    csv_path: str | None = None,
    smiles: list[str] | None = None,
    molecule_names: list[str] | None = None,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    smiles_column: str | None = None,
    scoring: str = "vina",
    exhaustiveness: int = 16,
    num_modes: int = 9,
    seed: int = 20260928,
    cpu: int = 1,
    max_molecules: int | None = None,
    output_root: str | None = None,
    background: bool = True,
) -> str:
    """Dock direct SMILES or a SMILES CSV with Vina against a receptor.

    Supply exactly one of csv_path or smiles; use molecule_names to label direct
    SMILES. Supply the binding-site box center and dimensions in Angstrom. The
    tool prepares ligands with Molscrub and Meeko and writes ranked affinities,
    docked SDF/PDBQT poses, protein-ligand PDB complexes, interaction details,
    2D interaction diagrams, and 3D pocket-pose renderings. Invalid SMILES are
    recorded and skipped.
    """
    return autodock_docking(
        csv_path=csv_path,
        receptor_path=receptor_path,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        smiles=smiles,
        molecule_names=molecule_names,
        smiles_column=smiles_column,
        scoring=scoring,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
        cpu=cpu,
        max_molecules=max_molecules,
        output_root=output_root,
        background=background,
    )


@tool
def autodock_job_status_tool(job_path: str, log_lines: int = 30) -> str:
    """Read status, artifacts, and recent logs for a Docking_autodock job."""
    return autodock_job_status(job_path=job_path, log_lines=log_lines)


@tool
def Target_preparation(
    pdb_path: str,
    backend: str = "auto",
    ph: float = 7.4,
    ph_tolerance: float = 1.0,
    remove_waters: bool = True,
    remove_heterogens: bool = True,
    keep_heterogens: str | None = None,
    fill_missing_sidechains: bool = True,
    fill_missing_loops: bool = False,
    minimize: bool = True,
    minimize_rmsd: float = 0.3,
    center_x: float | None = None,
    center_y: float | None = None,
    center_z: float | None = None,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    auto_box: bool = True,
    ligand_resname: str | None = None,
    ligand_chain: str | None = None,
    ligand_resid: int | None = None,
    box_padding: float = 4.0,
    output_root: str | None = None,
    background: bool = True,
) -> str:
    """Clean and prepare a raw PDB target for molecular docking.

    In auto mode the tool uses Schrodinger Protein Preparation Wizard when it
    is installed, and otherwise falls back to PDBFixer/OpenMM plus Meeko for
    AutoDock Vina. It removes waters and heterogens by default, adds missing
    atoms and hydrogens, and writes prepared PDB and PDBQT files. List required
    cofactors or metals as comma-separated residue names in keep_heterogens.
    By default, when no manual center is supplied, it identifies the largest
    drug-like HETATM residue and calculates a box before removing that ligand.
    Use ligand_resname, ligand_chain, and ligand_resid to select a specific
    co-crystal ligand. The calculated per-axis size is at least size_x/y/z and
    includes box_padding on both sides. It writes a Vina box file and, with the
    Schrodinger backend, a Glide grid archive.
    """
    return target_preparation(
        pdb_path=pdb_path,
        backend=backend,
        ph=ph,
        ph_tolerance=ph_tolerance,
        remove_waters=remove_waters,
        remove_heterogens=remove_heterogens,
        keep_heterogens=keep_heterogens,
        fill_missing_sidechains=fill_missing_sidechains,
        fill_missing_loops=fill_missing_loops,
        minimize=minimize,
        minimize_rmsd=minimize_rmsd,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        auto_box=auto_box,
        ligand_resname=ligand_resname,
        ligand_chain=ligand_chain,
        ligand_resid=ligand_resid,
        box_padding=box_padding,
        output_root=output_root,
        background=background,
    )


@tool
def target_preparation_job_status_tool(job_path: str, log_lines: int = 30) -> str:
    """Read status, artifacts, and recent logs for a Target_preparation job."""
    return target_preparation_job_status(job_path=job_path, log_lines=log_lines)


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


# Aggregate all tools for agent.py to invoke
bead_toolset = [
    Subsearch,
    Rev_subsearch,
    Similarity_prediction,
    ADMET_prediction,
    Mol_preparation,
    molecule_preparation_job_status_tool,
    Denovo_molgen,
    Constrained_molgen,
    generation_job_status_tool,
    Docking_schrodinger,
    docking_job_status_tool,
    Docking_autodock,
    autodock_job_status_tool,
    Target_preparation,
    target_preparation_job_status_tool,
    filter_then_similarity_tool,
]
