

import os
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from IPython.display import display
from PIL import Image
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem import rdMolEnumerator
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdFingerprintGenerator
from bead.utils.sim import Sim3D

def create_unique_filename(base_filename):
    index = 0
    name, ext = os.path.splitext(base_filename)
    filename = base_filename
    while os.path.exists(filename):
        index += 1
        filename = f"{name[:-2]}_{index}{ext}"
    return filename


def resolve_dataset_path(filename):

    aliases = {
        'mol_list.csv': 'mollist.csv',
        'mollist3.csv': 'mollist.csv',
        'target_molecule': 'mollist.csv',
        'target_molecules': 'mollist.csv',
        'target molecule': 'mollist.csv',
        'target molecules': 'mollist.csv',
        'dataset': 'mollist.csv',
        'database': 'mollist.csv',
        'molecule_dataset': 'mollist.csv',
        'molecule database': 'mollist.csv',
    }
    filename = str(filename)
    alias_key = filename.strip().lower()

    if os.path.exists(filename):
        return filename
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(current_dir, '..'))
    resolved_path = os.path.join(project_root, 'source', filename)
    if not os.path.exists(resolved_path) and alias_key in aliases:
        resolved_path = os.path.join(project_root, 'source', aliases[alias_key])
    
    if os.path.exists(resolved_path):
        return resolved_path
        
    raise FileNotFoundError(f"Could not find the dataset at '{filename}' or in the 'source' directory: {resolved_path}")


def align_bundle_coords(bndl):

    for m in bndl:
        Chem.SanitizeMol(m)
    mcs = rdFMCS.FindMCS(bndl,completeRingsOnly=True)
    q = Chem.MolFromSmarts(mcs.smartsString)
    rdDepictor.Compute2DCoords(q)
    for m in bndl:
        rdDepictor.GenerateDepictionMatching2DStructure(m,q)

def draw_mols_grid_with_highlights(
    mols,
    highlight_atom_lists=None,
    legends=None,
    sub_img_size=(300, 250),
    mols_per_row=3,
    max_mols=30,)

    mols = list(mols or [])[:max_mols]
    if not mols:
        return None

    if highlight_atom_lists is None:
        highlight_atom_lists = [[] for _ in mols]
    else:
        highlight_atom_lists = list(highlight_atom_lists)[:len(mols)]
        highlight_atom_lists += [[] for _ in range(len(mols) - len(highlight_atom_lists))]

    if legends is None:
        legends = ["" for _ in mols]
    else:
        legends = [str(item) for item in list(legends)[:len(mols)]]
        legends += ["" for _ in range(len(mols) - len(legends))]

    cell_w, cell_h = sub_img_size
    rows = (len(mols) + mols_per_row - 1) // mols_per_row
    grid = Image.new("RGB", (cell_w * mols_per_row, cell_h * rows), "white")

    for i, (mol, atoms, legend) in enumerate(zip(mols, highlight_atom_lists, legends)):
        if mol is None:
            continue
        atoms = [int(atom) for atom in atoms if 0 <= int(atom) < mol.GetNumAtoms()]
        img = Draw.MolToImage(mol, size=sub_img_size, highlightAtoms=atoms, legend=legend)
        if img.mode != "RGB":
            img = img.convert("RGB")
        x = (i % mols_per_row) * cell_w
        y = (i // mols_per_row) * cell_h
        grid.paste(img, (x, y))

    return grid


def display_and_save_image(img, basename):
    if img is None:
        return None
    try:
        display(img)
    except Exception as exc:
        print(f"Image display failed: {exc}")

    save_folder_path = os.path.join(os.getcwd(), "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
    image_path = create_unique_filename(os.path.join(save_folder_path, basename))
    img.save(image_path)
    print(f"Image has been saved to {image_path}")
    return image_path


def _safe_add_hs(mol, add_coords=True):
    if mol is None:
        return None
    mol.UpdatePropertyCache(strict=False)
    return Chem.AddHs(mol, addCoords=add_coords)


def Substruc_Prepare(mol):
    '''   
    This function performs structural enumeration, adds hydrogens, aligns coordinates, 
    and displays images for all substructures.

    Parameters
    ------------
    mol : str
        Input molecule file path (e.g., .mol file) or SMARTS string.
    
    Returns
    -------
    qry_addHs_bundle : list
        Returns a list of substructure molecules with added hydrogen atoms.
    '''
    
    qry_addHs_bundle = []
    if mol.endswith('.mol'):
        mol = resolve_dataset_path(mol)
        qry = Chem.MolFromMolFile(mol)
        if qry is None:
            print(f"Failed to load molecule file: {mol}")
            return []
    else:
        qry = Chem.MolFromSmiles(mol)
        if qry is None:
            qry = Chem.MolFromSmarts(mol)
        if qry is None:
            print(f"Failed to parse SMILES/SMARTS string: {mol}")
            return []
    
    qry_bundle = rdMolEnumerator.Enumerate(qry)
    
    if len(qry_bundle) > 1:
        
        for molecule in qry_bundle:
            mol1 = _safe_add_hs(molecule, add_coords=True)
            qry_addHs_bundle.append(mol1)
        align_bundle_coords(qry_bundle)
        print("Substructure list:")
        display(Draw.MolsToGridImage(qry_bundle, subImgSize=(200, 200), molsPerRow=len(qry_bundle)))
    
    elif qry_bundle is None or len(qry_bundle) == 0:
        mol1 = _safe_add_hs(qry, add_coords=True)
        qry_addHs_bundle.append(mol1)
        print("This substructure has only one conformation:")
        display(Draw.MolToImage(qry))
    else:
        molecule = qry_bundle[0]
        mol1 = _safe_add_hs(molecule, add_coords=True)
        qry_addHs_bundle.append(mol1)
        print("This substructure has only one conformation:")
        display(Draw.MolToImage(molecule))
    
    return qry_addHs_bundle

def Sub(substructs, mollist):

    mollist = resolve_dataset_path(mollist)
    
    data = pd.read_csv(mollist,index_col=0)
    molecule_list = data['canonical_smiles'].to_list()
    target_list = data['component_synonym'].to_list()
    activity_list = data['standard_value'].to_list()
    id_list = data.index.tolist()

    substruct_list = Substruc_Prepare(substructs)
    
    mols_no_H = []
    mols_addHs = []
    
    molecule_dict = {}

    for idx, molecule in zip(id_list, molecule_list):
        mol1 = Chem.MolFromSmiles(molecule)
        mols_no_H.append(mol1)
        mol2 = Chem.AddHs(mol1, addCoords=True)
        mols_addHs.append(mol2)

        molecule_dict[idx] = mol2
    
    matches_no_H = []
    matches_addHs = []
    matched_ats = []
    legends = []
    match_info = []
    
    for one in substruct_list:
        for aa, (idx,x) in enumerate(molecule_dict.items()):
            match = x.GetSubstructMatch(one)
            if match:
                x_no_H = Chem.RemoveHs(x)
                new_match = [atom.GetIdx() for atom in x_no_H.GetAtoms() if atom.GetIdx() in match]
                matches_no_H.append(x_no_H)
                matches_addHs.append(x)
                matched_ats.append(new_match)
                legends.append(f"{idx} | {data.loc[idx, 'component_synonym']} | activity_value: {data.loc[idx, 'standard_value']}")
            
                match_info.append({
                    'component_synonym': data.loc[idx, 'component_synonym'],
                    'standard_value': data.loc[idx, 'standard_value'],
                    'canonical_smiles': data.loc[idx, 'canonical_smiles'],
                    'id': idx
                })
    
    num_matches = len(matches_no_H)
    print(f"Number of matched molecules after substructure search: {num_matches}")
    if matches_no_H:
        img = draw_mols_grid_with_highlights(
            matches_no_H,
            highlight_atom_lists=matched_ats,
            sub_img_size=(300, 250),
            legends=legends,
            mols_per_row=5,
            max_mols=20,
        )
        if img is not None:
            display_and_save_image(img, "sub_matches_0.png")
    
    output_csv = 'search_result_0.csv'

    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    if match_info:
        df = pd.DataFrame(match_info)
        df.to_csv(file, index=True)
        print(f"Matching results have been saved to {file}")
    
    return file


def Subsearch(substructs, mollist):
    '''  
    Parameters
    ----------
    substructs : list/str
        Configuration SMILES/SMARTS string or file path for the substructure
    mollist : str
        Path to the CSV file containing the list of small molecules
    
    Returns
    -------
    file : str
        The path of the saved output CSV file
    '''
    
    mollist = resolve_dataset_path(mollist)
    
    data = pd.read_csv(mollist,index_col=0)
    molecule_list = data['canonical_smiles'].to_list()
    target_list = data['component_synonym'].to_list()
    activity_list = data['standard_value'].to_list()
    id_list = data.index.tolist()

    substruct_list = Substruc_Prepare(substructs)
    
    mols_no_H = []
    mols_addHs = []

    molecule_dict = {}

    for idx, molecule in zip(id_list, molecule_list):
        mol1 = Chem.MolFromSmiles(molecule)
        mols_no_H.append(mol1)
        mol2 = Chem.AddHs(mol1, addCoords=True)
        mols_addHs.append(mol2)

        molecule_dict[idx] = mol2
    
    matches_no_H = []
    matches_addHs = []
    matched_ats = []
    legends = []
    match_info = []
    
    for one in substruct_list:
        for aa, (idx,x) in enumerate(molecule_dict.items()):
            match = x.GetSubstructMatch(one)
            if match:
                x_no_H = Chem.RemoveHs(x)
                new_match = [atom.GetIdx() for atom in x_no_H.GetAtoms() if atom.GetIdx() in match]
                matches_no_H.append(x_no_H)
                matches_addHs.append(x)
                matched_ats.append(new_match)
                legends.append(f"{idx} | {data.loc[idx, 'component_synonym']} | activity_value: {data.loc[idx, 'standard_value']}")
            
                match_info.append({
                    'component_synonym': data.loc[idx, 'component_synonym'],
                    'standard_value': data.loc[idx, 'standard_value'],
                    'canonical_smiles': data.loc[idx, 'canonical_smiles'],
                    'id': idx
                })
    
    num_matches = len(matches_no_H)
    print(f"Number of matched molecules after substructure search: {num_matches}")
    if matches_no_H:
        img = draw_mols_grid_with_highlights(
            matches_no_H,
            highlight_atom_lists=matched_ats,
            sub_img_size=(300, 250),
            legends=legends,
            mols_per_row=3,
            max_mols=30,
        )
        if img is not None:
            display(img)
        
    output_csv = 'search_result_0.csv'
    script_dir = os.getcwd()
   
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    if match_info:
        df = pd.DataFrame(match_info)
        df.set_index('id', inplace=True)
        df.to_csv(file, index=True)
        print(f"Matching results have been saved to {file}")
    
    return file


def Similarity_prediction(
    idx,
    csv_path=None,
    min_similarity=0.7,
    max_similarity=0.95,
    num_confs=50,
    parallel=False,
    workers=None,
    rdkit_threads=0
):
    idx = int(idx)
    print(f"###The input ID is {idx}.###")
    
    if csv_path is None:
        csv_path = "mollist.csv"
    if not os.path.exists(csv_path):
        csv_path = resolve_dataset_path(csv_path)

    df = pd.read_csv(csv_path, index_col=0)
    df.index = df.index.astype(int)

    reference_df = df
    if idx not in reference_df.index:
        reference_path = resolve_dataset_path("mollist.csv")
        reference_df = pd.read_csv(reference_path, index_col=0)
        reference_df.index = reference_df.index.astype(int)
    if idx not in reference_df.index:
        raise ValueError(f"Reference molecule ID {idx} was not found in {csv_path} or mollist.csv.")

    print("**index1**: ", reference_df.loc[idx, "canonical_smiles"], " ", reference_df.loc[idx, "standard_value"])
    
    target_synonym = reference_df.loc[idx, 'component_synonym']
    target_df = df[df['component_synonym'] == target_synonym].copy()

    print("**target**: ", target_synonym, " candidate_count: ", len(target_df))

    if target_df.empty:
        raise ValueError(f"No molecules found associated with the target {target_synonym}.")

    query_smiles = reference_df.loc[idx, 'canonical_smiles']

    legends_target = []
    similarity_results = []
    legends_target_mols = []

    target_name = reference_df.loc[idx, 'component_synonym']
    activity = reference_df.loc[idx, 'standard_value'] 
    legends_target.append(f"{target_name} | activity_value: {activity}")
    
    print("Calculating 3D shape similarities. This might take a moment...")
    if parallel:
        worker_count = workers or min(os.cpu_count() or 1, len(target_df))
        worker_count = max(1, min(int(worker_count), len(target_df)))
        print(f"Using parallel Sim3D: workers={worker_count}, rdkit_threads={rdkit_threads}, num_confs={num_confs}")
        tasks = [
            (idx2, query_smiles, row['canonical_smiles'], num_confs, rdkit_threads)
            for idx2, row in target_df.iterrows()
        ]
        sim_by_idx = {}
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            for idx2, similarity, error in executor.map(_sim3d_worker, tasks):
                if error:
                    print(f"Failed to calculate 3D similarity for index {idx2}: {error}")
                sim_by_idx[idx2] = similarity
        similarity_results = [sim_by_idx[idx2] for idx2 in target_df.index]
    else:
        for idx2, row in target_df.iterrows():
            target_smiles = row['canonical_smiles']
        
            if target_smiles:
                try:
                    similarity = Sim3D(
                        query_smiles,
                        target_smiles,
                        num_confs=num_confs,
                        rdkit_threads=rdkit_threads
                    )
                    similarity_results.append(similarity)
                except Exception as e:
                    print(f"Failed to calculate 3D similarity for {target_smiles}: {e}")
                    similarity_results.append(0.0)
            else:
                similarity_results.append(0.0)
    
    target_df = target_df.copy()
    target_df['similarity'] = similarity_results
    
    target_df2 = target_df[
        (target_df.index != idx)
        & (target_df['similarity'] > min_similarity)
        & (target_df['similarity'] < max_similarity)
    ].copy()
    target_df2 = target_df2.sort_values(by='similarity', ascending=False)
    
    legends_target_mols = target_df2.apply(lambda row: f"{row.name} | {target_name} | {row['standard_value']} | {row['similarity']:.3f}", axis=1).tolist()
    
    mols_target = [Chem.MolFromSmiles(reference_df.loc[idx, 'canonical_smiles'])]
    img_target = draw_mols_grid_with_highlights(
        mols_target,
        sub_img_size=(300, 250),
        legends=legends_target,
        mols_per_row=1,
        max_mols=1,
    )
    if img_target is not None:
        display(img_target)
    
    mols_target_mols = [Chem.MolFromSmiles(row['canonical_smiles']) for idx, row in target_df2.iterrows()]
    img_target_mols = draw_mols_grid_with_highlights(
        mols_target_mols,
        sub_img_size=(300, 250),
        legends=legends_target_mols,
        mols_per_row=3,
        max_mols=100,
    )
    if img_target_mols is not None:
        display_and_save_image(img_target_mols, "similarity_grid_0.png")

    output_csv = 'sim_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    target_df2.to_csv(file, index=True)

    print(f"Similarity results have been saved in descending order to the file: {file}")

    return file


def Filter_then_similarity(
    reference_id,
    scaffold_file='Scaffold_1.mol',
    min_similarity=0.7,
    max_similarity=0.95,
    num_confs=50,
    parallel=True,
    workers=None,
    rdkit_threads=1
):
    reference_id = int(reference_id)
    source_csv = resolve_dataset_path('mollist.csv')
    df = pd.read_csv(source_csv, index_col=0)
    df.index = df.index.astype(int)

    if reference_id not in df.index:
        raise ValueError(f"Reference molecule ID {reference_id} was not found in {source_csv}.")

    target_name = df.loc[reference_id, 'component_synonym']
    target_df = df[df['component_synonym'] == target_name].copy()

    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    target_csv = os.path.join(save_folder_path, f"{target_name}_target_0.csv")
    target_csv = create_unique_filename(target_csv)
    target_df.to_csv(target_csv, index=True)
    print(f"Target-specific molecules have been saved to {target_csv}")
    print(f"Target {target_name} candidate count before reverse subsearch: {len(target_df)}")

    reverse_csv = Rev_subsearch(scaffold_file, target_csv)
    reverse_df = pd.read_csv(reverse_csv, index_col=0)
    print(f"Candidate count after reverse subsearch: {len(reverse_df)}")

    sim_csv = Similarity_prediction(
        reference_id,
        csv_path=reverse_csv,
        min_similarity=min_similarity,
        max_similarity=max_similarity,
        num_confs=num_confs,
        parallel=parallel,
        workers=workers,
        rdkit_threads=rdkit_threads
    )
    print(f"Filter-then-similarity workflow completed. Results saved to {sim_csv}")
    return sim_csv


def Rev_subsearch(mol, csv_list):

    print(f'###{mol}###')
    print(f'####{csv_list}####')
    
    csv_list = resolve_dataset_path(csv_list)
    
    df1 =  pd.read_csv(csv_list,index_col=0)
    try:
        file = Sub(mol, csv_list)
        df2 = pd.read_csv(file,index_col=0)
        print('subsearch performed')
    except:
        print('subsearch not performed')
        df2 = df1.iloc[0:0] 
    
    df3 = df1[~df1['canonical_smiles'].isin(df2['canonical_smiles'])]
    
    index_name = df3.index.name or 'id'
    df3_reset = df3.reset_index().rename(columns={df3.index.name or 'index': index_name})
    agg = {
        'standard_value': 'mean',
        'component_synonym': 'first',
        index_name: 'first'
    }
    if 'similarity' in df3_reset.columns:
        agg['similarity'] = 'first'
    df4 = df3_reset.groupby('canonical_smiles', as_index=False).agg(agg)
    df5 = df4.set_index(index_name)

    legends_target_mols = df5.apply(lambda row: f"{row.name} | {row['component_synonym']} | {row['standard_value']}", axis=1).tolist()
    mols_target_mols = [Chem.MolFromSmiles(row['canonical_smiles']) for idx, row in df5.iterrows()]
    img_target_mols = draw_mols_grid_with_highlights(
        mols_target_mols,
        sub_img_size=(300, 250),
        legends=legends_target_mols,
        mols_per_row=5,
        max_mols=100,
    )
    if img_target_mols is not None:
        display_and_save_image(img_target_mols, "reverse_result_grid_0.png")

    output_csv = 'reverse_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    
    df5.to_csv(file, index=True)
    
    return file
    

def merge_and_deduplicate(variables_list):
    csv_paths = []
    
    for i, var in enumerate(variables_list, 1):
        csv_path = Subsearch(var, 'mol_list.csv')
        csv_paths.append(csv_path)
        print(f"f{i}: {csv_path}")
    
    dataframes = []
    for path in csv_paths:
        try:
            df = pd.read_csv(path)
            dataframes.append(df)
            print(f"Successfully read: {path}")
        except Exception as e:
            print(f"Failed to read file {path}: {e}")
    
    if not dataframes:
        raise ValueError("No CSV files were successfully read")
    
    merged_df = pd.concat(dataframes, ignore_index=True)
    print(f"Total rows before merging: {len(merged_df)}")
    
    deduplicated_df = merged_df.drop_duplicates(subset=['id'], keep='first')
    print(f"Rows after deduplication: {len(deduplicated_df)}")
    print(f"Number of duplicate rows deleted: {len(merged_df) - len(deduplicated_df)}")

    output_csv = 'merge_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save1")
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    deduplicated_df.to_csv(file, index=False)
    
    return deduplicated_df
