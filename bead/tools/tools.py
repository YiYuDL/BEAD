
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from IPython.display import display
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem import rdMolEnumerator
from rdkit.Chem import rdFMCS

# Import the Sim3D function from the utils directory
from utils.sim import Sim3D


def create_unique_filename(base_filename):
    index = 0
    name, ext = os.path.splitext(base_filename)
    filename = base_filename
    while os.path.exists(filename):
        index += 1
        # Update filename directly using index instead of adding it in the middle of name
        filename = f"{name[:-2]}_{index}{ext}"
    return filename


def resolve_dataset_path(filename):
    """
    Dynamically resolves the path of a dataset. 
    It first checks if the file exists at the given path (e.g., an absolute path from previous steps).
    If not, it looks in the 'source' folder parallel to the 'tools' directory.
    """
    if os.path.exists(filename):
        return filename
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(current_dir, '..'))
    resolved_path = os.path.join(project_root, 'source', filename)
    
    if os.path.exists(resolved_path):
        return resolved_path
        
    raise FileNotFoundError(f"Could not find the dataset at '{filename}' or in the 'source' directory: {resolved_path}")


def align_bundle_coords(bndl):
    '''
    This function aligns the similar parts (e.g., MCS) of all molecules in the bundle to the same position,
    so that the aligned molecular structure diagrams can more clearly show their differences and commonalities.
    
    Parameters
    ----------
    bndl : list
        A list containing multiple RDKit molecule objects (`Chem.Mol`).
    
    Returns
    ----------
    None
    '''    
    for m in bndl:
        Chem.SanitizeMol(m)
    mcs = rdFMCS.FindMCS(bndl,completeRingsOnly=True)
    q = Chem.MolFromSmarts(mcs.smartsString)
    rdDepictor.Compute2DCoords(q)
    for m in bndl:
        rdDepictor.GenerateDepictionMatching2DStructure(m,q)

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
        qry = Chem.MolFromMolFile(mol)
        if qry is None:
            print(f"Failed to load molecule file: {mol}")
            return []
    else:
        # If it is a SMARTS format substructure, convert to molecule
        qry = Chem.MolFromSmarts(mol)
        if qry is None:
            print(f"Failed to parse SMARTS string: {mol}")
            return []
    
    qry_bundle = rdMolEnumerator.Enumerate(qry)
    
    # If the number of substructures is greater than 1, align coordinates
    if len(qry_bundle) > 1:
        
        for molecule in qry_bundle:
            mol1 = Chem.AddHs(molecule, addCoords=True)
            qry_addHs_bundle.append(mol1)
        align_bundle_coords(qry_bundle)
        print("Substructure list:")
        # Display images of multiple substructures
        display(Draw.MolsToGridImage(qry_bundle, subImgSize=(200, 200), molsPerRow=len(qry_bundle)))
    
    # If the number of substructures is less than or equal to 1
    elif qry_bundle is None or len(qry_bundle) == 0:
        mol1 = Chem.AddHs(qry, addCoords=True)
        qry_addHs_bundle.append(mol1)
        print("This substructure has only one conformation:")
        display(Draw.MolToImage(qry))
    
    return qry_addHs_bundle

def Sub(substructs, mollist):
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
    
    # Dynamically resolve path
    mollist = resolve_dataset_path(mollist)
    
    data = pd.read_csv(mollist,index_col=0)
    molecule_list = data['canonical_smiles'].to_list()
    target_list = data['component_synonym'].to_list()
    activity_list = data['standard_value'].to_list()
    id_list = data.index.tolist()

    substruct_list = Substruc_Prepare(substructs)
    
    mols_no_H = [] # Small molecule list (without H)
    mols_addHs = [] # Small molecule list (with H)
    
    molecule_dict = {}

    # Iterate through molecule_list, process each molecule and save in dictionary
    for idx, molecule in zip(id_list, molecule_list):
        mol1 = Chem.MolFromSmiles(molecule)
        mols_no_H.append(mol1)  # Store molecule without H
        mol2 = Chem.AddHs(mol1, addCoords=True)
        mols_addHs.append(mol2)  # Store molecule with H

        # Store hydrogenated molecules in dictionary using DataFrame index as key
        molecule_dict[idx] = mol2
    
    matches_no_H = [] # Matched molecules (without H)
    matches_addHs = [] # Matched molecules (with H)
    matched_ats = [] # Matched atoms
    legends = [] # Legends for each molecule in the image
    match_info = []  # Store matched molecule information
    
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
            
                # Store matched molecule, target, and activity values into match_info
                match_info.append({
                    'component_synonym': data.loc[idx, 'component_synonym'],
                    'standard_value': data.loc[idx, 'standard_value'],
                    'canonical_smiles': data.loc[idx, 'canonical_smiles'],
                    'id': idx
                })
    
    num_matches = len(matches_no_H)
    print(f"Number of matched molecules after substructure search: {num_matches}")
    if matches_no_H:
        # Display matched molecules: molecules without H, highlighting the substructure.
        display(Draw.MolsToGridImage(matches_no_H, highlightAtomLists=matched_ats, subImgSize=(300, 250), legends=legends, molsPerRow=5, maxMols=20))
    
    # Save the matching results to a CSV file
    output_csv = 'search_result_0.csv'

    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    
    # Check if save folder exists, create if it doesn't
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
    
    # Return the saved file path
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
    
    # Dynamically resolve path
    mollist = resolve_dataset_path(mollist)
    
    data = pd.read_csv(mollist,index_col=0)
    molecule_list = data['canonical_smiles'].to_list()
    target_list = data['component_synonym'].to_list()
    activity_list = data['standard_value'].to_list()
    id_list = data.index.tolist()

    substruct_list = Substruc_Prepare(substructs)
    
    mols_no_H = [] # Small molecule list (without H)
    mols_addHs = [] # Small molecule list (with H)

    molecule_dict = {}

    # Iterate through molecule_list, process each molecule and save in dictionary
    for idx, molecule in zip(id_list, molecule_list):
        mol1 = Chem.MolFromSmiles(molecule)
        mols_no_H.append(mol1)  # Store molecule without H
        mol2 = Chem.AddHs(mol1, addCoords=True)
        mols_addHs.append(mol2)  # Store molecule with H

        # Store hydrogenated molecules in dictionary using DataFrame index as key
        molecule_dict[idx] = mol2
    
    matches_no_H = [] # Matched molecules (without H)
    matches_addHs = [] # Matched molecules (with H)
    matched_ats = [] # Matched atoms
    legends = [] # Legends for each molecule in the image
    match_info = []  # Store matched molecule information
    
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
            
                # Store matched molecule, target, and activity values into match_info
                match_info.append({
                    'component_synonym': data.loc[idx, 'component_synonym'],
                    'standard_value': data.loc[idx, 'standard_value'],
                    'canonical_smiles': data.loc[idx, 'canonical_smiles'],
                    'id': idx
                })
    
    num_matches = len(matches_no_H)
    print(f"Number of matched molecules after substructure search: {num_matches}")
    if matches_no_H:
        # Display matched molecules: molecules without H, highlighting the substructure.
        display(Draw.MolsToGridImage(matches_no_H, 
                                      highlightAtomLists=matched_ats, subImgSize=(300, 250),
                                      legends=legends, molsPerRow=3, maxMols=30),)
        
    output_csv = 'search_result_0.csv'
    script_dir = os.getcwd()
   
    save_folder_path = os.path.join(script_dir, "save")
    # Check if save folder exists, create if it doesn't
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
    
    # Return the saved file path
    return file


def Similarity_prediction(idx):
    """
    Calculates the 3D similarity between the reference molecule and other molecules, 
    displays them, and sorts them in descending order of similarity.

    Parameters:
    idx: the input id for the reference molecule
    """
    idx = int(idx)
    print(f"###The input ID is {idx}.###")
    
    # Dynamically resolve path using the helper function
    csv_path = resolve_dataset_path('mollist3.csv')

    df = pd.read_csv(csv_path, index_col=0)

    print("**index1**: ", df.loc[idx, "canonical_smiles"], " ", df.loc[idx, "standard_value"])
    
    # Determine the target of the reference molecule
    target_synonym = df.loc[idx, 'component_synonym']
    
    target_df = df[df['component_synonym'] == target_synonym]

    print("**index2**: ", target_df.loc[idx, "canonical_smiles"], " ", target_df.loc[idx, "standard_value"])

    if target_df.empty:
        raise ValueError(f"No molecules found associated with the target {target_synonym}.")

    # Select the first molecule of the target as the reference molecule
    query_smiles = target_df.loc[idx, 'canonical_smiles']

    legends_target = []  # Annotations for target molecule (target, activity)
    similarity_results = []
    legends_target_mols = []  # Annotations for other molecules (target, activity, similarity)

    # Calculate annotations for the reference molecule (only display the first molecule)
    target_name = target_df.loc[idx, 'component_synonym']
    activity = target_df.loc[idx, 'standard_value'] 
    legends_target.append(f"{target_name} | activity_value: {activity}")
    
    print("Calculating 3D shape similarities. This might take a moment...")
    for idx2, row in target_df.iterrows():
        target_smiles = row['canonical_smiles']
    
        if target_smiles:
            try:
                # Calculate 3D shape similarity
                similarity = Sim3D(query_smiles, target_smiles)
                similarity_results.append(similarity)
            except Exception as e:
                # Fault tolerance in case of 3D conformation generation failure
                print(f"Failed to calculate 3D similarity for {target_smiles}: {e}")
                similarity_results.append(0.0)
        else:
            similarity_results.append(0.0)
    
    target_df['similarity'] = similarity_results
    
    target_df2 = target_df[(target_df['similarity'] > 0.5) & (target_df['similarity'] < 0.95)].copy()
    target_df2 = target_df2.sort_values(by='similarity', ascending=False)
    
    legends_target_mols = target_df2.apply(lambda row: f"{row.name} | {target_name} | {row['standard_value']} | {row['similarity']:.3f}", axis=1).tolist()
    
    # Display only the reference molecule
    mols_target = [Chem.MolFromSmiles(target_df.loc[idx, 'canonical_smiles'])]
    img_target = Draw.MolsToGridImage(
        mols_target, 
        highlightAtomLists=[[]]*len(mols_target),  # If highlighting atoms is not needed, set to an empty list
        subImgSize=(300, 250), 
        legends=legends_target, 
        molsPerRow=1, 
        maxMols=1
    )
    display(img_target)
    
    # Display other molecules
    mols_target_mols = [Chem.MolFromSmiles(row['canonical_smiles']) for idx, row in target_df2.iterrows()]
    img_target_mols = Draw.MolsToGridImage(
        mols_target_mols, 
        highlightAtomLists=[[]]*len(mols_target_mols),  
        subImgSize=(300, 250), 
        legends=legends_target_mols, 
        molsPerRow=3, 
        maxMols=100
    )
    display(img_target_mols)

    output_csv = 'sim_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    # Check if save folder exists, create if it doesn't
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


def Rev_subsearch(mol, csv_list):
    """
    Reverse search for molecules that do not contain the given substructure.

    Parameters:
    mol (str): Input substructure mol file
    csv_list (str): Input search query file
    """

    print(f'###{mol}###')
    print(f'####{csv_list}####')
    
    # Dynamically resolve path
    csv_list = resolve_dataset_path(csv_list)
    
    df1 =  pd.read_csv(csv_list,index_col=0) # Subsearch output result
    try:
        file = Sub(mol, csv_list)
        df2 = pd.read_csv(file,index_col=0)
        print('subsearch performed')
    except:
        print('subsearch not performed')
        df2 = df1.iloc[0:0] 
    
    df3 = df1[~df1['canonical_smiles'].isin(df2['canonical_smiles'])]
    
    df3_reset = df3.reset_index()
    df4 = df3_reset.groupby('canonical_smiles', as_index=False).agg({
        'standard_value': 'mean',  # Take the average of the 'standard_value' column
        'similarity': 'first',  # Take the first value of the 'similarity' column
        'component_synonym': 'first',  # Take the first value of the 'component_synonym' column
        'index': 'first'
    })
    df5 = df4.set_index('index')

    legends_target_mols = df5.apply(lambda row: f"{row.name} | {row['component_synonym']} | {row['standard_value']} | {row['similarity']:.3f}", axis=1).tolist()
    mols_target_mols = [Chem.MolFromSmiles(row['canonical_smiles']) for idx, row in df5.iterrows()]
    img_target_mols = Draw.MolsToGridImage(
        mols_target_mols, 
        highlightAtomLists=[[]]*len(mols_target_mols),  
        subImgSize=(300, 250), 
        legends=legends_target_mols, 
        molsPerRow=5, 
        maxMols=100
    )
    display(img_target_mols)

    output_csv = 'reverse_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save")
    # Check if save folder exists, create if it doesn't
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
    """
    Pass multiple variables into the Subsearch function, merge the returned CSV files, and deduplicate.
    
    Args:
        variables_list: A list containing multiple variables
    
    Returns:
        pd.DataFrame: Merged and deduplicated DataFrame
    """
    # Store all CSV file paths
    csv_paths = []
    
    # Iterate through the variables list and execute the Subsearch function
    for i, var in enumerate(variables_list, 1):
        # Call Subsearch function to get the CSV path. 
        # 'mol_list.csv' will now be dynamically resolved to the source folder by the helper function.
        csv_path = Subsearch(var, 'mol_list.csv')
        csv_paths.append(csv_path)
        print(f"f{i}: {csv_path}")
    
    # Read and merge all CSV files
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
    
    # Merge all DataFrames
    merged_df = pd.concat(dataframes, ignore_index=True)
    print(f"Total rows before merging: {len(merged_df)}")
    
    # Remove duplicate rows based on the 'id' column, keeping the first occurrence
    deduplicated_df = merged_df.drop_duplicates(subset=['id'], keep='first')
    print(f"Rows after deduplication: {len(deduplicated_df)}")
    print(f"Number of duplicate rows deleted: {len(merged_df) - len(deduplicated_df)}")

    output_csv = 'merge_result_0.csv'
    
    script_dir = os.getcwd()
    save_folder_path = os.path.join(script_dir, "save1")
    # Check if save folder exists, create if it doesn't
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
        print(f"Successfully created folder: {save_folder_path}")
    else:
        print(f"Folder already exists: {save_folder_path}")

    file = os.path.join(save_folder_path, output_csv)
    file = create_unique_filename(file)
    deduplicated_df.to_csv(file, index=False)
    
    return deduplicated_df
