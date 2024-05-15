import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from pathlib import Path

from rdkit import Chem
from chembl_structure_pipeline import standardizer
import rdkit.Chem.MolStandardize.rdMolStandardize as rdMolStandardize
from rdkit.Chem import PandasTools
from rdkit.Chem import inchi as rd_inchi
from molvs import standardize_smiles
from molvs import Standardizer

def string_to_int(s):
    mapping = {"Active": 1, "Inactive": 0}
    return mapping.get(s, None)

pathlist = Path(r"D:\OneDrive\Documentos\LabMol\IC-Citotoxicidade\DeepCytosafe\THP1\RAW").glob('**/*.csv')
savepath = r"D:\OneDrive\Documentos\LabMol\IC-Citotoxicidade\DeepCytosafe"

for path in tqdm(pathlist, desc="Processing files"):
    path_name = path.name
    path_str = str(path)
    
    df  = pd.read_csv(path_str)

    df = df.loc[:, ['PUBCHEM_EXT_DATASOURCE_SMILES', 'PUBCHEM_ACTIVITY_OUTCOME', 'PUBCHEM_CID']]
    df = df.rename(columns={'PUBCHEM_EXT_DATASOURCE_SMILES':'Molecule', 'PUBCHEM_ACTIVITY_OUTCOME':'Outcome'})

    df = df.dropna(subset=['Molecule'])
    df = df.dropna(subset=['Outcome'])
    
    df['Outcome'] = df['Outcome'].apply(string_to_int)

    # Standardize the SMILES
    df['final_smiles'] = [Chem.MolToSmiles(Chem.MolFromMolBlock(standardizer.standardize_molblock(Chem.MolToMolBlock(Chem.MolFromSmiles(smile, sanitize=True))))) for smile in df['Molecule']]
    df = df.reset_index(drop=True)

    df = df.dropna(subset=['final_smiles'])
    df = df.loc[:, ['final_smiles', 'Outcome', 'PUBCHEM_CID']]

    # Calculate the InChI
    inchi_list = []
    for smiles in df['final_smiles']:
        mol = Chem.MolFromSmiles(smiles)
        inchi = Chem.inchi.MolToInchi(mol)
        inchi_list.append(inchi)

    # Adicionar a coluna de InChI no dataframe
    df['InChI'] = inchi_list

    df_active = df.query('Outcome == 0')
    df_inactive = df.query('Outcome == 1')

    df_active = df_active.drop_duplicates(subset=['InChI'], inplace=False)
    df_inactive = df_inactive.drop_duplicates(subset=['InChI'], inplace=False)

    df_no_dup_concord = pd.concat([df_active, df_inactive], axis=0)

    final_drop_dup = df_no_dup_concord.drop_duplicates(subset=['InChI'], keep=False, inplace=False)

    df = final_drop_dup
    df = df.reset_index(drop=True)

    df.to_csv(os.path.join(savepath, f"Curated_{path_name}"), index=False)
    print(f"File {path_name} curated has been saved.")