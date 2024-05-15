import pandas as pd
import numpy as np
from tqdm import tqdm

import os
from pathlib import Path
from rdkit import Chem
from chembl_structure_pipeline import standardizer
from rdkit.Chem.MolStandardize.metal import MetalDisconnector
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
    
    df0  = pd.read_csv(path_str)

    df0 = df0.loc[:, ['PUBCHEM_EXT_DATASOURCE_SMILES', 'PUBCHEM_ACTIVITY_OUTCOME', 'PUBCHEM_CID']]
    df0 = df0.rename(columns={'PUBCHEM_EXT_DATASOURCE_SMILES':'Molecule', 'PUBCHEM_ACTIVITY_OUTCOME':'Outcome'})

    df0 = df0.dropna(subset=['Molecule'])
    df0 = df0.dropna(subset=['Outcome'])
    
    df0['Outcome'] = df0['Outcome'].apply(string_to_int)

    rdMol = [Chem.MolFromSmiles(smile, sanitize=True) for smile in df0['Molecule']]
    molBlock = [Chem.MolToMolBlock(mol) for mol in rdMol]
    stdMolBlock = [standardizer.standardize_molblock(mol_block) for mol_block in molBlock]
    molFromMolBlock = [Chem.MolFromMolBlock(std_molblock) for std_molblock in stdMolBlock]
    mol2smiles = [Chem.MolToSmiles(m) for m in molFromMolBlock]
    df0['final_smiles'] = mol2smiles
    df0 = df0.reset_index(drop=True)

    df0 = df0.dropna(subset=['final_smiles'])
    df0 = df0.loc[:, ['final_smiles', 'Outcome', 'PUBCHEM_CID']]

    # Calculate the InChI
    inchi_list = []
    for smiles in df0['final_smiles']:
        mol = Chem.MolFromSmiles(smiles)
        inchi = Chem.inchi.MolToInchi(mol)
        inchi_list.append(inchi)

    # Adicionar a coluna de InChI no dataframe
    df0['InChI'] = inchi_list

    df0_active = df0.query('Outcome == 0')
    df0_inactive = df0.query('Outcome == 1')

    df0_active = df0_active.drop_duplicates(subset=['InChI'], inplace=False)
    df0_inactive = df0_inactive.drop_duplicates(subset=['InChI'], inplace=False)

    df_no_dup_concord = pd.concat([df0_active, df0_inactive], axis=0)

    final_drop_dup = df_no_dup_concord.drop_duplicates(subset=['InChI'], keep=False, inplace=False)

    df0 = final_drop_dup
    df0 = df0.reset_index(drop=True)

    df0.to_csv(os.path.join(savepath, f"Curated_{path_name}"), index=False)
    print(f"File {path_name} curated has been saved.")