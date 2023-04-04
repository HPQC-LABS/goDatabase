import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from datetime import  datetime
df_5k = pd.read_json('/Users/vandanrevanur/personal/codes/chemistry/goDatabase/data/df_5k.json', orient='split')

xyz_info = df_5k.xyz_pbe_relaxed
smiles_info = df_5k.canonical_smiles

atom_count_info = {}
version = 0
author = 'Vandan Revanur'

def convert_smiles_to_molecular_formula(input_smiles: str) -> str:
    mol = Chem.MolFromSmiles(input_smiles)
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
    return molecular_formula

with open("gw5000.txt", "w") as text_file:
    for i,row in df_5k.iterrows():
        molecular_formula = None
        atom_count = None
        try:
            mol = Chem.MolFromSmiles(smiles)
            atom_count = 0
            mol = Chem.AddHs(mol)
            for a in mol.GetAtoms():
                atom_count+=1

            atom_count_info[smiles] = atom_count
            molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
        except:
            pass

        smiles = row.canonical_smiles.strip('\n').strip('\t')
        smiles_info = f'# SMILES : {smiles}'

        common_name = '# common name: None'
        smarts = '# Smarts: None'
        cas = '# cas: None'
        xyz = row.xyz_pbe_relaxed
        reference_code = f'# Reference code: {row.refcode_csd}'

        source_paper = '# Source paper: https://doi.org/10.1038/s41597-020-0385-y (Stuke et al.)'
        source_data_host = '# Source data host: https://doi.org/10.14459/2019mp1507656'
        source_download_link = '# Source download link: https://dataserv.ub.tum.de/index.php/s/m1507656 , Filename: df_5k.json'
        source_method = '# Source method: PBE_TS-vdW'
        date_added = f'# Date added: 2 April 2023 ({author})'
        date_modified = f'# Date modified: {datetime.today()} ({author})'

        meta_data_comments = f'{source_paper}\n{source_data_host}\n{source_download_link}\n{source_method}\n{date_added}\n{date_modified}\n'

        name = f'NAME = {molecular_formula}:GW500.v{version}'
        n_atoms = f'# Number of atoms: {atom_count}'
        inchi = f'# InChI = {row.inchi}'
        molecule_info = f'{name}\n{n_atoms}\n{common_name}\n{inchi}\n{smiles_info}\n{smarts}\n{reference_code}\n{cas}\n{meta_data_comments}\n{xyz}\n'

        text_file.write(molecule_info)

    # print(smiles)
    # print(xyz)
# print(atom_count_info)
print(sorted(atom_count_info, key=atom_count_info.get))