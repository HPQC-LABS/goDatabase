import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from datetime import  datetime
import re

df_5k = pd.read_json('/Users/vandanrevanur/personal/codes/chemistry/goDatabase/data/df_5k.json', orient='split')

xyz_info = df_5k.xyz_pbe_relaxed
smiles_info = df_5k.canonical_smiles

atom_count_info = {}
version = 0
author = 'Vandan Revanur'

def convert_smiles_to_molecular_formula(inchi: str) -> str:
    # try:
    #     mol = Chem.MolFromInchi(inchi)
    #     molecular_formula = rdMolDescriptors.CalcMolFormula(mol)
    # except:
    #     molecular_formula = None




    ids_slash = [m.start() for m in re.finditer('/', inchi)]
    start_idx = ids_slash[0] + 1
    end_idx = ids_slash[1]

    molecular_formula = inchi[start_idx: end_idx]
    return  molecular_formula


def calc_n_atoms(inchi, xyz):
    try:
        mol = Chem.MolFromInchi(inchi)
        atom_count = 0
        mol = Chem.AddHs(mol)
        for a in mol.GetAtoms():
            atom_count += 1

        atom_count_info[smiles] = atom_count

    except:
        atom_count = xyz.count('\n')-2 # 2 corresponds to the 2 lines in the header of the XYZ string
    return atom_count


source_paper = '# Source paper: https://doi.org/10.1038/s41597-020-0385-y (Stuke et al.)'
source_data_host = '# Source data host: https://doi.org/10.14459/2019mp1507656'
source_download_link = '# Source download link: https://dataserv.ub.tum.de/index.php/s/m1507656 , Filename: df_5k.json'
source_method = '# Source method: PBE_TS-vdW'
date_added = f'# Date added: 2 April 2023 ({author})'
date_modified = f'# Date modified: {datetime.today()} ({author})'
common_name = '# common name: None'
smarts = '# Smarts: None'
cas = '# cas: None'
meta_data_comments = f'{source_paper}\n{source_data_host}\n{source_download_link}\n{source_method}\n{date_added}\n{date_modified}\n'

with open("gw5000.txt", "w") as text_file:
    for i,row in df_5k.iterrows():
        smiles = row.canonical_smiles.strip('\n').strip('\t')
        inchi = row.inchi
        molecular_formula = convert_smiles_to_molecular_formula(inchi)
        xyz = row.xyz_pbe_relaxed
        atom_count = calc_n_atoms(inchi, xyz)
        reference_code = f'# Reference code: {row.refcode_csd}'
        name = f'NAME = {molecular_formula}:GW500.v{version}'
        smiles_info = f'# SMILES : {smiles}'
        n_atoms = f'# Number of atoms: {atom_count}'
        inchi_info = f'# InChI = {inchi}'
        molecule_info = f'{name}\n{n_atoms}\n{common_name}\n{inchi}\n{smiles_info}\n{smarts}\n{reference_code}\n{cas}\n{meta_data_comments}\n{xyz}\n'
        text_file.write(molecule_info)
