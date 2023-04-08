import pandas as pd
import re
from mendeleev import element
from collections import Counter
from tqdm import tqdm
from dateutil import tz
from datetime import datetime


s_time = datetime.now(tz=tz.gettz('Europe/Stockholm'))
print(f'Started at time : {s_time}')

df_5k = pd.read_json('/Users/vandanrevanur/personal/codes/chemistry/goDatabase/data/gw5000/df_5k.json', orient='split')

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


def calc_custom_molecular_formula(xyz):
    atoms = []
    for line in xyz.split('\n'):
        str_regex = '[a-zA-Z]+'
        atom_name = re.findall(str_regex, line)

        atoms.extend(atom_name)

    counter = Counter(atoms)

    atom_names = list(counter.keys())
    n_atoms = list(counter.values())

    molecular_formula_info = []
    for atom, n in zip(atom_names, n_atoms):
        molecular_formula_info.append(atom)
        molecular_formula_info.append(str(n))

    custom_molecular_formula = ''.join(molecular_formula_info)
    return custom_molecular_formula


def calculate_n_electrons_in_mol(xyz):
    molecular_formula = calc_custom_molecular_formula(xyz)
    # corner case: C20H19AsCl2P

    n_electrons = 0
    num_regex = "[0-9]+"
    pos_of_numbers = [
        sub_str.start() for sub_str in re.finditer(num_regex, molecular_formula)
    ]
    str_regex = '[a-zA-Z]+'
    pos_of_atomic_symbols = [
        sub_str.start() for sub_str in re.finditer(str_regex, molecular_formula)
    ]

    for pos_atom_symbol, pos_n_atoms in zip(pos_of_atomic_symbols, pos_of_numbers):
        atom_symbol = molecular_formula[pos_atom_symbol:pos_n_atoms]
        n_atoms = int(molecular_formula[pos_n_atoms])
        atom = element(atom_symbol)
        atom_number = atom.atomic_number
        n_electrons += atom_number*n_atoms

    return n_electrons

def remove_xyz_header(xyz):
    end_idx_header = [m.start() for m in re.finditer('\n', xyz)][2]
    xyz_header_removed = xyz[end_idx_header+1:]
    return xyz_header_removed

source_paper = '# Source paper: https://doi.org/10.1038/s41597-020-0385-y (Stuke et al.)'
source_data_host = '# Source data host: https://doi.org/10.14459/2019mp1507656'
source_download_link = '# Source download link: https://dataserv.ub.tum.de/index.php/s/m1507656 , Filename: df_5k.json'
source_method = '# Source method: PBE_TS-vdW'
date_added = f'# Date added: 2023-02-05  14:14:48 ({author})'
date_modified = f'# Date modified: {datetime.today()} ({author})'
common_name = '# common name: Unknown'
smarts = '# Smarts: Unknown'
cas = '# cas: Unknown'
spin_multiplicity = 'Singlet'
electronic_state = 'Ground'
charge = 0
point_group = 'Unknown'
spm_info = f'Spin multiplicity = {spin_multiplicity}'
electron_state_info = f'Electronic state = {electronic_state}'
charge_info = f'Charge = {charge}'
point_group_info = 'Point group = Unknown'


with open("gw5000.txt", "w") as text_file:
    for i,row in tqdm(df_5k.iterrows(), desc='Generating gw5k dataset: ', total=df_5k.shape[0]):
        smiles = row.canonical_smiles.strip('\n').strip('\t')
        inchi = row.inchi
        xyz = row.xyz_pbe_relaxed
        xyz = remove_xyz_header(xyz)

        molecular_formula = convert_smiles_to_molecular_formula(inchi)
        n_electrons_in_mol = calculate_n_electrons_in_mol(xyz)
        n_elec_info = f'Number of electrons = {n_electrons_in_mol}'
        properties = f'# Properties: {n_elec_info}, {spm_info}, {electron_state_info}, {charge_info}, {point_group_info}'
        meta_data_comments = f'{source_paper}\n{source_data_host}\n{source_download_link}\n{source_method}\n{properties}\n{date_added}\n{date_modified}\n'

        atom_count = row.number_of_atoms
        # atom_count = calc_n_atoms(inchi, xyz)
        reference_code = f'# Reference code: {row.refcode_csd}'
        name = f'NAME = {molecular_formula}:GW5000.v{version}'
        smiles_info = f'# SMILES : {smiles}'
        n_atoms = f'# Number of atoms: {atom_count}'
        inchi_info = f'# InChI = {inchi}'
        molecule_info = f'{name}\n{n_atoms}\n{common_name}\n{inchi_info}\n{smiles_info}\n\n{smarts}\n\n{reference_code}\n{cas}\n{meta_data_comments}\n{xyz}\n'
        text_file.write(molecule_info)


e_time = datetime.now(tz=tz.gettz('Europe/Stockholm'))
print(f'Endded at time : {e_time}')
print(f'time diff: {e_time - s_time}')