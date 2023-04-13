import pandas as pd
import re
from collections import Counter
from tqdm import tqdm
from dateutil import tz
from datetime import datetime
from dataclasses import dataclass
import os
def convert_inchi_to_molecular_formula(inchi: str) -> str:
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

    atom_number_map = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
    'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25,'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33,
    'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44,
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
    'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77,
    'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88,
    'Ac': 89, 'Th': 90, 'Pa': 91,'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99,
    'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
    'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

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
        atom_number = atom_number_map[atom_symbol]
        n_electrons += atom_number*n_atoms

    return n_electrons

def remove_xyz_header(xyz):
    end_idx_header = [m.start() for m in re.finditer('\n', xyz)][2]
    xyz_header_removed = xyz[end_idx_header+1:]
    return xyz_header_removed


def write_geom_database(geom_info):
    mol_info_filepath = os.path.join(geom_info.mol_info_dir_path,f'df_{geom_info.number_of_molecules//1000}k.json')
    output_geom_database_filepath =  f'gw{geom_info.number_of_molecules}.txt'
    spin_multiplicity = geom_info.spin_multiplicity
    electronic_state= geom_info.electronic_state
    charge = geom_info.charge
    point_group = geom_info.point_group
    source_paper =  f'# Source paper: {geom_info.source_paper_link} ({geom_info.source_paper_authors})'
    source_data_host = f'# Source data host: {geom_info.source_data_host}'
    source_download_link = f'# Source download link: {geom_info.source_download_link} , Filename: df_{geom_info.number_of_molecules//1000}k.json'
    source_method =  f'# Source method: {geom_info.source_method}'
    version = geom_info.version
    common_name =  f'# Common name: {geom_info.common_name}'
    smarts =  f'# Smarts: {geom_info.smarts}'
    cas = f'# CAS: {geom_info.cas}'
    number_of_molecules = geom_info.number_of_molecules

    date_added: str = f'# Date added: {geom_info.date_of_database_creation} ({geom_info.author_creator})'
    date_modified: str = f'# Date modified: {datetime.today().isoformat(sep=" ", timespec="seconds")} ({geom_info.author_modifier})'

    s_time = datetime.now(tz=tz.gettz('Europe/Stockholm'))

    print(f'Started at time : {s_time.isoformat(sep=" ", timespec="seconds")}')
    df_molecules = pd.read_json(mol_info_filepath,orient='split')

    spm_info = f'Spin multiplicity = {spin_multiplicity}'
    electron_state_info = f'Electronic state = {electronic_state}'
    charge_info = f'Charge = {charge}'
    point_group_info = f'Point group = {point_group}'

    molecule_names = []
    output_info = []

    inchi_col_vals = df_molecules.inchi.values
    smiles_col_vals = df_molecules.canonical_smiles.map(lambda x: x.strip('\n').strip('\t'))
    xyz_col_vals = df_molecules.xyz_pbe_relaxed.map(lambda x: remove_xyz_header(x))
    n_atoms_col_vals = df_molecules.number_of_atoms
    csd_code_col_vals = df_molecules.refcode_csd


    for inchi,smiles, xyz, atom_count, csd_code in tqdm(
            zip(inchi_col_vals,smiles_col_vals,xyz_col_vals, n_atoms_col_vals, csd_code_col_vals),
            desc='Generating gw5k dataset: ', total=df_molecules.shape[0]):

        molecular_formula = convert_inchi_to_molecular_formula(inchi)
        n_electrons_in_mol = calculate_n_electrons_in_mol(xyz)
        n_elec_info = f'Number of electrons = {n_electrons_in_mol}'
        properties = f'# Properties: {n_elec_info}, {spm_info}, {electron_state_info}, {charge_info}, {point_group_info}'
        meta_data_comments = f'{source_paper}\n{source_data_host}\n{source_download_link}\n{source_method}\n{properties}\n{date_added}\n{date_modified}\n'
        reference_code = f'# Reference code: {csd_code}'

        if molecular_formula in molecule_names:
            n_times_molecule = molecule_names.count(molecular_formula)
            name = f'NAME = {molecular_formula}({n_times_molecule + 1}):GW{number_of_molecules}.v{version}'
        else:
            name = f'NAME = {molecular_formula}:GW{number_of_molecules}.v{version}'

        molecule_names.append(molecular_formula)

        smiles_info = f'# SMILES : {smiles}'
        n_atoms = f'# Number of atoms: {atom_count}'
        inchi_info = f'# {inchi}'
        molecule_info = f'{name}\n{n_atoms}\n{common_name}\n{inchi_info}\n{smiles_info}\n\n{smarts}\n\n{reference_code}\n{cas}\n{meta_data_comments}\n{xyz}\n'

        output_info.append(molecule_info)

    with open(output_geom_database_filepath, "w") as text_file:
        for mol_str in output_info:
            text_file.write(mol_str)

    e_time = datetime.now(tz=tz.gettz('Europe/Stockholm'))

    print(f'Endded at time : {e_time.isoformat(sep=" ", timespec="seconds")}')
    print(f'time diff: {e_time - s_time}')


@dataclass
class GeomInfo():
    number_of_molecules : int = 62000
    mol_info_dir_path: str = f'/Users/vandanrevanur/personal/codes/chemistry/goDatabase/data/gw5000'
    version: int = 0
    author_creator: str = 'Vandan Revanur'
    author_modifier: str = 'Vandan Revanur'
    date_of_database_creation: str ='2023-02-05 14:14:48'
    source_paper_authors = 'Stuke et al.'
    source_paper_link: str = 'https://doi.org/10.1038/s41597-020-0385-y'
    source_data_host: str = 'https://doi.org/10.14459/2019mp1507656'
    source_download_link: str = 'https://dataserv.ub.tum.de/index.php/s/m1507656'
    source_method: str = 'PBE_TS-vdW'
    common_name: str = 'Unknown'
    smarts: str = 'Unknown'
    cas: str = 'Unknown'
    spin_multiplicity: str = 'Singlet'
    electronic_state: str = 'Ground'
    charge: int = 0
    point_group: str = 'Unknown'


if __name__ == '__main__':
    geom_info = GeomInfo()
    write_geom_database(geom_info)