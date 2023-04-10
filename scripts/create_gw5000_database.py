import pandas as pd
import re
from mendeleev import element
from collections import Counter
from tqdm import tqdm
from dateutil import tz
from datetime import datetime

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


def write_geom_database(mol_info_filepath, out_geom_database_filepath):
    s_time = datetime.now(tz=tz.gettz('Europe/Stockholm'))

    print(f'Started at time : {s_time.isoformat(sep=" ", timespec="seconds")}')
    df_5k = pd.read_json(mol_info_filepath,orient='split')

    version = 0
    author = 'Vandan Revanur'
    source_paper = '# Source paper: https://doi.org/10.1038/s41597-020-0385-y (Stuke et al.)'
    source_data_host = '# Source data host: https://doi.org/10.14459/2019mp1507656'
    source_download_link = '# Source download link: https://dataserv.ub.tum.de/index.php/s/m1507656 , Filename: df_5k.json'
    source_method = '# Source method: PBE_TS-vdW'
    date_added = f'# Date added: 2023-02-05  14:14:48 ({author})'
    date_modified = f'# Date modified: {datetime.today().isoformat(sep=" ", timespec="seconds")} ({author})'
    common_name = '# Common name: Unknown'
    smarts = '# Smarts: Unknown'
    cas = '# CAS: Unknown'
    spin_multiplicity = 'Singlet'
    electronic_state = 'Ground'
    charge = 0
    point_group = 'Unknown'
    spm_info = f'Spin multiplicity = {spin_multiplicity}'
    electron_state_info = f'Electronic state = {electronic_state}'
    charge_info = f'Charge = {charge}'
    point_group_info = f'Point group = {point_group}'

    molecule_names = []
    output_info = []

    inchi_col_vals = df_5k.inchi.values
    smiles_col_vals = df_5k.canonical_smiles.map(lambda x: x.strip('\n').strip('\t'))
    xyz_col_vals = df_5k.xyz_pbe_relaxed.map(lambda x: remove_xyz_header(x))
    n_atoms_col_vals = df_5k.number_of_atoms
    csd_code_col_vals = df_5k.refcode_csd


    for inchi,smiles, xyz, atom_count, csd_code in tqdm(
            zip(inchi_col_vals,smiles_col_vals,xyz_col_vals, n_atoms_col_vals, csd_code_col_vals),
            desc='Generating gw5k dataset: ', total=df_5k.shape[0]):

        molecular_formula = convert_inchi_to_molecular_formula(inchi)
        n_electrons_in_mol = calculate_n_electrons_in_mol(xyz)
        n_elec_info = f'Number of electrons = {n_electrons_in_mol}'
        properties = f'# Properties: {n_elec_info}, {spm_info}, {electron_state_info}, {charge_info}, {point_group_info}'
        meta_data_comments = f'{source_paper}\n{source_data_host}\n{source_download_link}\n{source_method}\n{properties}\n{date_added}\n{date_modified}\n'
        reference_code = f'# Reference code: {csd_code}'

        if molecular_formula in molecule_names:
            n_times_molecule = molecule_names.count(molecular_formula)
            name = f'NAME = {molecular_formula}({n_times_molecule + 1}):GW5000.v{version}'
        else:
            name = f'NAME = {molecular_formula}:GW5000.v{version}'

        molecule_names.append(molecular_formula)

        smiles_info = f'# SMILES : {smiles}'
        n_atoms = f'# Number of atoms: {atom_count}'
        inchi_info = f'# {inchi}'
        molecule_info = f'{name}\n{n_atoms}\n{common_name}\n{inchi_info}\n{smiles_info}\n\n{smarts}\n\n{reference_code}\n{cas}\n{meta_data_comments}\n{xyz}\n'

        output_info.append(molecule_info)

    with open(out_geom_database_filepath, "w") as text_file:
        for mol_str in output_info:
            text_file.write(mol_str)

    e_time = datetime.now(tz=tz.gettz('Europe/Stockholm'))

    print(f'Endded at time : {e_time.isoformat(sep=" ", timespec="seconds")}')
    print(f'time diff: {e_time - s_time}')


if __name__ == '__main__':
    output_geom_database_filepath = "gw5000.txt"
    mol_info_filepath = '/Users/vandanrevanur/personal/codes/chemistry/goDatabase/data/gw5000/df_5k.json'
    write_geom_database(mol_info_filepath, output_geom_database_filepath)