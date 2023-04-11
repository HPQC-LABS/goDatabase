#!/usr/bin/env python
# Authors: Vandan Revanur (vandan@hpqc.org)

import requests  # First run "pip install requests"
import psi4
import re

def get_xyz(db_data:str, name:str)-> str:
    p = re.compile(f'NAME = {name}(.*?)NAME', flags=re.DOTALL)
    result = p.search(db_data)
    xyz_info = result.group(1)
    return xyz_info

geometry_database =   'https://raw.githubusercontent.com/HPQC-LABS/goDatabase/master/benchmark_sets/gw5000.txt'
molecule_name =       'C10H15NO2S2'
geometry_identifier = 'GW5000'
geometry_version =    '0'
geometry_name = molecule_name+':'+geometry_identifier+'.v'+geometry_version

page = requests.get(geometry_database)
data = page.text
molecule_xyz = get_xyz(data, geometry_name)

mol_psi4 = psi4.geometry(molecule_xyz)
print(mol_psi4)
