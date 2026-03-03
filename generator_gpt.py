import os
from collections import defaultdict

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pathlib import Path
from pot_from_gpt import generate_pot_from_data
from sys_from_gpt import generate_sys_from_data
from sys_data_gpt import generate_sys_data
from Inp_files import generate_inps
import numpy as np


def set_types_gpt(structure):
    sga = SpacegroupAnalyzer(structure)
    symm_struct = sga.get_symmetrized_structure()

    equiv_indices = symm_struct.equivalent_indices
    equiv_sites = symm_struct.equivalent_sites

    specie_counter = defaultdict(int)
    specie_total = defaultdict(int)

    # --- сначала считаем сколько групп у каждого элемента ---
    for sites in equiv_sites:
        specie = sites[0].specie.symbol
        specie_total[specie] += 1

    # --- теперь идём строго по equiv_groups ---
    for type_idx, (inds, sites) in enumerate(zip(equiv_indices, equiv_sites)):

        specie = sites[0].specie.symbol
        specie_counter[specie] += 1

        if specie_total[specie] == 1:
            type = specie
        else:
            type = f"{specie}_{specie_counter[specie]}"

        for idx in inds:
            structure[idx].properties["type"] = type
            structure[idx].properties["type_idx"] = type_idx + 1

    return structure


def set_magmoms(structure: Structure, out: Outcar):
    magmoms = [mag['tot'] for mag in out.magnetization]

    for i, _ in enumerate(structure.sites):
        if abs(magmoms[i]) < 0.1:
            structure.sites[i].properties['spin'] = 0
        else:
            structure.sites[i].properties['spin'] = round(magmoms[i], 2)
    return structure

to_calc = []


def parse_vasp2spr(path: Path):
    pos = Poscar.from_file(Path(f"{path}/CONTCAR"))
    out = Outcar(path / "OUTCAR")

    structure = set_magmoms(pos.structure, out)

    data = generate_sys_data(structure)
    pot_text = generate_pot_from_data(data)
    with open(f'{path}/{data["system_name"]}.pot', 'w') as f:
        f.write(pot_text)
    # print(pot_text)

    sys_text = generate_sys_from_data(data)
    with open(f'{path}/{data["system_name"]}.sys', 'w') as f:
        f.write(sys_text)
    # print(sys_text)

    magmoms = [data['prim_structure'].sites[i[0]].spin for i in data['symmetrized'].equivalent_indices]
    scf, jxc = generate_inps(data["system_name"], magmoms, path)


    with open(f'{path}/SCF.inp', 'w') as f:
        f.write(scf)

    with open(f'{path}/JXC.inp', 'w') as f:
        f.write(jxc)


wd = Path('SPR_KKR/Al/Tsharp')
# wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/Al/Tsharp_new_best_parser')
parse_vasp2spr(wd)
exit()
# alloys = [i for i in os.listdir(wd) if i != 'tested']
# for alloy in alloys:
#     groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
#     for group in groups:
#         if 'In' == alloy and 'XA' == group:
#             continue
#
#         path = wd / alloy / group
#         print(f'{alloy} - {group}', end='')
#
#         parse_vasp2spr(path)
#         print(' - complete')
#
