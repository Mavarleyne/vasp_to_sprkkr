import json
import os, sys
from collections import defaultdict

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pathlib import Path
from Potential import generate_pot_from_data
from System import generate_sys_from_data
from Sys_data import generate_sys_data
from Inp_files import generate_inps
import numpy as np


def print_progress(current, total, bar_length=80):
    """
    Отображает прогресс выполнения в консоли.

    Args:
        current (int): Текущее значение.
        total (int): Общее количество итераций.
        bar_length (int): Длина полосы прогресса в символах.
    """
    fraction = current / total
    filled_length = int(bar_length * fraction)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    percent = fraction * 100
    sys.stdout.write(f'\r|{bar}| {percent:.1f}% ({current}/{total})')
    sys.stdout.flush()

    if current == total:
        print()  # перенос строки после завершения


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
    # print(f'sites: {len(structure.sites)}, mags: {len(magmoms)}\n{structure.formula}')

    for i, _ in enumerate(structure.sites):
        if abs(magmoms[i]) < 0.1:
            structure.sites[i].properties['spin'] = 0
        else:
            structure.sites[i].properties['spin'] = round(magmoms[i], 2)
    return structure

to_calc = []


def parse_vasp2spr(path: Path):
    pos = Poscar.from_file(Path(f"{path}/POSCAR"))
    out = Outcar(path / "OUTCAR")
    path = path / 'sprkkr'
    path.mkdir(exist_ok=True)

    structure = set_magmoms(pos.structure, out)

    data = generate_sys_data(structure)
    pot_text = generate_pot_from_data(data)
    with open(f'{path.as_posix()}/{data["system_name"]}.pot', 'w') as f:
        f.write(pot_text)
    # print(pot_text)

    sys_text = generate_sys_from_data(data)
    with open(f'{path.as_posix()}/{data["system_name"]}.sys', 'w') as f:
        f.write(sys_text)
    # print(sys_text)

    magmoms = [data['prim_structure'].sites[i[0]].spin for i in data['symmetrized'].equivalent_indices]
    scf, jxc = generate_inps(data["system_name"], magmoms, path)

    with open(f'{path.as_posix()}/SCF.inp', 'w') as f:
        f.write(scf)

    with open(f'{path.as_posix()}/JXC.inp', 'w') as f:
        f.write(jxc)


wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/samples_all')
# paths = json.loads((wd / 'after_formation.json').read_text(encoding='utf-8'))
paths = (wd / 'filtration_step3.dat').read_text().split('\n')
print(f"\n{' VASP to SPR-KKR parsing ':=^100}")
for i, path in enumerate(paths):
    parse_vasp2spr(Path(path))
    print_progress(i + 1, len(paths))

# for i, (alloy, _) in enumerate(paths.items()):
#     for group in paths[alloy].keys():
#         for order in paths[alloy][group].keys():
#             path = Path(paths[alloy][group][order])
#             parse_vasp2spr(path)

# els = []
# for i, (alloy, _) in enumerate(paths.items()):
#     for group in paths[alloy].keys():
#         for order in paths[alloy][group].keys():
#             path = Path(paths[alloy][group][order])
#             print(path.as_posix())
#             form = path.parts[-3]
#             form = form.replace('4', ' ').replace('8', ' ')
#             form = form.split()
#             els.extend(form)
# print(set(i for i in els))
# print(len(set(els)))
# wd = Path('SPR_KKR/Al/Tsharp')
# wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/Al/Tsharp_new_best_parser')
# parse_vasp2spr(wd)
# exit()
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

# wd = Path('SPR_KKR')
# for path in wd.rglob('CONTCAR'):
#     if 'In' in path.as_posix() or 'Sn' in path.as_posix():
#         continue
#     parse_vasp2spr(path.parent)

# wd = '/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ'
# for i in ['Al', 'Ga', 'Si', 'Ge']:
#     for j in ['L21', 'XA', 'TC', 'TP', 'Tsharp']:
#         print(f'{wd}/{i}/{j}')


