from pathlib import Path

from pymatgen.core import Structure, Element
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import datetime
import numpy as np
from collections import defaultdict
from typing import List, Optional, Tuple, Dict
from Wigner_Seitz_radius import get_rws, get_rws_physical
from brave_from_pearson import Pearson, international_numbers_to_AP


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


def generate_sys_data(structure):
    """
    Анализирует структуру, возвращает все данные, необходимые
    для генерации .sys и .pot файлов в формате SPR-KKR/xband.

    Возвращает словарь с ключами:
        bravais
        system_name
        prim_structure
        symmetrized (from primitive)
        atom_types  [(label, Z, avg_mag, site_indices), ...]
        site2type   [int,...]  — сопоставление каждого сайта типу
        rws_dict    {elem: rws}
        prim_matrix (3x3)
        a_au        (float)
        a_ang       (float)
        spacegroup  (int)
        magmoms     reduced list of mags
    """

    sga = SpacegroupAnalyzer(structure, symprec=1e-3)
    prim = sga.get_primitive_standard_structure(keep_site_properties=True)
    prim = set_types_gpt(prim)
    conv = sga.get_conventional_standard_structure(keep_site_properties=True)
    sga_prim = SpacegroupAnalyzer(prim)
    symmetrized = sga_prim.get_symmetrized_structure()

    pearson = sga.get_pearson_symbol()[:2]

    br = [str(i) for i in Pearson.from_symbol(pearson)]
    brave = f'{br[1]:>13}        {br[2]:<12}{br[3]:<15}{br[4]:<7}{br[5]:<6}'

    # 'BRAVAIS            9        tetragonal  body-centered  4/mmm  D_4h'
    # 'BRAVAIS           13        cubic       face-centered  m3m    O_h '

    nat = Poscar(prim).natoms
    # print((prim.composition.elements, nat))
    system_name = "".join([f'{el.symbol}{n}' if n != 1 else f'{el.symbol}' for el, n in zip(prim.composition.elements, nat)])
    # print(system_name)
    # решётка
    a_ang = conv.lattice.a
    ANG_TO_AU = 1.889726125
    a_au = a_ang * ANG_TO_AU
    prim_matrix = prim.lattice.matrix / a_ang

    # радиусы Вигнера–Зейтца
    rws_dict = get_rws(prim)

    rmt_class, rws_class = get_rws_physical(
        prim,
        symmetrized
    )

    return {
        "bravais": brave,
        "system_name": system_name,
        "prim_structure": prim,
        "conv_structure": conv,
        "symmetrized": symmetrized,
        "rws_dict": rws_dict,
        "prim_matrix": prim_matrix,
        "a_au": a_au,
        "a_ang": a_ang,
        "spacegroup": (sga.get_space_group_number(), international_numbers_to_AP[sga.get_space_group_number()]),
        "rmt_class": rmt_class,
        "rws_class": rws_class
    }


if __name__ == '__main__':
    p = Path('SPR_KKR/Al/Tsharp/CONTCAR')
    s = Structure.from_file(p)

    sga = SpacegroupAnalyzer(s)
    new = sga.get_refined_structure()
    print(f'Old:\n{s.frac_coords}')
    print(f'Refined:\n{new.frac_coords}')

    data = generate_sys_data(s)
    print(data['rws_class'])
