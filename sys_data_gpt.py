from pathlib import Path

from pymatgen.core import Structure, Element
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import datetime
import numpy as np
from collections import defaultdict
from typing import List, Optional, Tuple, Dict
from Wigner_Seitz_radius import get_rws
from brave_from_pearson import Pearson, international_numbers_to_AP


def split_atom_types_by_magmom(prim: Structure,
                               magmoms: Optional[List[float]],
                               tol: float = 1e-3
                               ):
    """
    Возвращает:
      - types: [(label, Z, avg_mag, [site_indices])]
      - site2type: mapping site_index -> type_index (1-based)
    Нумерация с суффиксами делается per-element: Fe_1, Fe_2, ...
    Если элемент не разбивается (только 1 кластер), label = 'Fe' (без _1).
    """
    n_sites = len(prim.sites)
    n_magmoms = len(magmoms)

    # --- редукция магнитных моментов ---
    if n_magmoms != n_sites:
        if n_magmoms % n_sites != 0:
            raise ValueError(
                f"Нельзя редуцировать magmoms: {n_magmoms=} не кратно {n_sites=}"
            )
        block = n_magmoms // n_sites
        magmoms = [np.mean(magmoms[i * block:(i + 1) * block]) for i in range(n_sites)]

    # --- группировка по элементам ---
    elem_groups = defaultdict(list)
    for i, site in enumerate(prim.sites):
        elem_groups[site.specie.symbol].append(i)

    types = []
    site2type = [None] * n_sites
    tindex = 1
    elem_counter = defaultdict(int)

    # --- кластеризация по магнитным моментам ---
    for elem, inds in elem_groups.items():
        mm = [magmoms[i] for i in inds]
        clusters = []
        assigned = [False] * len(inds)
        for j, val in enumerate(mm):
            if assigned[j]:
                continue
            cluster = [j]
            assigned[j] = True
            for k in range(j + 1, len(mm)):
                if abs(mm[k] - val) <= tol:
                    cluster.append(k)
                    assigned[k] = True
            clusters.append(cluster)

        # --- формирование атомных типов (Fe, Fe_1, Fe_2 и т.д.) ---
        for ci, cluster in enumerate(clusters, start=1):
            site_inds = [inds[idx] for idx in cluster]
            avg_mag = float(np.mean([magmoms[idx] for idx in site_inds]))
            Z = Element(elem).Z
            if len(clusters) == 1:
                label = elem
            else:
                elem_counter[elem] += 1
                label = f"{elem}_{elem_counter[elem]}"
            types.append((label, Z, avg_mag, site_inds.copy()))
            for site_idx in site_inds:
                site2type[site_idx] = tindex
            tindex += 1

    return types, site2type, magmoms



def generate_sys_data(structure, magmoms=None):
    """
    Анализирует структуру, возвращает все данные, необходимые
    для генерации .sys и .pot файлов в формате SPR-KKR/xband.

    Возвращает словарь с ключами:
        bravais
        system_name
        prim_structure
        atom_types  [(label, Z, avg_mag, site_indices), ...]
        site2type   [int,...]  — сопоставление каждого сайта типу
        rws_dict    {elem: rws}
        prim_matrix (3x3)
        a_au        (float)
        a_ang       (float)
        spacegroup  (int)
        magmoms     reduced list of mags
    """
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    import numpy as np

    sga = SpacegroupAnalyzer(structure, symprec=1e-3)
    prim = sga.get_primitive_standard_structure()
    conv = sga.get_conventional_standard_structure()

    pearson = sga.get_pearson_symbol()[:2]
    br = [str(i) for i in Pearson.from_symbol(pearson)]
    brave = f'{br[1]:>13}        {br[2]:<12}{br[3]:<15}{br[4]:<7}{br[5]:<6}'

    # 'BRAVAIS            9        tetragonal  body-centered  4/mmm  D_4h'
    # 'BRAVAIS           13        cubic       face-centered  m3m    O_h '

    nat = Poscar(prim).natoms
    system_name = "".join([f'{el.symbol}{n}' if n != 1 else f'{el.symbol}' for el, n in zip(prim.composition.elements, nat)])

    # типы атомов
    atom_types, site2type, magmom = split_atom_types_by_magmom(prim, magmoms)

    # решётка
    a_ang = conv.lattice.a
    ANG_TO_AU = 1.889726125
    a_au = a_ang * ANG_TO_AU
    prim_matrix = prim.lattice.matrix / a_ang

    # радиусы Вигнера–Зейтца
    rws_dict = get_rws(prim)

    return {
        "bravais": brave,
        "system_name": system_name,
        "prim_structure": prim,
        "conv_structure": conv,
        "atom_types": atom_types,
        "site2type": site2type,
        "rws_dict": rws_dict,
        "prim_matrix": prim_matrix,
        "a_au": a_au,
        "a_ang": a_ang,
        "spacegroup": sga.get_space_group_number(),
        "magmoms": magmom
    }
