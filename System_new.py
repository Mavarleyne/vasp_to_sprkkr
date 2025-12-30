from collections import namedtuple, Counter
from pathlib import Path
from datetime import datetime
from pymatgen.core import Structure, Element
from pymatgen.io.vasp import Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SgA
from pymatgen.io.vasp.inputs import Poscar
import numpy as np
from spglib import spglib

import utils
from Wigner_Seitz_radius import get_rws
from brave_from_pearson import Pearson, international_numbers_to_AP

struc = '''Mn4 Cr2 Pt2                             
   1.00000000000000     
     0.0000000000000004    2.7099026210485557    3.6053512551112084
     5.4198052420971115    0.0000000000000000    0.0000000000000004
     0.0000000000000004    2.7099026210485557   -3.6053512551112079
  Mn/             Cr/             Pt/           
               4               2               2
Direct
  0.5000000000000000  0.2500000000000000  0.0000000000000000
  0.0000000000000000  0.7500000000000000  0.5000000000000000
  0.5000000000000000  0.7500000000000000  0.0000000000000000
  0.0000000000000000  0.2500000000000000  0.5000000000000000
  0.0000000000000000  0.5000000000000000  0.0000000000000000
  0.5000000000000000  0.0000000000000000  0.5000000000000000
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.5000000000000000  0.5000000000000000  0.5000000000000000

  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00
'''

s = Poscar.from_str(struc)
# print(7.24 * 7.24 * 13.6)
# print(s.structure.volume)


def get_magmoms_from_outcar(outcar_path: Path):
    outcar = Outcar(outcar_path)
    return [i['tot'] for i in outcar.magnetization]


def get_magnetic_primitive(structure: Structure, magmoms: list):
    magnetic_structure = Structure(structure.lattice.matrix,
                                   structure.species,
                                   structure.frac_coords,
                                   site_properties={"magmom": magmoms})
    # Использование spglib для магнитной симметрии
    cell = (magnetic_structure.lattice.matrix,
            magnetic_structure.frac_coords,
            [magnetic_structure.species.index(sp) for sp in magnetic_structure.species],
            magmoms)

    dataset = spglib.get_magnetic_symmetry_dataset(cell, symprec=1e-5)
    print(f"Magnetic space group: {dataset['uni_number']}")

    # Получение примитивной ячейки
    primitive_cell = spglib.find_primitive(cell, symprec=1e-5)
    if primitive_cell:
        print(primitive_cell)
        prim_lattice, prim_positions, prim_numbers, prim_magmoms = primitive_cell
        prim_structure = Structure(
            lattice=prim_lattice,
            species=[structure.species[n] for n in prim_numbers],
            coords=prim_positions,
            coords_are_cartesian=False,
            site_properties={"magmom": prim_magmoms}
        )
        print(f"Primitive magnetic structure: {prim_structure}")
        # prim_structure.to(filename="POSCAR_primitive", fmt="poscar")


def get_sites_coords(path: str):
    sample = '''number of sites NQ
  4
 IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ
  1   1    0.000000000000    0.500000000000    0.470380000000      2.724114898250   3    1   1
  2   1    0.500000000000    0.000000000000    0.470380000000      2.724114898250   3    1   1
  3   2    0.000000000000    0.000000000000    0.940760000000      2.708975319341   3    1   2
  4   3    0.000000000000    0.000000000000    0.000000000000      2.923957339840   3    1   3
'''
    return


def get_sites_equivalence(path):
    sample = '''number of sites classes NCL
  3
ICL WYCK NQCL IQECL (equivalent sites)
  1   -    2  1  2
  2   -    1  3
  3   -    1  4
'''
    return


def get_sites_occupation(path: str):
    sample = '''number of atom types NT
  3
 IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)
  1  25     Mn       2 1.000  1  2
  2  24     Cr       1 1.000  3
  3  78     Pt       1 1.000  4
'''
    return


def generate_sys(path: Path):
    # --- Настройки ---
    output_filename = "structure.sys"
    structure = Structure.from_file(f'{path}/POSCAR')  # или .vasp, POSCAR и т.д.
    # structure = Structure.from_str(struc, 'Poscar')
    sga = SgA(structure)
    pearson = sga.get_pearson_symbol()[:2]
    rwss = get_rws(structure)  # Радиус Вигнера-Сейца по умолчанию (в а.е.)

    structure = sga.get_conventional_standard_structure()
    sga = SgA(structure)

    # --- Преобразование ---
    bohr = 1.889726125  # 1 Å = 1.8897 Bohr

    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles

    lattice_param_a = a * bohr
    b_over_a = b / a
    c_over_a = c / a

    structure_prim = SgA(structure).get_primitive_standard_structure()
    primitive_vectors = structure_prim.lattice.matrix / structure.lattice.\
        a

    # Получаем список уникальных элементов
    unique_elements = list({site.specie for site in structure_prim})

    # Подсчитываем количество атомов каждого элемента
    # element_counts = Counter(site.specie for site in structure_prim)

    # Сортируем по убыванию количества атомов
    # unique_elements.sort(key=lambda e: element_counts[e], reverse=True)
    # print(unique_elements)
    # labels = utils.get_labels_for_sites(path)
    # elements = list(map(Element, list(labels.values())))
    # labels = list(labels.keys())
    # if len(labels) == len(unique_elements):
    #     element_to_it = {el.symbol: i + 1 for i, el in enumerate(unique_elements)}
    # else:
    #     element_to_it = {el: i + 1 for i, el in enumerate(labels)}
    labels, elements, unique_indexes = utils.get_labels_for_sites(path)
    num_sites = {i: labels.count(i) for i in labels}

    s_inf = utils.get_sites_info(path)

    # --- Запись .sys файла ---
    with open(output_filename, "w") as f:
        print('#' * 100)
        filestring = ''
        filestring += f"system data-file created by pymatgen2sprkkr\n"
        filestring += f"<generated-system>\n"
        filestring += f"xband-version\n5.0\n"
        filestring += f"dimension\n3D\n"
        filestring += f"{' '.join([str(i) for i in Pearson.from_symbol(pearson)][1:])}\n"
        filestring += f"space group number (ITXC and AP)\n" \
                      f"{sga.get_space_group_number()} {international_numbers_to_AP[sga.get_space_group_number()]}\n"
        filestring += f"structure type\nUNKNOWN\n"

        filestring += f"lattice parameter A  [a.u.]\n" \
                      f"{lattice_param_a:.12f}\n"
        filestring += f"ratio of lattice parameters  b/a  c/a\n" \
                      f"{b_over_a:.12f} {c_over_a:.12f}\n"
        filestring += f"lattice parameters  a b c  [a.u.]\n" \
                      f"{(a * bohr):.12f} {(b * bohr):.12f} {(c * bohr):.12f}\n"
        filestring += f"lattice angles  alpha beta gamma  [deg]\n" \
                      f"{alpha:.12f} {beta:.12f} {gamma:.12f}\n"

        filestring += f"primitive vectors     (cart. coord.) [A]\n"
        for vec in primitive_vectors:
            filestring += f"  {vec[0]:.12f}  {vec[1]:.12f}  {vec[2]:.12f}\n"

        filestring += f"number of sites NQ\n  {len(structure_prim)}\n"
        filestring += f" IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"
        for i, site in enumerate(structure_prim.sites):
            pos = site.coords / structure.lattice.a
            element = site.specie.symbol
            ito = unique_indexes[i]
            rws = rwss[site.specie.symbol]
            filestring += f"{s_inf[i]['site_index']:3d} {s_inf[i]['type_index']:3d} {pos[0]:16.12f} {pos[1]:16.12f} {pos[2]:16.12f} {rws:20.12f}   3    1  {s_inf[i]['type_index']:2d}\n"

        filestring += f"number of sites classes NCL \n  {s_inf[-1]['type_index']}\n"
        filestring += f"ICL WYCK NQCL IQECL (equivalent sites)\n"
        for i in range(len(labels)):
            filestring += f"  {i + 1}   -    1  {i + 1}\n"

        filestring += f"number of atom types NT\n  {s_inf[-1]['type_index']}\n"
        filestring += f" IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"

        ind = 1
        cont = False
        for i, site in enumerate(s_inf):
            if cont:
                cont = False
                continue

            IT = ind
            ZT = Element(site['element']).Z
            TXTT = site['type']
            NAT = 1
            IQAT = site['site_index']

            if i != len(s_inf)-1 and site['magmom'] == s_inf[i + 1]['magmom']:
                IQAT = f"{site['site_index']} {s_inf[i + 1]['site_index']}"
                NAT = 2
                cont = True

            filestring += f' {IT:2}  {ZT:2}  {TXTT:4}  {NAT:3}  1     {IQAT}\n'
            ind += 1

        # print(filestring)
        return filestring


def write_sys(filestring: str, output_path: Path):
    with open(output_path, 'w') as f:
        f.write(filestring)


if __name__ == '__main__':
    path = Path('for_spr/Ta4Fe8Ru4/216/FM')
    struc = Structure.from_file(f'{path}/POSCAR')
    with open(f'{path}/magmoms', 'r') as f:
        magmoms = list(map(float, f.read().split()))
    get_magnetic_primitive(struc, magmoms)
    pass

