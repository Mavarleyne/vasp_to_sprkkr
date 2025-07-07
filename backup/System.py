from collections import namedtuple, Counter
from pathlib import Path
from datetime import datetime
from pymatgen.core import Structure, Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SgA
from pymatgen.io.vasp.inputs import Poscar
import numpy as np
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
print(7.24*7.24*13.6)
print(s.structure.volume)
def generate_sys(poscar: Path):
    # --- Настройки ---
    output_filename = "structure.sys"
    structure = Structure.from_file(poscar)  # или .vasp, POSCAR и т.д.
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
    primitive_vectros = structure_prim.lattice.matrix / structure.lattice.a

    # Получаем список уникальных элементов
    unique_elements = list({site.specie for site in structure_prim})

    # Подсчитываем количество атомов каждого элемента
    element_counts = Counter(site.specie for site in structure_prim)

    # Сортируем по убыванию количества атомов
    unique_elements.sort(key=lambda e: element_counts[e], reverse=True)
    # print(unique_elements)
    element_to_it = {el: i + 1 for i, el in enumerate(unique_elements)}

    # --- Запись .sys файла ---
    with open(output_filename, "w") as f:
        print('#' * 100)
        filestring = ''
        filestring += "system data-file created by pymatgen2sprkkr\n"
        filestring += "<generated-system>\n"
        filestring += "xband-version\n5.0\n"
        filestring += "dimension\n3D\n"
        filestring += f"{' '.join([str(i) for i in Pearson.from_symbol(pearson)][1:])}\n"
        filestring += f"space group number (ITXC and AP)\n{sga.get_space_group_number()} {international_numbers_to_AP[sga.get_space_group_number()]}\n"
        filestring += "structure type\nUNKNOWN\n"

        filestring += "lattice parameter A  [a.u.]\n%.12f\n" % lattice_param_a
        filestring += "ratio of lattice parameters  b/a  c/a\n%.12f %.12f\n" % (b_over_a, c_over_a)
        filestring += "lattice parameters  a b c  [a.u.]\n%.12f %.12f %.12f\n" % (a * bohr, b * bohr, c * bohr)
        filestring += "lattice angles  alpha beta gamma  [deg]\n%.12f %.12f %.12f\n" % (alpha, beta, gamma)

        filestring += "primitive vectors     (cart. coord.) [A]\n"
        for vec in primitive_vectros:
            filestring += "  %.12f  %.12f  %.12f\n" % tuple(vec)

        filestring += "number of sites NQ\n  %d\n" % len(structure_prim)
        filestring += " IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"
        for i, site in enumerate(structure_prim.sites):
            pos = site.coords / structure.lattice.a
            # print(site.coords)
            # print(structure.lattice.abc)
            # print(site.frac_coords)
            element = site.specie
            ito = element_to_it[element]
            rws = rwss[site.specie.symbol]
            filestring += "%3d %3d %16.12f %16.12f %16.12f %20.12f   3    1  %2d\n" % (i + 1, i + 1, *pos, rws, ito)

        filestring += "number of sites classes NCL \n  %d\n" % len(structure_prim)
        filestring += "ICL WYCK NQCL IQECL (equivalent sites)\n"
        for i in range(len(structure_prim)):
            filestring += "  %d   -    1  %d\n" % (i + 1, i + 1)

        filestring += "number of atom types NT\n  %d\n" % len(unique_elements)
        filestring += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
        for el, it in element_to_it.items():
            z = el.Z
            label = el.symbol
            iq_sites = [i + 1 for i, site in enumerate(structure_prim) if site.specie == el]
            filestring += " %2d %3d  %4s   %2d  1.000  %s\n" % (
                it, z, label, len(iq_sites), " ".join(map(str, iq_sites)))
        print(filestring)
        return filestring


if __name__ == '__main__':
    # generate_sys(Path("for_spr/Co8Ru4Ir4/119/FM/POSCAR"))
    pass
