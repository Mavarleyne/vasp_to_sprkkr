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


from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Пример структуры (замените на вашу)
struct = Structure.from_file(Path('SPR_KKR/Al/Tsharp/CONTCAR'))
# print(struct)
prim = SpacegroupAnalyzer(struct).get_primitive_standard_structure()

# Создайте анализатор
analyzer = SpacegroupAnalyzer(prim, symprec=0.01)  # symprec для стабильности

# Получите симметризованную структуру
# sym_struct = analyzer.get_symmetrized_structure()
sym_struct = analyzer.get_symmetrized_structure()


# Теперь доступны Wyckoff-данные:
print(prim)
print(sym_struct)
print("Буквы Вайкоффа (symbols):", sym_struct.wyckoff_symbols)  # Например: ['4a', '4b']
print("Метки с мультипликативностью (labels):\n", sym_struct.frac_coords)
print("Эквивалентные сайты:", sym_struct.equivalent_indices)

# Если нужно полный SpaceGroup объект:
# sg = analyzer.get_symmetry_dataset()
# print("Пространственная группа:", sg.symbol)  # Например: 'Fm-3m'
exit()
sites = {'IQ': [],
         'QX': [],
         'Qy': [],
         'QZ': []}

occupation = {'IQ': [], # Индекс атома
              'IREFQ': [], # Какому референс-потенциалу принадлежит
              'IMQ': [], # Магнитная подрешетка
              'NOQ': [], # Число состояний в данной позиции
              'ITOQ': [], # Тип атома (ссылается на блок TYPES)
              'CONC': []} # Концентрация (1.0 — чистый элемент)

reference_system = {'IREF': [],
                    'VREF': [], # Усредненный потенциал (для старта)
                    'RMTREF': []} # Радиус сферы МТ (может быть 0 на старте)

mag_direction = {'IQ': [],
                 'QMTET': [], # Угол θ (в радианах или градусах — зависит от настроек)
                 'QMPHI': []} # Угол φ (азимут)

mesh_info = {'IM': [], # Индекс сетки (по атомам)
             'R(1)': [], # 	Начальный радиус сетки
             'DX': [], # Шаг сетки
             'JRMT': [], # Индекс точки, соответствующей RMT
             'RMT': [], # Радиус МТ-сферы
             'JRWS': [], # Индекс для границы Wigner-Seitz
             'RWS': []} # Радиус WS-сферы

types = {'IT': [], # Номер типа атома
         'TXTT': [], # Химическое обозначение
         'ZT': [], # Заряд ядра
         'NCORT': [], # Кол-во электронов в core
         'NVALT': [], # Кол-во валентных
         'NSEMCORSHLT': []} # Полусердечные состояния (0 — игнорируются)

path = Path('for_spr/Ti4Fe8Cu4/119/FiM')

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
primitive_vectors = structure_prim.lattice.matrix / structure.lattice.a

print(SgA(structure).get_space_group_number())
print(SgA(structure_prim).get_space_group_number())
# print(primitive_vectors)


for i, site in enumerate(structure_prim.sites):
    pos = site.coords / structure.lattice.a
    print(site.specie.symbol)
    print(pos)

exit()

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
