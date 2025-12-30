import os
from pathlib import Path
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element


def get_structures(poscar_path: Path):
    '''
    return raw structure, converted to prim, conerted to conventional
    :param poscar_path:
    :return: raw structure, converted to prim, conerted to conventional
    '''
    structure = Structure.from_file(poscar_path)
    structure_primitive = SpacegroupAnalyzer(structure).get_primitive_standard_structure()
    structure_conventional = SpacegroupAnalyzer(structure).get_conventional_standard_structure()

    return structure, structure_primitive, structure_conventional


def get_magmoms_dict(path: Path):
    struc = Structure.from_file(f'{path}/POSCAR')

    magmoms = {s.symbol: [] for s in struc.elements}
    # print(magmoms)
    if os.path.isfile(f'{path}/OUTCAR'):
        out = Outcar(f'{path}/OUTCAR')
        mags = out.magnetization
        # magmoms = [str(i['tot']) for i in out.magnetization]
        for i, site in enumerate(struc.sites):
            magmoms[site.specie.symbol].append(str(mags[i]['tot']))
    else:
        with open(f'{path}/magmoms') as f:
            temp = f.readlines()[0].split()

            for i, site in enumerate(struc.sites):
                magmoms[site.specie.symbol].append(temp[i])

    for key, val in magmoms.items():
        magmoms[key] = val[0::2]

    return magmoms


def get_magmoms_list(path: Path):
    with open(f'{path}/magmoms') as f:
        temp = f.readlines()[0].split()
        return temp


def get_magmoms_str(path: Path):
    if os.path.isfile(f'{path}/OUTCAR'):
        out = Outcar(f'{path}/OUTCAR')
        magmoms = [str(i['tot']) for i in out.magnetization]
    else:
        with open(f'{path}/magmoms') as f:
            magmoms = f.readlines()[0].split()

    if len(magmoms) > 4:
        magmoms = magmoms[0::2]
    temp = ','.join(magmoms)
    return temp


def get_labels_for_sites(path: Path):
    '''
    get dict with labels for sites with unique mag moment {label: element}
    :param path:
    :return dict of {label: element}:
    '''
    labels = []
    elements = []
    unique_indexes = []
    index = 0
    mags = get_magmoms_dict(path)
    # if len(mags.values()) == len(set(mags.values())):
    # print(mags)
    for index, (key, val) in enumerate(get_magmoms_dict(path).items()):
        if len(val) == 1:
            labels.append(f'{key}')
            elements.append(Element(key))
            unique_indexes.append(index + 1)
        elif len(val) == 2 and len(set(val)) == 1:
            for i, v in enumerate(val):
                labels.append(f'{key}')
                elements.append(Element(f'{key}'))
                unique_indexes.append(index + 1)
        else:
            for i, v in enumerate(val):
                labels.append(f'{key}_{i+1}')
                elements.append(Element(f'{key}'))
                unique_indexes.append(index + 1)
    # return {label: el for label, el in zip(labels, elements)}
    return labels, elements, unique_indexes


def get_type_indexes(path: Path, magmom: list):
    # repeat = 0
    # indexes = []
    # for i in range(len(magmom)):
    #     if i == 0:
    #         pass
    #
    #     if magmom[i] == magmom[i - 1]:
    #         repeat += 1
    #     indexes.append(i - repeat + 1)
    if os.path.isfile(f'{path}/OUTCAR'):
        out = Outcar(f'{path}/OUTCAR')
        magmoms = [str(i['tot']) for i in out.magnetization]
        magmoms = [magmoms[i] for i in range(0, len(magmoms), len(magmoms) % 4)]

    indexes = []
    repeat = 0
    for i in range(len(magmom)):
        if i == 0:
            pass

        if magmom[i] != magmom[i - 1]:
            repeat += 1

        indexes.append(repeat)
    return indexes


def get_types(path: Path):
    mag_d = get_magmoms_dict(path)
    types = []
    for key, item in mag_d.items():
        if len(set(item)) == 1 and len(item) > 1:
            types += [key] * len(item)
        elif len(item) == 1:
            types += [key]
        else:
            types += [f'{key}_{i+1}' for i, _ in enumerate(item)]
    return types


def get_sites_info(path: Path):
    sites = []
    keys = ['site_index', 'type_index', 'magmom', 'element', 'type']

    norm, prim, conv = get_structures(Path(f'{path}/POSCAR'))
    print(prim)
    # mag_l = get_magmoms_list(path)
    mag_l = list(dict.fromkeys([float(f'{float(i):.2f}') for i in get_magmoms_list(path)]))
    # mag_l = list(dict.fromkeys(mag_l))
    # print(mag_l)

    mag_d = get_magmoms_dict(path)
    # print(mag_d)
    a = list(mag_d.values())
    mag_l2 = [element for item in a for element in item]
    # mag_l2 = a[0] + a[1] + a[2]

    types = get_types(path)
    type_indexes = get_type_indexes(path, mag_l2)
    # return
    # types = [f'{i.specie.symbol}' for i in prim if len(mag_d[i.specie.symbol]) == 1 else f'{i.specie.symbol}_{}']
    print()

    for i, site in enumerate(prim.sites):
        # print(site)
        sites.append(dict.fromkeys(keys))
        sites[i]['site_index'] = i + 1
        sites[i]['type_index'] = type_indexes[i]
        sites[i]['magmom'] = mag_l2[i]
        sites[i]['element'] = site.specie.symbol
        sites[i]['type'] = types[i]
        # print(sites[i])
    return sites



if __name__ == '__main__':
    # _, s2, s3 = map(Poscar, get_structures(Path('for_spr/Ti4Fe8Cu4/216/FiM/POSCAR')))
    # # print(s1)
    # print(s2)
    # print(s3)
    # print(get_labels_for_sites(Path('for_spr/Ti4Fe8Cu4/216/FiM')))
    get_sites_info(Path('for_spr/Ti4Fe8Cu4/216/FiM'))
    get_sites_info(Path('for_spr/Mn8Cr4Pt4/139/FiM'))

    exit()
    wd = Path('for_spr')
    alloys = ['Mn8Cr4Pt4', 'Ti4Fe8Cu4']
    for alloy in alloys:
        groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
        for group in groups:
            orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
            for order in orders:
                path = f'{wd}/{alloy}/{group}/{order}'
                # print(get_labels_for_sites(path))