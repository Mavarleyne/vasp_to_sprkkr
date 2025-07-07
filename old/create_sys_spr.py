import os, sys
from ase import Atoms
from ase2sprkkr.potentials.potentials import Potential
import numpy as np
from pathlib import Path
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar
from sysfile_cont import sysfile_content
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


scf_inp = '''
###############################################################################
#  SPR-KKR input file    NAME_SCF.inp 
###############################################################################

CONTROL  DATASET     = NAME 
         ADSI        = SCF 
         POTFIL      = NAME.pot 
         PRINT = 0    

MODE     SP-SREL 

TAU      BZINT= POINTS  NKTAB= 4000 

ENERGY   GRID={5}  NE={30} 
         ImE=0.0 Ry   EMIN=-0.2 Ry

SCF      NITER=200 MIX=0.05 VXC=PBE
         TOL=0.00001  MIXOP=0.20  ISTBRY=1 
         QIONSCL=0.80 
         NOSSITER
         MSPIN={MAGMOMS} 
'''

basis = np.array([[0, 0.5, 0.5],
                  [0.5, 0, 0.5],
                  [0.5, 0.5, 0]])


def counter(base_dir: str):
    count = 0
    for alloy in os.listdir(base_dir):
        alloy_path = os.path.join(base_dir, alloy)
        if not os.path.isdir(alloy_path):
            continue

        for group in os.listdir(alloy_path):
            group_path = os.path.join(alloy_path, group)
            if not os.path.isdir(group_path):
                continue

            for order in os.listdir(group_path):
                order_path = os.path.join(group_path, order)
                if not os.path.isdir(order_path):
                    continue
                count += 1
    return count


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


def get_sites_from_pos(pos: Poscar):
    sites_dict = pos.as_dict()['structure']['sites']
    sites = []
    for line in sites_dict:
        sites.append(line['xyz'])
    return np.array(sites)


def create_spr_inputs(path: Path, out_path: Path):
    """
    Generate and create scf.inp, jxc.inp, *.sys, *.pot
    :param path: Pathlib-obj to folder with POSCAR of target system
    :param out_path: Pathlib-obj to save-folder
    :return: None
    """
    if os.path.isfile(f'{path}\\CONTCAR') and os.path.getsize(f'{path}\\CONTCAR') > 10:
        pos = Poscar.from_file(f'{path}\\CONTCAR')
    else:
        pos = Poscar.from_file(f'{path}\\POSCAR')
    space = SpacegroupAnalyzer(pos.structure)
    print(space.get_space_group_number())
    print(Poscar(space.get_conventional_standard_structure()))
    space = SpacegroupAnalyzer(space.get_conventional_standard_structure()).get_primitive_standard_structure()
    pos_2 = Poscar(space)
    print(SpacegroupAnalyzer(pos_2.structure).get_space_group_number())
    pos_2.write_file(f'{path}\\POSCAR_prim')
    # print(space.get_space_group_number())
    # print(space.cart_coords)
    from System import struc
    space = Poscar.from_str(struc).structure
    sites = np.array(space.frac_coords)
    matrix = space.lattice.matrix
    # print(space.composition.to_pretty_string())
    # print(space.num_sites)
    # exit()
    # name = ''.join([f'{pos.site_symbols[i]}{pos.natoms[i]}' for i in range(len(pos.natoms))])
    name = space.composition.to_pretty_string()

    with open(f'{path}\\magmoms') as f:
        atoms_mag = list(map(float, f.read().split()))

    if pos.structure.num_sites == 8:
        mspin = [atoms_mag[i] for i in range(0, len(atoms_mag), 2)]
    elif pos.structure.num_sites == 4:
        mspin = [atoms_mag[i] for i in range(0, len(atoms_mag), 1)]
    else:
        mspin = None
        print('num of sites not 4 and 8')
    # print(mspin)

    print(out_path)
    print(space)
    atoms = Atoms(symbols=name,
                  cell=matrix,
                  positions=sites,
                  # scaled_positions=sites,
                  # magmoms=mspin,
                  pbc=True)
    # print('#'*100)
    # print(atoms.cell.get_bravais_lattice().pearson_symbol)
    # exit()
    inp = sysfile_content(atoms)
    # with open(f'{out_path}/{name}.sys', 'w') as f:
        # f.write(inp)
    print(inp)

    pot = Potential.from_atoms(atoms)
    # pot.save_to_file(f'{out_path}\\{name}.pot')


    with open(Path(f'{out_path}/SCF.inp'), 'w') as f:
        global scf_inp

        temp = scf_inp.replace('NAME', name).replace('MAGMOMS', ','.join([str(i) for i in mspin]))
        f.write(temp)



if __name__ == '__main__':
    wd = Path(f'for_spr')
    # TODO: fix problem with Zr4Mn4Co4W4 mart
    total = counter(wd)
    current = 0
    alloys = [i for i in os.listdir(f'{wd}') if os.path.isdir(f'{wd}/{i}')]
    for alloy in alloys:
        groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
        for group in groups:
            orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
            min_e = 0
            for order in orders:
                path = f'{wd}/{alloy}/{group}/{order}'
                current += 1
                # print_progress(current, total, 80)

                # pos = Poscar.from_file(Path(f'for_spr/{i}/{j}/{k}/0/POSCAR'))
                # print(pos)
                # print('#'*100)
                create_spr_inputs(Path(f'{wd}/{alloy}/{group}/{order}'), Path(f'{wd}/{alloy}/{group}/{order}'))
                # exit()



