import os
import shutil

from ase import Atoms
from ase2sprkkr.potentials.potentials import Potential
# from ase2sprkkr.sprkkr.sysfile import write_sysfile, sysfile_content
from ase.spacegroup import crystal
import numpy as np
from pathlib import Path
from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms
from ase.spacegroup import get_spacegroup
from ase2sprkkr.physics.lattice_data import LatticeData
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core import Structure
from ase.atom import Atom
from sysfile_cont import sysfile_content


scf_inp = '''
###############################################################################
#  SPR-KKR input file    $name_SCF.inp 
###############################################################################
 
CONTROL  DATASET     = $name 
         ADSI        = SCF 
         POTFIL      = $name.pot 
         PRINT = 0    
 
MODE     SP-SREL 
 
TAU      BZINT= POINTS  NKTAB= 4000 
 
ENERGY   GRID={5}  NE={30} 
         ImE=0.0 Ry   EMIN=-0.2 Ry
 
SCF      NITER=200 MIX=0.05 VXC=PBE
         TOL=0.00001  MIXOP=0.20  ISTBRY=1 
         QIONSCL=0.80 
         NOSSITER 
'''

basis = np.array([[0, 0.5, 0.5],
                  [0.5, 0, 0.5],
                  [0.5, 0.5, 0]])


def create_spr_inputs(path: Path, out_path: Path):
    """
    Generate and create scf.inp, jxc.inp, *.sys, *.pot
    TODO select type of return_data
    :param path: str-path to folder with POSCAR and INCAR of target system
    :return: None or input_file strs
    """
    # pos = Poscar.from_file(r'poscars\POSCAR_L21')
    # TODO: если 8 атомов, поменять местами 2 и 4 атомы
    pos = Poscar.from_file(f'{path}')
    sites_dict = pos.as_dict()['structure']['sites']
    sites = []
    for line in sites_dict:
        sites.append(line['xyz'])
    sites = np.array(sites)
    matrix = pos.structure.lattice.matrix

    # print(pos.structure)
    # print(pos.structure.lattice.matrix)
    # print(pos.natoms)
    # print(pos.site_symbols)
    # print(sites)
    # print(matrix)
    # print(pos.structure.lattice.abc)
    # space = SpacegroupAnalyzer(pos.structure)
    # print(space.get_primitive_standard_structure())
    # print(matrix)
    # exit()
    positions = sites
    # print(positions)
    name = ''.join([f'{pos.site_symbols[i]}{pos.natoms[i]}' for i in range(len(pos.natoms))])
    atoms = Atoms(symbols=name,
                  cell=matrix,
                  positions=positions,
                  pbc=True)
    # atoms_2 = SPRKKRAtoms.promote_ase_atoms(atoms)
    # print(atoms_2.sites)
    # convert = np.array([[0.,  0.5, 0.5],
    #                     [0.5, 0.,  0.5],
    #                     [0.5, 0.5, 0. ]])
    # a_lat = np.dot(pos.structure.lattice.matrix, np.linalg.inv(basis))
    # print(a_lat)
    # for i in atoms.get_cell():
    #     print(i)

    # print(sysfile_content(atoms))

    # a, b, c = matrix[1, 0]*2, matrix[0, 1]*2, matrix[0, 2]*2
    # print(a, b, c)
    # with open(f'for_spr/sys_template_{spacegroup}', 'r') as f:
    #     sys_template = f.read()
    # sys_template = sys_template.replace('system_name', name)
    # if 'austenite' in str(path):
    #     volume = a * b * c
    #     for lat in ['A', 'B', 'C']:
    #         sys_template = sys_template.replace(f'{lat}_LATTICE', f'{volume**(1/3) * 1.8897259886}')
    #     sys_template = sys_template.replace('CA_LAT', f'1.0')
    # elif 'martensite' in str(path):
    #     sys_template = sys_template.replace('A_LATTICE', f'{a*1.8897259886}')
    #     sys_template = sys_template.replace('B_LATTICE', f'{b*1.8897259886}')
    #     sys_template = sys_template.replace('C_LATTICE', f'{c*1.8897259886}')
    #     sys_template = sys_template.replace('CA_LAT', f'{c/a}')


    #TODO calculate radius vigner zeits ne znayu kak na angliyskom
    # for i in range(1, 5):
    #     sys_template = sys_template.replace(f'WS_radius_{i}', '0.000000000000')

    # if len(pos.site_symbols) == 2:
    #     numbers = [Atom(pos.site_symbols[0]).number, Atom(pos.site_symbols[0]).number,
    #                Atom(pos.site_symbols[1]).number, Atom(pos.site_symbols[1]).number]
    #     symbols = [f'{pos.site_symbols[0]}_1', f'{pos.site_symbols[0]}_2',
    #                f'{pos.site_symbols[1]}_1', f'{pos.site_symbols[1]}_2']
    #     for i in range(1, 5):
    #         sys_template = sys_template.replace(f'atom_num_{i}', f'{numbers[i-1]}')
    #         sys_template = sys_template.replace(f'atom_name_{i}', f'{symbols[i-1]}')
    #
    # elif len(pos.site_symbols) == 3:
    #
    #     symbols = []
    #     for i, s in zip(pos.natoms, pos.site_symbols):
    #         if i == 2:
    #             symbols.append(f'{s}_1')
    #             symbols.append(f'{s}_2')
    #         else:
    #             symbols.append(s)
    #
    #     elements = [s.split('_')[0] if '_' in s else s for s in symbols]
    #
    #     for i, s in enumerate(symbols):
    #         sys_template = sys_template.replace(f'atom_num_{i + 1}', f'{Atom(elements[i]).number}')
    #         sys_template = sys_template.replace(f'atom_name_{i + 1}', f'{s}')
    #
    # elif len(pos.site_symbols) == 4:
    #     for i, s in enumerate(pos.site_symbols):
    #         sys_template = sys_template.replace(f'atom_num_{i + 1}', f'{Atom(s).number}')
    #         sys_template = sys_template.replace(f'atom_name_{i + 1}', f'{s}')
    # else:
    #     print('INCORRECT SET OF ELEMENTS')
    #     exit()


    # print(sysfile_content(atoms))
    # print(sys_template)
    with open(f'{out_path}/{name}.sys', 'w') as f:
        f.write(sysfile_content(atoms))

    pot = Potential.from_atoms(atoms)
    pot.save_to_file(f'{out_path}/{name}.pot')

    with open(f'{out_path}/SCF.inp', 'w') as f:
        global scf_inp
        scf_inp = scf_inp.replace('$name', name)
        f.write(scf_inp)
    #
    # print(f'{out_path}\\{name}.txt')


for i in ['216']:
    for j in os.listdir(f'for_spr/{i}'):
        # if j != 'Zr8Nb8':
        #     continue
        for k in os.listdir(f'for_spr/{i}/{j}'):
            # if k != 'martensite':
            #     continue
            # pos = Poscar.from_file(Path(f'for_spr/{i}/{j}/{k}/0/POSCAR'))
            # print(pos)
            # print('#'*100)
            create_spr_inputs(Path(f'for_spr/{i}/{j}/{k}/0'), int(i), Path(f'for_spr/{i}/{j}/{k}/0'))
# create_spr_inputs(Path(f'poscars/POSCAR_XA_PRIM'), 216, Path(f'test/Ni2MnGa_aust_auto'))
# create_spr_inputs(Path(f'poscars/POSCAR_XA_prim_mart'), 225, Path(f'test/Ni2MnGa_mart_auto'))


