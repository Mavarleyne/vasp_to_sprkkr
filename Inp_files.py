import os.path
import shutil
from pathlib import Path
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from utils import get_magmoms_str
import Potential
import System_new
from Potential import generate_pot

scf_inp_template = '''###############################################################################
#  SPR-KKR input file    *SCF.inp 
###############################################################################

CONTROL  DATASET     = {system} 
         ADSI        = SCF 
         POTFIL      = {system}.pot 
         PRINT = 0    

MODE     SP-SREL 

TAU      BZINT= POINTS  NKTAB= {nktab} 

ENERGY   GRID={5}  NE={30} 
         ImE=0.0 Ry   EMIN={emin} Ry

SCF      NITER={niter} MIX={mix} VXC=PBE
         TOL=0.00001  MIXOP={mixop}  ISTBRY=1 
         QIONSCL=0.80 
         NOSSITER
         MSPIN={{magmoms}} 
'''

replaces_scf = {'{system}': None,
                '{nktab}': '4000',
                '{emin}': '-0.2',
                '{magmoms}': None,
                '{niter}': '200',
                '{mix}': '0.05',
                '{mixop}': '0.20',
                '{tol}': '0.00001'}

jxc_inp_template = '''###############################################################################
#  SPR-KKR input file    *JXC.inp 
#  created by xband on Sat 14 Jun 22:38:43 BST 2025
###############################################################################
 
CONTROL  DATASET     = {system} 
         ADSI        = JXC 
         POTFIL      = {system}.pot 
         PRINT = 0    
 
MODE     SP-SREL 
 
TAU      BZINT= POINTS  NKTAB= {nktab} 
 
ENERGY   GRID={5}  NE={30} 
         EMIN={emin} Ry
 
TASK     JXC   CLURAD={cluster_radius}
'''

replaces_jxc = {'{system}': None,
                '{nktab}': '2000',
                '{emin}': '-0.2',
                '{cluster_radius}': '2.5'}

templates = {'scf': scf_inp_template,
             'jxc': jxc_inp_template}

replaces = {'scf': replaces_scf,
            'jxc': replaces_jxc}


def create_inp_file(type: str, tags_to_replace: dict, output_path: Path=None):
    if type not in ['scf', 'jxc']:
        raise KeyError('Incorrect type of input file, try "scf" or "jxc"')
    inp = templates[type]
    for key, val in tags_to_replace.items():
        inp = inp.replace(key, val)
    # with open(f'{output_path}/{type.upper()}.inp', 'w') as f:
    #     f.write(inp)
    return inp


def get_tags_scf(tags_to_replace: dict, name, magmoms):
    tags_to_replace['{system}'] = name
    tags_to_replace['{magmoms}'] = magmoms
    mags = magmoms.split()
    if len(mags) == 4:
        tags_to_replace['{nktab}'] = '4000'
    elif len(mags) == 8:
        tags_to_replace['{nktab}'] = '2000'
    elif len(mags) == 16:
        tags_to_replace['{nktab}'] = '1000'
    return tags_to_replace


def get_tags_jxc(tags_to_replace: dict, name):
    tags_to_replace['{system}'] = name
    return tags_to_replace


def generate_inps(name: str, magmoms: list):
    mags = ' '.join([str(i) for i in magmoms])
    scf = get_tags_scf(replaces_scf, name, mags)
    jxc = get_tags_jxc(replaces_jxc, name)
    # create_inp_file('scf', scf)
    # create_inp_file('jxc', jxc)
    return create_inp_file('scf', scf), create_inp_file('jxc', jxc)


if __name__ == '__main__':
    wd = Path('for_spr')
    # alloys = ['Mn8Cr4Pt4', 'Ti4Fe8Cu4']
    alloys = [i for i in os.listdir(wd) if i != 'tested']
    for alloy in alloys:
        groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
        for group in groups:
            orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
            for order in orders:

                path = f'{wd}/{alloy}/{group}/{order}'
                # if os.path.exists(f'{path}/vampire'):
                #     shutil.rmtree(f'{path}/vampire')
                # for i in os.listdir(path):
                #     if i not in ['magmoms', 'POSCAR']:
                #         os.remove(f'{path}/{i}')
                # continue
                # pot, name = generate_pot(path)
                # Potential.write_pot(pot, f'{path}/{name}.pot')
                #
                # sys = System_new.generate_sys(f'{path}')
                # System_new.write_sys(sys, f'{path}/{name}.sys')

                # main(name, path)


    # main()
    # print(get_magmoms_dict(Path('for_spr/Ta4Fe8Ru4/216/FM')))
    # mags = get_magmoms_list(Path('for_spr/Ta4Fe8Ru4/216/FM'))
    # print(get_magmoms_str(Path('for_spr/Ta4Fe8Ru4/216/FM')))
    # struc =
