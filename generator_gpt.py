import os
from pymatgen.io.vasp import Poscar, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pathlib import Path
from pot_from_gpt import generate_pot_file
from sys_data_gpt import generate_sys_data
from sys_from_gpt import generate_sys_from_data
from Inp_files import generate_inps

# pos = Poscar.from_file(Path("for_spr/Ti4Fe8Cu4/119/FiM/POSCAR"))
# magmoms = [-0.534, 1.755, 2.137, -0.035]  # при необходимости

to_calc = []

wd = Path('SPR_KKR')
# alloys = ['Mn8Cr4Pt4', 'Ti4Fe8Cu4']
alloys = [i for i in os.listdir(wd) if i != 'tested']
for alloy in alloys:
    groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
    for group in groups:
        # orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
        # for order in orders:

        # path = Path(f'/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/{alloy}/{group}')
        path = Path(f'{wd}/{alloy}/{group}')
        print(str(path).replace('\\', '/'))
        # for i in os.listdir(path):
        #     if i not in ['magmoms', 'POSCAR']:
        #         os.remove(f'{path}/{i}')
        # continue
        # continue
        pos = Poscar.from_file(Path(f"{path}/CONTCAR"))
        sga = SpacegroupAnalyzer(pos.structure)
        # print(f'{alloy}_{group}: {sga.get_space_group_number()}')
        # continue

        if os.path.isfile(f'{path}/OUTCAR'):
            magmoms = [i['tot'] for i in Outcar(f'{path}/OUTCAR').magnetization]

        elif os.path.isfile(f'{path}/magmoms'):
            with open(f'{path}/magmoms') as f:
                magmoms = list(map(float, f.read().replace('\n', '').split()))
        else:
            print('magmoms are not defined')
            continue

        structure = pos.structure

        data = generate_sys_data(structure, magmoms)
        # print(data['magmoms'])

        pot_text = generate_pot_file(data)
        with open(f'{path}/{data["system_name"]}.pot', 'w') as f:
            f.write(pot_text)
        # print(pot_text)
        # print('#'*100)

        sys_text = generate_sys_from_data(data)
        with open(f'{path}/{data["system_name"]}.sys', 'w') as f:
            f.write(sys_text)
        # print(sys_text)
        print('#' * 100)

        scf, jxc = generate_inps(data["system_name"], data["magmoms"])
        print(scf)

        with open(f'{path}/SCF.inp', 'w') as f:
            f.write(scf)

        with open(f'{path}/JXC.inp', 'w') as f:
            f.write(jxc)
