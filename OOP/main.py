import os
from pathlib import Path
from Potential import generate_pot
from System import generate_sys


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


wd = Path(f'for_spr')
total = counter(wd)
current = 0
alloys = [i for i in os.listdir(f'{wd}') if os.path.isdir(f'{wd}/{i}')]
for alloy in alloys:
    groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
    for group in groups:
        orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
        # min_e = 0
        for order in orders:
            path = f'{wd}/{alloy}/{group}/{order}'
            current += 1
            sys = generate_sys(f'{path}/POSCAR')
            with open(f'{path}/system.sys', 'w') as f:
                f.write(sys)

            pot = generate_pot(f'{path}/POSCAR')
            with open(f'{path}/potential.pot', 'w') as f:
                f.write(pot)




            # TODO: CLEANER
            # files = [i for i in os.listdir(f'{wd}/{alloy}/{group}/{order}')]
            # for file in files:
            #     path = f'{wd}/{alloy}/{group}/{order}/{file}'
            #     if file not in ['POSCAR', 'magmoms']:
            #         try:
            #             shutil.rmtree(path)
            #         except:
            #             os.remove(path)