import os
from pathlib import Path

jxc_qsub = '''#!/bin/bash
#PBS -d .
#PBS -l nodes=4:ppn=20
#PBS -N Jxc_100
#PBS -j oe
#PBS -l walltime=240:00:00


export HOME_DIR=`pwd`
echo $HOME_DIR

COMMANDS


wait
'''

def check_convergence(path: str):
    try:
        with open(f'{path}*SCF.out') as f:
            if 'SCF - cycle converged !!!!!!!!!' in f.read():
                return 'success'
            else:
                return 'not converge'
    except:
        return 'Excepption with opening *SCF.out'


def prepare_jxc_potential_files(directory: Path):
    """
    В указанной директории:
    - переименовывает файл *.pot → *.pot_old, если *.pot_old не существует
    - переименовывает файл *.pot_new → *.pot
    """
    path = directory

    # Переименование *.pot → *.pot_old (если нет *.pot_old)
    for pot_file in path.glob("*.pot"):
        pot_old = pot_file.with_name(pot_file.stem + ".pot_old")
        pot_new = pot_file.with_name(pot_file.stem + ".pot_new")
        if not pot_old.exists() and pot_new.exists():
            pot_file.rename(pot_old)
            pot_new.rename(pot_file)
            print(f"🔄 {pot_file.name} → {pot_old.name}")
            print(f"✅ {pot_new.name} → {pot_file.name}")
        else:
            print(f"⚠️ Пропущено: {pot_old.name} уже существует, {pot_new.name} не существует")

    # Переименование *.pot_new → *.pot
    # for new_file in path.glob("*.pot_new"):
    #     pot_target = new_file.with_name(new_file.stem.replace(".pot_new", "") + ".pot")
    #     new_file.rename(pot_target)
    #     print(f"✅ {new_file.name} → {pot_target.name}")


def create_jxc_inp(path: str):
    for file in os.listdir(path):
        if file[-4:] == '.pot':
            name = file[:-4]
            break
    print(name)
    prepare_jxc_potential_files(f'{path}')
    with open(f'{path}/*JXC.inp', 'w') as f:
        inp = jxc_inp.replace('#NAME#', name)
        current_date = datetime.datetime.now().strftime('%d-%m-%Y')
        inp = inp.replace('#DATE#', current_date)
        f.write(inp)

command = '(cd $HOME_DIR/PATH_TO_JXC_INP;	    /share/SPRKKR-8.6/kkrgen8.6MPI *JXC.inp > *JXC.out 2>&1) &'
commands = []

wd = 'for_spr'
alloys = ['Mn8Cr4Pt4', 'Ti4Fe8Cu4']
for alloy in alloys:
    groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
    for group in groups:
        orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
        for order in orders:
            path = f'{alloy}/{group}/{order}'
            commands.append(command.replace('PATH_TO_JXC_INP', path))

print(jxc_qsub.replace('COMMANDS', '\n'.join(commands)))
