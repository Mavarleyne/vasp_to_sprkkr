import os

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
