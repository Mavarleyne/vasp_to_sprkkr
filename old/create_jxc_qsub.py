import os
import numpy as np
from pathlib import Path

from ase2sprkkr import SPRKKRAtoms
from ase2sprkkr.common.unique_values import UniqueValuesMapping
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
import shutil as sh
from ase import Atoms

jxc_qsub = '''
#!/bin/bash
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
command = '(cd $HOME_DIR/PATH_TO_JXC_INP	    /share/SPRKKR-8.6/kkrgen8.6MPI *JXC.inp > *JXC.out 2>&1) &'
commands = []
index = 0

path = Path('for_spr')
for i in ['216', '225']:
    for j in os.listdir(Path(f'{path}/{i}/')):
        for k in os.listdir(Path(f'{path}/{i}/{j}')):
            p = f'{i}/{j}/{k};'
            commands.append(command.replace('PATH_TO_JXC_INP', f'{p:30}'))
            index += 1
            if index == 100:
                break
            # print(command.replace('PATH_TO_JXC_INP', f'{i}/{j}/{k}/0'))

print(jxc_qsub.replace('COMMANDS', '\n'.join(commands)))