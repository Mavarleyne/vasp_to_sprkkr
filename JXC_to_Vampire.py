# -*- coding: utf-8 -*-
# !/usr/bin/env python

import random
import os
import shutil

import numpy as np
import scipy.constants as constants
# from sprkkr_files import vampire_run, command_vampire
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar


da_max = 2.15
macrocell_size = np.array([4, 4, 4])
T_min = 0
T_max = 1000
T_step = 10
MC_step = 1000

vampire_run = '''#!/bin/bash
#PBS -d .
#PBS -l nodes=1:ppn=20 #:iband
#PBS -N NAME
#PBS -j oe
#PBS -l walltime=2000:00:00

currDir=`pwd`

LD_LIBRARY_PATH=/share/intel/mkl/lib/intel64/:$LD_LIBRARY_PATH
. /share/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpivars.sh
export I_MPI_FALLBACK_DEVICE=disable
export I_MPI_FABRICS=shm #:ofa
export I_MPI_PIN=disable
export LD_LIBRARY_PATH

ulimit -s unlimited

COMMANDS


wait'''

command_vampire = '(cd PATH_TO_VAMP_INP;      /share/vampire/bin/vampire-serial-intel) &'

mat_sample = '''#---------------------------------------------------
# Material {}
#---------------------------------------------------
material[index]:material-name={}
material[index]:damping-constant=1.0
material[index]:atomic-spin-moment={} !muB
material[index]:uniaxial-anisotropy-constant=0.0
material[index]:material-element={}
material[index]:initial-spin-direction = 0.0,0.0,1.0
material[index]:uniaxial-anisotropy-direction = 0.0 , 0.0, 1.0
'''


#### Считываем тип атомов ####
def read_atoms(path='*JXC.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    IQ = []
    conc = []
    flag = False
    for line in lines:
        inp = line.split()
        if flag and len(inp) >= 9:
            # for i in range(int(inp[len(inp) - 1]) - int(inp[8]) + 1):

            for i in range(int(inp[nat_i])):
                IQ += [inp[1]]
                conc += [float(inp[7])]

        if len(inp) > 3:
            if inp[0] == 'type' and inp[1] == 'TXTT' and inp[2] == 'NL':
                for i in range(len(inp)):
                    if inp[i] == 'NAT':
                        nat_i = i
                flag = True
        if len(inp) == 1:
            if inp[0] == 'mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm':
                break
    conc = np.array(conc)
    return IQ, conc


#### Считывание файла JXC.out ####
def read_cell(path='*JXC.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    #### Считываем вектора трансляций, базис ####
    Data_lattice = []
    flag = False
    for line in lines:
        inp = line.split()
        if len(inp) == 0:
            continue
        if flag == True and len(inp) == 5 and inp[0] == '(':
            Data_lattice += [float(inp[1].replace(',', '')), float(inp[2].replace(',', '')), float(inp[3])]
        if inp[0] == '<INIT_MOD_LATTICE>' and len(inp) >= 1:
            flag = True
        if inp[len(inp) - 1] == '2*pi/a' and len(inp) > 1:
            flag = False
        if inp[0] == 'lattice' and inp[1] == 'constant' and inp[2] == 'ALAT' and len(inp) >= 1:
            a_lat = float(inp[3])

    primitive_vectors = np.array(Data_lattice[:9])
    primitive_vectors.shape = (3, 3)

    # lat = Lattice(primitive_vectors)
    # print(lat.angles)
    # print(lat.matrix)
    basis = np.array(Data_lattice[9:len(Data_lattice)])
    basis.shape = (int(len(basis) / 3), 3)
    # species = [i.split('_')[0] for i in read_atoms(path)[0]]
    # print(primitive_vectors)
    # print(species)
    # print(basis)
    # struc = Structure(lattice=primitive_vectors * a_lat,
    #                   species=species,
    #                   coords=basis,
    #                   to_unit_cell=True)
    # sga = SpacegroupAnalyzer(struc).get_space_group_number()
    # print(sga)
    # print(Poscar(struc))
    # print(SpacegroupAnalyzer(Poscar(struc).structure).get_space_group_number())
    # exit()

    return primitive_vectors * a_lat * 0.52917721090, basis
    # return struc.lattice.matrix * a_lat * 0.52917721090,  struc.frac_coords


#### Считываем обменные интегралы ####
def read_J(d_a, path='*JXC.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    J = []
    for i in range(len(lines)):
        inp = lines[i].split()

        if len(inp) == 0:
            continue

        if inp[0] == 'IQ' and inp[6] == 'JQ' and len(inp) == 12:
            IQ = float(inp[5]);
            JQ = float(inp[11])
            inp_J = lines[i + 3].split()
            if float(inp_J[8]) < d_a and abs(float(inp_J[10]) * 1000) > 0.01:
                J += [[IQ - 1, JQ - 1, float(inp_J[2]), float(inp_J[3]), float(inp_J[4]),
                       float(inp_J[10]) * constants.e * 2]]
                J += [[JQ - 1, IQ - 1, -float(inp_J[2]), -float(inp_J[3]), -float(inp_J[4]),
                       float(inp_J[10]) * constants.e * 2]]

    J = np.array(J)

    sorted_idx = np.lexsort(J.T)
    sorted_data = J[sorted_idx, :]
    row_mask = np.append([True], np.any(np.diff(sorted_data, axis=0), 1))
    out = sorted_data[row_mask]
    out = sorted(out, key=lambda x: (x[0], x[1]))
    J = np.array(out)

    return J


#### Считывание магнитных моментов из SCF.out ####
def read_magmom(num_atoms, composition: list, path='*SCF.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    magmom = np.zeros(num_atoms)
    flag = False
    n = 0
    for i, atom in zip(range(num_atoms), composition):
        for line in lines:
            inp = line.split()
            if len(inp) > 5 and inp[1] == 'E=' and inp[4] == 'IT=' and inp[6] == atom:
                flag = True
            if len(inp) > 9 and flag and inp[0] == 'sum':
                magmom[n] = float(inp[4])
                flag = False
        n = n + 1
    return magmom


def write_ucf_and_input(path: str):
    with open(f'{path}/vampire.UCF', "w") as file:
        file.write('# Unit cell size (Angstrom):\n')
        file.write('1.0 1.0 1.0\n')
        cell = read_cell(f'{path}/*SCF.out')
        np.savetxt(file, cell[0], fmt='%6f', header='Unit cell lattice vectors:')
        file.write('# Atoms\n')
        file.write('{} {}\n'.format(cell[1].shape[0], cell[1].shape[0]))  # Число атомов в элементарной ячейке, число материалов
        for i in range(cell[1].shape[0]):
            file.write('{} {} {} {} {}\n'.format(i, cell[1][i, 0], cell[1][i, 1], cell[1][i, 2], i))
        file.write('# Interactions\n')

        interactions = np.column_stack((np.arange(0, read_J(da_max, f'{path}/*JXC.out').shape[0]), read_J(da_max, f'{path}/*JXC.out')))
        file.write('{} isotropic\n'.format(interactions.shape[0]))
        np.savetxt(file, interactions, fmt='%d %d %d %d %d %d %.4e')

    with open(f'{path}/input', "w") as file:
        file.write('#------------------------------------------\n')
        file.write('# Creation attributes:\n')
        file.write('#------------------------------------------\n')
        file.write('create:full\n')
        file.write('create:periodic-boundaries-x\n')
        file.write('create:periodic-boundaries-y\n')
        file.write('create:periodic-boundaries-z\n')
        file.write('#------------------------------------------\n')
        file.write('material:file=vampire.mat\n')
        file.write('material:unit-cell-file = "vampire.UCF"\n')
        file.write('#------------------------------------------\n')
        file.write('# System Dimensions:\n')
        file.write('#------------------------------------------\n')

        x_size = (cell[0][0, 0] ** 2 + cell[0][0, 1] ** 2 + cell[0][0, 2] ** 2) ** 0.5
        y_size = (cell[0][1, 0] ** 2 + cell[0][1, 1] ** 2 + cell[0][1, 2] ** 2) ** 0.5
        z_size = (cell[0][2, 0] ** 2 + cell[0][2, 1] ** 2 + cell[0][2, 2] ** 2) ** 0.5

        file.write('dimensions:unit-cell-size-x = {} !A\n'.format(x_size))
        file.write('dimensions:unit-cell-size-y = {} !A\n'.format(y_size))
        file.write('dimensions:unit-cell-size-z = {} !A\n'.format(z_size))

        file.write('dimensions:system-size-x = {} !nm\n'.format(0.1 * x_size * macrocell_size[0]))
        file.write('dimensions:system-size-y = {} !nm\n'.format(0.1 * y_size * macrocell_size[1]))
        file.write('dimensions:system-size-z = {} !nm\n'.format(0.1 * z_size * macrocell_size[2]))
        file.write('#------------------------------------------\n')
        file.write('# Simulation attributes:\n')
        file.write('#------------------------------------------\n')
        file.write('sim:minimum-temperature = {}\n'.format(T_min))
        file.write('sim:maximum-temperature = {}\n'.format(T_max))
        file.write('sim:temperature-increment = {}\n'.format(T_step))

        file.write('sim:equilibration-time-steps = {}\n'.format(MC_step))
        file.write('sim:loop-time-steps = {}\n'.format(MC_step))
        file.write('sim:time-steps-increment = 1\n')
        file.write('#------------------------------------------\n')
        file.write('# Program and integrator details\n')
        file.write('#------------------------------------------\n')
        file.write('sim:program=curie-temperature\n')
        file.write('sim:integrator = monte-carlo\n')
        file.write('#------------------------------------------\n')
        file.write('# Data output\n')
        file.write('#------------------------------------------\n')
        file.write('#output:real-time\n')
        file.write('output:temperature\n')
        file.write('#output:material-magnetisation\n')
        file.write('output:magnetisation-length\n')
        file.write('output:mean-total-energy\n')


def generate_vampire_inputs(wd: str):
    alloys = [i for i in os.listdir(f'{wd}') if os.path.isdir(f'{wd}/{i}')]
    for alloy in alloys:
        groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
        for group in groups:
            if alloy == 'Ta4Ti8Mo4' and group == '216':
                continue
            orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
            for order in orders:
                path = f'{wd}/{alloy}/{group}/{order}'
                # if os.path.isfile(f'{path}/vampire/output'):
                #     continue

                if not os.path.isdir(f'{path}/vampire'):
                    os.mkdir(f'{path}/vampire')

                for file in ['*JXC.out', '*SCF.out']:
                    src = f'{path}/{file}'
                    dst = f'{path}/vampire/{file}'
                    shutil.copy2(src, dst)
                path = f'{wd}/{alloy}/{group}/{order}/vampire'
                shutil.copy2('/home/buche/VaspTesting/Danil/magnetocaloric_nn/vampire/vampire.mat',
                             f'{path}/vampire.mat')
                # print(len(read_atoms(f'{path}/*JXC.out')[0]))
                # print(read_atoms(f'{path}/*JXC.out')[0])
                # print(read_atoms(f'{path}/*JXC.out')[1])
                # print(read_magmom(num, f'{path}/*SCF.out'))
                print(path)
                labels = read_atoms(f'{path}/*JXC.out')[0]
                rwss = read_atoms(f'{path}/*JXC.out')[0]
                num = len(read_atoms(f'{path}/*JXC.out')[0])
                mags = read_magmom(num, labels, f'{path}/*SCF.out')
                mags = [abs(float(mag)) for mag in mags]
                # print(mags)

                with open(f'{path}/vampire.mat', 'w') as f:
                    f.write(f'material:num-materials = {num}\n')
                    for i in range(num):
                        # mat =  ''

                        mat = mat_sample.format(i + 1, labels[i], mags[i], labels[i].split('_')[0]).replace('[index]', f'[{i+1}]')
                        if mags[i] < 0.1:
                            mat = mat.replace(f'material[{i+1}]:atomic-spin-moment', f'# material[{i+1}]:atomic-spin-moment')

                        f.write(mat)
                # for i in range(num):
                #     print(mat_sample.format(i+1, labels[i], mags[i], labels[i].split('_')[0]), end='')
                write_ucf_and_input(path)


def generate_run(wd: str):
    comms = ''
    alloys = [i for i in os.listdir(f'{wd}') if os.path.isdir(f'{wd}/{i}')]
    for alloy in alloys:
        groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
        for group in groups:
            if alloy == 'Ta4Ti8Mo4' and group == '216':
                continue
            orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
            for order in orders:
                path = f'{wd}/{alloy}/{group}/{order}/vampire'
                comms += command_vampire.replace('PATH_TO_VAMP_INP', f'{path};') + '\n'

    with open(f'{wd}/vampire_qsub', 'w') as f:
        f.write(vampire_run.replace('COMMANDS', comms))
    # print(vampire_run.replace('COMMANDS', comms))


def get_curve(wd):
    alloys = [i for i in os.listdir(f'{wd}') if os.path.isdir(f'{wd}/{i}')]
    for alloy in alloys:
        groups = [i for i in os.listdir(f'{wd}/{alloy}') if os.path.isdir(f'{wd}/{alloy}/{i}')]
        for group in groups:
            if alloy == 'Ta4Ti8Mo4' and group == '216':
                continue
            orders = [i for i in os.listdir(f'{wd}/{alloy}/{group}') if os.path.isdir(f'{wd}/{alloy}/{group}/{i}')]
            for order in orders:
                path = f'{wd}/{alloy}/{group}/{order}/vampire'
                if not os.path.isfile(f'{path}/output'):
                    continue
                with open(f'{path}/output') as f:
                    temp = f.readlines()
                    flag = False
                    curve = []
                    for line in temp:
                        if line[0] == '0':
                            flag = True
                            curve.append([float(i) for i in line.split()[:2]])
                        if flag:
                            curve.append([float(i) for i in line.split()[:2]])
                    curve = np.array(curve)
                    print(curve)


if __name__ == '__main__':
    wd = '/home/buche/VaspTesting/Danil/magnetocaloric_nn/for_spr'
    generate_vampire_inputs(wd)
    # read_cell('/home/buche/VaspTesting/Danil/magnetocaloric_nn/for_spr/Ti4Fe8Cu4/119/FiM/vampire/JXC.out')
    generate_run(wd)
    # get_curve(wd)

'''
#!/bin/bash
#PBS -d .
#PBS -l nodes=1:ppn=5 #:iband
#PBS -N Vampire_for_spr
#PBS -j oe
#PBS -l walltime=2000:00:00

currDir=`pwd`

LD_LIBRARY_PATH=/share/intel/mkl/lib/intel64/:$LD_LIBRARY_PATH
. /share/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpivars.sh
export I_MPI_FALLBACK_DEVICE=disable
export I_MPI_FABRICS=shm #:ofa
export I_MPI_PIN=disable
export LD_LIBRARY_PATH

ulimit -s unlimited

(cd /home/buche/VaspTesting/Danil/magnetocaloric_nn/for_spr/Ti4Fe8Cu4/119/FiM/vampire;     /share/vampire/bin/vampire-serial-intel) &
(cd /home/buche/VaspTesting/Danil/magnetocaloric_nn/for_spr/Sc4Co4Ni8/139/FiM/vampire;     /share/vampire/bin/vampire-serial-intel) &
(cd /home/buche/VaspTesting/Danil/magnetocaloric_nn/for_spr/Mn8Cr4Pt4/225/FiM/vampire;     /share/vampire/bin/vampire-serial-intel) &
(cd /home/buche/VaspTesting/Danil/magnetocaloric_nn/for_spr/Mn8Cr4Pt4/139/FiM/vampire;     /share/vampire/bin/vampire-serial-intel) &
(cd /home/buche/VaspTesting/Danil/magnetocaloric_nn/for_spr/Ti4Ni8Ru4/119/FM/vampire;   /share/vampire/bin/vampire-serial-intel) &


wait
'''