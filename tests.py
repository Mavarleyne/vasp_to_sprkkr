# -*- coding: utf-8 -*-
# !/usr/bin/env python

import random
import os
import numpy as np
import scipy.constants as constants

da_max = 2.5
macrocell_size = np.array([4, 4, 4])
T_min = 0
T_max = 1500
T_step = 10
MC_step = 1000


#### Считывание файла JXC.out ####
def read_cell(path='JXC.out'):
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

    basis = np.array(Data_lattice[9:len(Data_lattice)])
    basis.shape = (int(len(basis) / 3), 3)

    return primitive_vectors * a_lat * 0.52917721090, basis


#### Считываем тип атомов ####
def read_atoms(path='JXC.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    IQ = []
    conc = []
    flag = False
    for line in lines:
        inp = line.split()
        if flag == True and len(inp) >= 9:
            for i in range(int(inp[len(inp) - 1]) - int(inp[8]) + 1):
                IQ += [inp[1]]
                conc += [float(inp[7])]
        if len(inp) > 3:
            if inp[0] == 'type' and inp[1] == 'TXTT' and inp[2] == 'NL':
                flag = True
        if len(inp) == 1:
            if inp[0] == 'mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm':
                break
    conc = np.array(conc)
    return IQ, conc


#### Считываем обменные интегралы ####
def read_J(d_a, path='JXC.out'):
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
def read_magmom(num_atoms, path='SCF.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    magmom = np.zeros(num_atoms)
    flag = False
    n = 0
    for i in range(num_atoms):
        for line in lines:
            inp = line.split()
            if len(inp) > 5 and inp[1] == 'E=' and inp[4] == 'IT=' and int(inp[5]) == i + 1:
                flag = True
            if len(inp) > 9 and flag == True and inp[0] == 'sum':
                magmom[n] = float(inp[4])
                flag = False
        n = n + 1
    return magmom


file = open('vampire.UCF', "w")
file.write('# Unit cell size (Angstrom):\n')
file.write('1.0 1.0 1.0\n')
cell = read_cell()
np.savetxt(file, cell[0], fmt='%6f', header='Unit cell lattice vectors:')
file.write('# Atoms\n')
file.write('{} {}\n'.format(cell[1].shape[0], cell[1].shape[0]))  # Число атомов в элементарной ячейке, число материалов
for i in range(cell[1].shape[0]):
    file.write('{} {} {} {} {}\n'.format(i, cell[1][i, 0], cell[1][i, 1], cell[1][i, 2], i))
file.write('# Interactions\n')

interactions = np.column_stack((np.arange(0, read_J(da_max).shape[0]), read_J(da_max)))
file.write('{} isotropic\n'.format(interactions.shape[0]))
np.savetxt(file, interactions, fmt='%d %d %d %d %d %d %.4e')
file.close()

file = open('input', "w")
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
file.close()
