# -*- coding: utf-8 -*-
# !/usr/bin/env python

import random
import os
import shutil
import warnings
from pathlib import Path
from typing import Dict

import numpy as np
import scipy.constants as constants
# from sprkkr_files import vampire_run, command_vampire
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar


da_max = 10  # lst val: 2.15
macrocell_size = np.array([36, 36, 36])
macrocell_atoms = 20000
T_min = 0
T_max = 1800
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


class Material:
    def __init__(self, idx, label, spin, element, sites):
        self.idx = idx
        self.label = label
        self.spin = spin
        self.element = element
        self.sites = sites

    def to_mat_block(self):
        if abs(self.spin) < 0.15:
            spin_str = 'non-magnetic=keep'
        else:
            spin_str = f'spin-moment = {self.spin: .6f}'
        return (
            f"material[{self.idx}]\n"
            f"material-name = {self.label}\n"
            f"element = {self.element}\n"
            f"{spin_str}\n"
            f"number-of-sites = {self.sites}\n"
        )


class Materials:
    def __init__(self, path: Path):
        if isinstance(path, str):
            path = Path(path)
        self._materials = []
        scfs = [i for i in path.glob('*SCF.out')]
        if len(scfs) > 1:
            raise Warning('Too many scf.out files, got only first')
        scf = scfs[0].read_text().split('\n')

        flag = False
        flag_mat = False
        table_start = None
        table_end = None
        for i, line in enumerate(scf):
            if 'type TXTT   NL VAL COR mesh    RMT     RWS   NAT  CONC  on sites' in line and table_start is None:
                table_start = i
            if table_start is not None and table_end is None and len(line.split()) < 5:
                table_end = i
                break

        for i in range(len(scf) - 1, 0, -1):
            if 'results extrapolated to corrected Fermi energy:' in scf[i]:
                start_mags = i
                break

        for i in range(start_mags, len(scf), 1):
            line = scf[i]
            if flag and 'sum' in line:
                magmom = abs(float(line.split()[4]))
                flag = False
                flag_mat = True

            if 'DOS      NOS     P_spin   m_spin    P_orb    m_orb     B_val      B_core' in line:
                flag = True
                it = int(scf[i - 1].split()[-2])
                name = scf[i - 1].split()[-1]
                element = name.split('_')[0]

            if flag_mat:
                for j in range(table_start, table_end):
                    if name in scf[j]:
                        sites = [int(i) for i in scf[j].split()[10:]]
                        break

                self.add(Material(idx=it,
                                  label=name,
                                  spin=magmom,
                                  element=element,
                                  sites=sites))
                flag_mat = False

    def __repr__(self):
        out = ['Materials object:']
        for mat in self._materials:
            out.append(f'idx: {mat.idx}, label: {mat.label}, spin: {mat.spin}, el: {mat.element}, sites: {mat.sites}')
        return '\n'.join(out)

    def __str__(self):
        out = ['Materials object:']
        for mat in self._materials:
            out.append(f'idx: {mat.idx}, label: {mat.label}, spin: {mat.spin}, el: {mat.element}, sites: {mat.sites}')
        return '\n'.join(out)

    def parse_scf_out(self, path: Path):

        scfs = [i for i in path.glob('*SCF.out')]
        if len(scfs) > 1:
            raise Warning('Too many scf.out files, got only first')
        scf = scfs[0].read_text().split('\n')

        flag = False
        flag_mat = False
        table_start = None
        table_end = None
        for i, line in enumerate(scf):
            if 'type TXTT   NL VAL COR mesh    RMT     RWS   NAT  CONC  on sites' in line:
                table_start = i
            if table_start is not None and table_end is None and len(line.split()) < 5:
                table_end = i

            if flag and 'sum' in line:
                magmom = float(line.split()[4])
                flag = False
                flag_mat = True

            if 'DOS      NOS     P_spin   m_spin    P_orb    m_orb     B_val      B_core' in line:
                flag = True
                it = int(scf[i - 1].split()[-2])
                name = scf[i - 1].split()[-1]
                element = name.split('_')[0]

            if flag_mat:
                for j in range(table_start, table_end):
                    if name in scf[j]:
                        sites = scf[j].split()[10:]
                        break

                self.add(Material(idx=it,
                                  label=name,
                                  spin=magmom,
                                  element=element,
                                  sites=sites))

    def add(self, material: Material):
        self._materials.append(material)

    def __iter__(self):
        return iter(self._materials)

    def __len__(self):
        return len(self._materials)

    # --- генерация .mat ---
    def to_mat(self):
        return "\n".join(m.to_mat_block() for m in self._materials)


#### Считываем тип атомов ####
def read_atoms(path='*JXC.out'):
    '''

    :param path:
    :return: indexes, concentrates
    '''
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
                IQ.append(inp[1])
                conc.append(float(inp[7]))

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


def get_ba_ca(path: Path):
    sys_path = [i for i in path.parent.glob('*.sys')]
    # print(path)
    if not sys_path:
        sys_path = [i for i in path.parents[1].glob('*.sys')]
        if not sys_path:
            print('.sys doesn\'t exist')
            raise FileExistsError
            # return '.sys doesn\'t exist'

    sys_text = sys_path[0].read_text().split('\n')
    flag = False
    for line in sys_text:
        if flag:
            ba, ca = (float(i.strip()) for i in line.split())
            return ba, ca

        if 'ratio of lattice parameters' in line:
            flag = True

    print('ba, ca didn\'t define')
    raise ValueError
    # return 'ba, ca didn\'t define'


#### Считывание файла JXC.out ####
def read_cell(path='*JXC.out'):
    '''
    :param path:
    :return: tuple of basis and coords
    '''
    # print(path)
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    ba, ca = get_ba_ca(Path(path))
    # ba, ca = 1, 1

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
    primitive_vectors[:, 1] /= ba
    primitive_vectors[:, 2] /= ca

    lat = Lattice(primitive_vectors)

    # print(lat.matrix)
    basis = np.array(Data_lattice[9:len(Data_lattice)])
    basis.shape = (int(len(basis) / 3), 3)
    basis[:, 1] /= ba
    basis[:, 2] /= ca
    # species = [i.split('_')[0] for i in read_atoms(path)[0]]
    # print(np.round(primitive_vectors * a_lat * 0.52917721090, 5), np.round(basis, 3))
    # exit()

    a = a_lat * (primitive_vectors[0, 0]**2 + primitive_vectors[0, 1]**2 + primitive_vectors[0, 2]**2)**0.5
    b = a_lat * (primitive_vectors[1, 0]**2 + primitive_vectors[1, 1]**2 + primitive_vectors[1, 2]**2)**0.5 * ba
    c = a_lat * (primitive_vectors[2, 0]**2 + primitive_vectors[2, 1]**2 + primitive_vectors[2, 2]**2)**0.5 * ca
    lens = np.array([a, b, c]) * 0.52917721090
    return np.round(primitive_vectors, 5), np.round(basis, 3), lens
    # return np.round(primitive_vectors * a_lat * 0.52917721090, 5), np.round(basis, 3)


#### Считываем обменные интегралы ####
def read_J(d_a, dr_max_count: int = 100, path: Path = '*JXC.out'):
    if isinstance(path, str):
        path = Path(path)
    temp = path.read_text().split('\n')

    if 'VERSION' in temp[13]:
        kkr_ver = temp[13].split()[3]

    if kkr_ver == '8.6.0':
        line_length = 13
        dr_pos = 10
        jij_pos = 12
        n1_pos = 4
        flag_line = 'IT   IQ   JT   JQ    N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [meV]  J_ij [eV]'
    elif kkr_ver == '7.7.3':
        line_length = 11
        dr_pos = 8
        jij_pos = 10
        n1_pos = 2
        flag_line = 'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]'
    elif kkr_ver == '6.3.1':
        line_length = 11
        dr_pos = 8
        jij_pos = 10
        n1_pos = 2
        flag_line = 'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]'
    else:
        raise RuntimeError('Unsupportable version of SPR-KKR, try to set other version')

    J = []
    dr_count = 0
    prev = 0
    # curr = 0
    if dr_max_count == -1:
        limit = False
    else:
        limit = True
    for i in range(len(temp)):
        inp = temp[i].split()

        if len(inp) == 0:
            continue

        if inp[0] == 'IQ' and inp[6] == 'JQ' and len(inp) == 12:
            IT = float(inp[5])
            JT = float(inp[11])

            if IT == 1 and JT == 1:
                print('exchange 0-0 exist')
            inp_J = temp[i + 3].split()
            if float(inp_J[dr_pos]) < d_a and abs(float(inp_J[jij_pos]) * 1000) > 1e-2:
                curr = float(inp_J[dr_pos])
                if abs(prev - curr) > 0.001:
                    # print(curr)
                    dr_count += 1

                if dr_count == dr_max_count + 1 and limit:
                    break
                    # pass
                # print(constants.e)
                # exit()
                J.append([IT - 1, JT - 1,
                          float(inp_J[n1_pos]),
                          float(inp_J[n1_pos + 1]),
                          float(inp_J[n1_pos + 2]),
                          float(inp_J[jij_pos]) * constants.e * 2])
                J.append([JT - 1, IT - 1,
                         -float(inp_J[n1_pos]),
                         -float(inp_J[n1_pos + 1]),
                         -float(inp_J[n1_pos + 2]),
                          float(inp_J[jij_pos]) * constants.e * 2])
                prev = curr
                # print(curr, float(inp_J[10]))
    # print(dr_count)
    J = np.array(J)

    sorted_idx = np.lexsort(J.T)
    sorted_data = J[sorted_idx, :]
    # sorted_data = J
    row_mask = np.append([True], np.any(np.diff(sorted_data, axis=0), 1))
    out = sorted_data[row_mask]
    # out = sorted(out, key=lambda x: (x[0], x[1]))
    J = np.array(out)

    return J


#### Считывание магнитных моментов из SCF.out ####
def read_magmom(num_atoms, composition: list, path='*SCF.out'):
    '''

    :param num_atoms:
    :param composition: list of atom_labels
    :param path:
    :return: ndarray of float magmom from SCF.out
    '''
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


def get_macrosize(cell: np.ndarray, required_n_atoms: int = macrocell_atoms):
    size = round((required_n_atoms / cell.shape[0]) ** (1/3))
    return size, size, size


def write_ucf_and_input(path: str, dr_max: int, materials: Materials):
    with open(f'{path}/vampire.UCF', "w") as file:
        file.write('# Unit cell size (Angstrom):\n')

        cell = read_cell(f'{path}/*SCF.out')

        # x_size = (cell[0][0, 0] ** 2 + cell[0][0, 1] ** 2 + cell[0][0, 2] ** 2) ** 0.5
        # y_size = (cell[0][1, 0] ** 2 + cell[0][1, 1] ** 2 + cell[0][1, 2] ** 2) ** 0.5
        # z_size = (cell[0][2, 0] ** 2 + cell[0][2, 1] ** 2 + cell[0][2, 2] ** 2) ** 0.5
        x_size = cell[2][0]
        y_size = cell[2][1]
        z_size = cell[2][2]
        # coef = np.array([[x_size, y_size, z_size],
        #                  [x_size, y_size, z_size],
        #                  [x_size, y_size, z_size]])

        file.write(f'{x_size} {y_size} {z_size}\n')
        print(f"{path.split('/')[-3]}_{path.split('/')[-2]}: {x_size} {y_size} {z_size}")
        np.savetxt(file, cell[0], fmt='%6f', header='Unit cell lattice vectors:')
        file.write('# Atoms\n')
        file.write(f'{cell[1].shape[0]} {len(materials)}\n')  # Число атомов в элементарной ячейке, число материалов
        for i in range(cell[1].shape[0]):
            for mat in materials:
                print(mat)
                if i+1 in mat.sites:
                    num_material = mat.idx - 1
            file.write(f'{i} {cell[1][i, 0]} {cell[1][i, 1]} {cell[1][i, 2]} {num_material}\n')
        file.write('# Interactions\n')
        Jij = read_J(da_max, dr_max, f'{path}/*JXC.out')
        interactions = np.column_stack((np.arange(0, Jij.shape[0]), Jij))
        file.write(f'{interactions.shape[0]} isotropic\n')
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

        # x_size = (cell[0][0, 0] ** 2 + cell[0][0, 1] ** 2 + cell[0][0, 2] ** 2) ** 0.5
        # y_size = (cell[0][1, 0] ** 2 + cell[0][1, 1] ** 2 + cell[0][1, 2] ** 2) ** 0.5
        # z_size = (cell[0][2, 0] ** 2 + cell[0][2, 1] ** 2 + cell[0][2, 2] ** 2) ** 0.5
        x_size = cell[2][0]
        y_size = cell[2][1]
        z_size = cell[2][2]

        file.write(f'dimensions:unit-cell-size-x = {x_size} !A\n')
        file.write(f'dimensions:unit-cell-size-y = {y_size} !A\n')
        file.write(f'dimensions:unit-cell-size-z = {z_size} !A\n')

        macrocell_size = get_macrosize(cell[1])

        file.write(f'dimensions:system-size-x = {0.1 * x_size * macrocell_size[0]} !nm\n')
        file.write(f'dimensions:system-size-y = {0.1 * y_size * macrocell_size[1]} !nm\n')
        file.write(f'dimensions:system-size-z = {0.1 * z_size * macrocell_size[2]} !nm\n')

        file.write('#------------------------------------------\n')
        file.write('# Exchange parameter:\n')
        file.write('#------------------------------------------\n')
        file.write('#exchange:ab-initio = true\n')

        file.write('#------------------------------------------\n')
        file.write('# Simulation attributes:\n')
        file.write('#------------------------------------------\n')
        file.write(f'sim:minimum-temperature = {T_min}\n')
        file.write(f'sim:maximum-temperature = {T_max}\n')
        file.write(f'sim:temperature-increment = {T_step}\n')

        file.write(f'sim:equilibration-time-steps = {MC_step // 10}\n')
        file.write(f'sim:loop-time-steps = {MC_step}\n')
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


def write_materials(materials: Materials):

    mat_file = [f'material:num-materials = {len(materials)}']
    for mat in materials:
        mat_file.append(f'#---------------------------------------------------')
        mat_file.append(f'# Material {mat.idx}')
        mat_file.append(f'#---------------------------------------------------')
        mat_file.append(f'material[{mat.idx}]:material-name={mat.label}')
        mat_file.append(f'material[{mat.idx}]:damping-constant=1.0')
        if abs(mat.spin) < 0.15:
            mat_file.append(f'material[{mat.idx}]:non-magnetic=keep')
        else:
            mat_file.append(f'material[{mat.idx}]:atomic-spin-moment={mat.spin} !muB')
        mat_file.append(f'material[{mat.idx}]:uniaxial-anisotropy-constant=0.0')
        mat_file.append(f'material[{mat.idx}]:material-element={mat.element}')
        mat_file.append(f'material[{mat.idx}]:initial-spin-direction = 0.0,0.0,1.0')
        mat_file.append(f'material[{mat.idx}]:uniaxial-anisotropy-direction = 0.0 , 0.0, 1.0')

    return '\n'.join(mat_file)


def generate_vampire_inputs_recursive(root_path: Path, depth: int, dr_max: int):
    '''

    :param root_path:
    :param depth: different between length of root and target path
    :param dr_max: count of coordination spheres
    :return:
    '''
    # depth = 2

    for system_path in root_path.rglob('*JXC_auto.out'):
        # if 'Al/L21' not in system_path.as_posix():
        #     continue
        if len(system_path.relative_to(root_path).parts) != depth + 1:
            continue
        # print(system_path)
        path = system_path.parent
        # if dr_max:
        #     path = path / dr_max

        # if (path / 'vampire').exists():
        #     shutil.rmtree(f'{path}/vampire')
        #
        # continue
        vamp_dir_name = 'vampire_new'
        # vamp_dir_name = str(dr_max)
        if not (path / vamp_dir_name).exists():
            (path / vamp_dir_name).mkdir()

        for file in ['*JXC.out', '*SCF.out']:
            src = path / file
            dst = path / vamp_dir_name / file
            shutil.copy2(src, dst)
        path = path / vamp_dir_name
        # shutil.copy2('/home/buche/VaspTesting/Danil/magnetocaloric_nn/vampire/vampire.mat',
        #              f'{path}/vampire.mat')

        # print(path)
        labels = read_atoms(f'{path}/*JXC.out')[0]
        rwss = read_atoms(f'{path}/*JXC.out')[0]
        num = len(read_atoms(f'{path}/*JXC.out')[0])
        mags = read_magmom(num, labels, f'{path}/*SCF.out')
        mags = [abs(float(mag)) for mag in mags]
        # print(mags)

        materials = Materials(path)
        # print(materials)
        (path / 'vampire.mat').write_text(write_materials(materials))
        # print(write_materials(materials))
        write_ucf_and_input(path.as_posix(), dr_max, materials)


def get_curve_recursively(root: Path):
    '''
    get Magnetisation against temperature curve
    from certain folder recursively (look for "output" files
    :param root: root folder
    :return: dict[full_path_str, np.ndarray curve]
    '''
    curves = {}
    for path in root.rglob('output'):
        if not check_status_vampire(path.parent):
            continue

        if 'test' in path.as_posix():
            continue
        out = path.read_text().split('\n')
        flag = False
        curve = []
        for line in out[:-1]:
            if len(line.split()) <= 2:
                continue
            if line[0] == '0':
                flag = True
                curve.append([float(i) for i in line.split()[:2]])
            if flag:
                curve.append([float(i) for i in line.split()[:2]])
        curves[path.parent.as_posix()] = np.array(curve)
        # print(curve)
    # return dict(sorted(curves.items(), key=lambda item: int(item[0].split('/')[-1])))
    return curves


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.interpolate import make_interp_spline, PchipInterpolator


def check_status_vampire(p: Path):
    path = (p / 'log')
    if path.exists():
        temp = path.read_text()
        if 'Simulation ended gracefully.' in temp:
            return True
        else:
            return False
    return False


def plot_exchange_from_jxc(path: Path):
    path = path / '*JXC.out'
    temp = path.read_text().split('\n')
    jxc = []
    prev_dr = 0
    flag = False
    for line in temp:
        if flag and len(line.split()) == 11:
            # print(line.split())
            dr = float(line.split()[8])
            Jij = float(line.split()[10]) * 1000  # meV

            if dr == prev_dr:
                continue
            else:
                prev_dr = dr

            jxc.append([dr, Jij])
            flag = False

        if 'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]' in line:
            flag = True
    jxc = np.array(jxc)
    # print(jxc)

    f = plt.figure(figsize=(6, 4.5))
    ax = f.add_subplot(111)
    plt.tick_params(axis='both',  # Применяем параметры к обеим осям
                    which='major',  # Применяем параметры к вспомогательным делениям
                    direction='in',  # Рисуем деления внутри и снаружи графика
                    # length = 10,    #  Длина делений
                    # width = 2,     #  Ширина делений
                    # color = 'm',    #  Цвет делений
                    pad=10,  # Расстояние между черточкой и ее подписью
                    labelsize=12,  # Размер подписи
                    labelcolor='k',  # Цвет подписи
                    bottom=True,  # Рисуем метки снизу
                    top=True,  # сверху
                    left=True,  # слева
                    right=True)  # и справа

    ax.set_xlabel(r'$d/a$')
    ax.set_ylabel(r'$J_{ij}$, мэВ')

    ax.plot(jxc[:, 0], jxc[:, 1], '-',
            linewidth=1,
            label='Fe-Fe',
            marker='o',
            mec='k',
            mew=0.5,
            markersize=7)

    jxc = path
    jxc = reversed(jxc.read_text().split('\n'))
    for line in jxc:
        if 'Curie temperature within mean field approximation  T_C' in line:
            T_curie = line.split()[-2]

    f.text(0.70, 0.5, r'$T_C^{MFA}$ = ' + T_curie.replace('.', ',') + ' K', backgroundcolor='#C8C8C8')
    ax.legend(loc='best')
    ax.grid()
    plt.show()
    # if float(T_curie) == 0:
    #     continue
    # ax.set_xlim(0, 7)
    f.savefig(path.parent ,
              transparent=False,
              bbox_inches='tight',
              dpi=300,
              pad_inches=0.01)


def plot_all_mags(curves: Dict[str, np.ndarray], save_path: Path):
    fig, ax = plt.subplots(figsize=(6, 4))

    for path, curve in curves.items():
        T = curve[1:, 0]
        M = curve[1:, 1]
        T_smooth = np.linspace(T.min(), T.max(), 500)
        # y_smooth = X_Y_Spline(x_smooth)
        print(T)
        pchip = PchipInterpolator(T, M)
        M_smooth = pchip(T_smooth)

        label = path.split('/')[-1]
        parts = path.split('/')
        label = f'{parts[-3]}_{parts[-2]}'
        ax.plot(T_smooth, M_smooth, lw=1, label=label)

    ax.set_xlabel("Температура (K)", fontsize=12)
    ax.set_ylabel("Намагниченность (норм.)", fontsize=12)
    # ax.set_title(rf"Температура Кюри — критическая аппроксимация"
    #              "\n" rf"$M = A\,(1 - T/T_C)^{{\beta}}$", fontsize=13)
    ax.legend(fontsize=6, ncol=2)
    ax.grid(True, alpha=0.4)
    ax.set_ylim(bottom=-0.02)
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"График сохранён: {save_path}")


def get_tc_mfa(path_to_jxc: Path, spheres: int = None, kkr_ver: str = None):
#     T_C = 2/3*sum_J_AA [Joul] /kB [Joul/K]

    path = path_to_jxc / '*JXC.out'
    temp = path.read_text().split('\n')

    if spheres is None:
        temp.reverse()
        for line in temp:
            if 'Curie temperature within mean field approximation' in line:
                Tc = float(line.split()[-2])
                return Tc

    jxc = []
    flag = False
    count_spheres = 0
    prev_dr = 0

    if 'VERSION' in temp[13]:
        kkr_ver = temp[13].split()[3]

    if kkr_ver == '8.6.0':
        line_length = 13
        dr_pos = 10
        jij_pos = 11
        flag_line = 'IT   IQ   JT   JQ    N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [meV]  J_ij [eV]'
    elif kkr_ver == '7.7.3':
        line_length = 11
        dr_pos = 8
        jij_pos = 10
        flag_line = 'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]'
    elif kkr_ver == '6.3.1':
        line_length = 11
        dr_pos = 8
        jij_pos = 10
        flag_line = 'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]'
    else:
        raise RuntimeError('Unsupportable version of SPR-KKR, try to set other version')

    for line in temp:
        if flag and len(line.split()) == line_length:
            # print(line.split())
            dr = float(line.split()[dr_pos])
            Jij = float(line.split()[jij_pos])  # meV

            if dr != prev_dr:
                count_spheres += 1
                prev_dr = dr

            if count_spheres > spheres and spheres:
                break

            jxc.append([dr, Jij])
            flag = False

        if flag_line in line:
            flag = True
    jxc = np.array(jxc)

    k = 8.617333262e-5
    # print(np.sum(jxc[:, 1]))
    # exit()
    Tc = 2/3 * np.sum(jxc[:, 1]) / k
    return Tc


if __name__ == '__main__':
    # wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/Fe/')
    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ')

    # curves = get_curve_recursively(wd)

    # curves = dict(sorted(curves.items(), key=lambda item: item[0].split('/')[-3]))
    # plot_all_mags(curves, (wd / 'All_mags.png'))
    generate_vampire_inputs_recursive(wd, 2, 1)
    # J = read_J(da_max, path=wd / 'Al' / 'L21' / '*JXC.out')
    # for line in J:
    #     if line[0] == 0 and line[1] == 0:
    #         print('good exchanfe')
    # print(J)
    # finds = sorted([i.parent.as_posix() for i in wd.rglob('*.UCF')])
    # (wd / 'vampire_work_paths').write_text('\n'.join(finds))
    exit()
    # plot_exchange_from_jxc(wd)
    # exit()
    # print(len(wd.relative_to(wd.parent).parts))
    #

    # exit()
    # jij = read_J(da_max, 2, 'JXC.out')
    # print(np.round(jij[:-1]), jij[-1])
    # for i in jij:
    #     print(np.round(i[:-1]), i[-1])
    # wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ')
    # wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/Fe_new_parser/')
    #
    # for i in range(1, 49):
    #     stat = check_status_vampire(wd / str(i))
    #     print(f'{i:<2} - {stat}')
    # exit()

    # wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/Fe_new_parser/')
    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/Fe/')
    # beta_type = 'fixed'
    # beta_type = 'free'
    # # generate_vampire_inputs_recursive(wd, 0, 45)
    n_t = []
    for i in range(1, 24):
        # generate_vampire_inputs_recursive(wd, 0, i)
        # print((wd / str(i)).as_posix())
        path = wd / str(i) / 'log'
        log = path.read_text().split('\n')

        for line in log:
            if 'Simulation run time' in line:
                n_t.append([i, float(line.split()[-1])])
                break
    print(n_t)
    t_n = np.array(n_t)
    print(t_n)
    exit()
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.plot(t_n[:, 0], t_n[:, 1])
    ax.set_xticks(np.arange(0, 45, 1))
    ax.set_xlim(0, 45)
    plt.show()
    exit()
    # generate_run_recursively(wd)
    # exit()
    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/')

    curves = get_curve_recursively(wd)
    # curves = dict(sorted(curves.items(), key=lambda item: int(item[0].split('/')[-1])))

    plot_all_mags(curves, (wd / 'All_mags.png'))

    # test_curve = list(curves.values())[0]
    # test_curve[:, 1] = test_curve[:, 1]**(1/0.365)
    # x = list(curves.values())[0][:, 0]
    # y = list(curves.values())[0][:, 1]
    # y = test_curve[:, 1]
    # idx = np.where(y > 0.05)
    # x = x[idx]
    # y = y[idx]
    # test_coefs = np.polyfit(x, y, 1)
    #
    # x_fit = np.linspace(0, 1400, 500)
    # y_fit = test_coefs[0] * x_fit + test_coefs[1]
    #
    # tc = x_fit[np.where(y_fit > 0)][-1]
    # print(tc)

    # for i in test_curve:
    #     if i[1] < 0.01:
    #         print(i)
    #         break
    # exit()
    # f, ax = plt.subplots()
    # ax.plot(test_curve[:, 0], test_curve[:, 1], 'or-')
    # ax.plot(x_fit, y_fit, 'b-')
    # plt.grid()
    # plt.show()

    # exit()
    # for i in curves.keys():
    #     print(i)
    # exit()
    # # print(list(curves.values())[0])
    # print(dict(sorted(curves.items(), key=lambda item: int(item[0].split('/')[-1]))))
    # exit()
    tc = []
    tc_mean = []
    from Critical_fit_TC import curie_fit_fixed_beta
    for i, (path, curve) in enumerate(curves.items(), start=1):
        # curie_critical_fit(curve, f'{path}/curve.png')
        # tc.append([i, curie_max_derivative(curve, f'{path}/curve.png', plot=False)])
        print(f'{"": #^100}')
        tc_mfa = get_tc_mfa(Path(path), None, None)
        tc_fit = curie_fit_fixed_beta(curve, tc_mfa, f'{path}/curve.png')

        # if tc_fit[1] > 0.4:
        #     continue
        tc.append([i, tc_fit])
        # mean_tc = sum([obj[1] for obj in tc])/len(tc)
        # tc_mean.append([i, mean_tc])

    tc = np.array(tc_mean)
    # tc = get_mean_field_Tc(wd)
    print(tc)
    exit()

    # X_Y_Spline = make_interp_spline(tc[:-1, 0], tc[:-1, 1], k=len(tc)-1)
    x_smooth = np.linspace(tc[:, 0].min(), tc[:, 0].max(), 500)
    # y_smooth = X_Y_Spline(x_smooth)

    pchip = PchipInterpolator(tc[:, 0], tc[:, 1])
    y_smooth = pchip(x_smooth)

    fig, ax = plt.subplots(figsize=(8, 6))

    # ax.plot(tc[:-1, 0], tc[:-1, 1], marker='o', linewidth=2, markersize=6)
    plt.plot(x_smooth, y_smooth, linewidth=2, zorder=5)
    # plt.plot(tc[:-1, 0], tc[:-1, 1], '-', linewidth=2, zorder=5)
    plt.scatter(tc[:, 0], tc[:, 1], color='red', linewidth=2, marker='o', zorder=10)

    ax.set_xlabel("Количество координационных сфер", fontsize=13)
    ax.set_ylabel("Средняя температура Кюри, $T_C$ (K)", fontsize=13)
    ax.set_title("Monte-Carlo", fontsize=14)

    ax.set_xticks(np.arange(0, 50, 5))
    ax.grid(True, linestyle='--', alpha=0.6, zorder=1)
    ax.tick_params(labelsize=11)

    fig.tight_layout()
    fig.savefig(f"{wd}/Monte_Carlo_mean_tc.png", dpi=300)
    plt.show()

    # dst = wd / f'curves'
    # dst.mkdir(exist_ok=True)
    # for i in wd.rglob('curve.png'):
    #     shutil.copy2(i, dst / f'{i.parent.stem}_mag.png')

    print('end')
    # plt.close(fig)
    # for i in range(1, 45):
    #     generate_vampire_inputs_recursive(wd, 0, i)





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