import shutil
from datetime import datetime
from pathlib import Path
from System import generate_sys
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Potcar, PotcarSingle, Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SgA
from brave_from_pearson import Pearson, international_numbers_to_AP
from collections import Counter
import numpy as np
import os
from Wigner_Seitz_radius import get_rws


class Potential:
    def __init__(self, path: str | Path, structure: Structure):
        '''

        :param structure:
        '''
        self.structure = structure
        self.conventional_structure = SgA(structure).get_conventional_standard_structure()
        self.primitive_structure = SgA(structure).get_primitive_standard_structure()


    __paws = {
        'H': 'H', 'He': 'He', 'Li': 'Li_sv', 'Be': 'Be', 'B': 'B', 'C': 'C', 'N': 'N', 'O': 'O', 'F': 'F', 'Ne': 'Ne',
        'Na': 'Na_pv', 'Mg': 'Mg', 'Al': 'Al', 'Si': 'Si', 'P': 'P', 'S': 'S', 'Cl': 'Cl', 'Ar': 'Ar', 'K': 'K_sv',
        'Ca': 'Ca_sv', 'Sc': 'Sc_sv', 'Ti': 'Ti_sv', 'V': 'V_sv', 'Cr': 'Cr_pv', 'Mn': 'Mn_pv', 'Fe': 'Fe', 'Co': 'Co',
        'Ni': 'Ni', 'Cu': 'Cu', 'Zn': 'Zn', 'Ga': 'Ga_d', 'Ge': 'Ge_d', 'As': 'As', 'Se': 'Se', 'Br': 'Br', 'Kr': 'Kr',
        'Rb': 'Rb_sv', 'Sr': 'Sr_sv', 'Y': 'Y_sv', 'Zr': 'Zr', 'Nb': 'Nb_sv', 'Mo': 'Mo_sv', 'Tc': 'Tc_pv',
        'Ru': 'Ru_pv',
        'Rh': 'Rh_pv', 'Pd': 'Pd', 'Ag': 'Ag', 'Cd': 'Cd', 'In': 'In_d', 'Sn': 'Sn_d', 'Sb': 'Sb', 'Te': 'Te', 'I': 'I',
        'Xe': 'Xe', 'Cs': 'Cs_sv', 'Ba': 'Ba_sv', 'La': 'La', 'Ce': 'Ce', 'Pr': 'Pr_3', 'Nd': 'Nd_3', 'Pm': 'Pm_3',
        'Sm': 'Sm_3', 'Eu': 'Eu_2', 'Gd': 'Gd_3', 'Tb': 'Tb_3', 'Dy': 'Dy_3', 'Ho': 'Ho_3', 'Er': 'Er_3', 'Tm': 'Tm_3',
        'Yb': 'Yb_2', 'Lu': 'Lu_3', 'Hf': 'Hf_pv', 'Ta': 'Ta_pv', 'W': 'W_sv', 'Re': 'Re', 'Os': 'Os', 'Ir': 'Ir',
        'Pt': 'Pt', 'Au': 'Au', 'Hg': 'Hg', 'Tl': 'Tl_d', 'Pb': 'Pb_d', 'Bi': 'Bi_d', 'Po': 'Po_d', 'At': 'At',
        'Rn': 'Rn', 'Fr': 'Fr_sv', 'Ra': 'Ra_sv', 'Ac': 'Ac', 'Th': 'Th', 'Pa': 'Pa', 'U': 'U', 'Np': 'Np', 'Pu': 'Pu',
        'Am': 'Am', 'Cm': 'Cm'
    }

    def get_paw(self, symbol):
        if symbol in self.__paws.keys():
            return self.__paws[symbol]
        else:
            raise KeyError('Such PAW_potential doesn\'t exist, try to set it manually by self.set_paw(symbol, paw)')

    def set_paw(self, symbol, paw):
        '''
        set your non-default paw
        :param symbol: Chemical element
        :param paw: Name of needed paw (Na_pv, V_sv etc.)
        :return:
        '''
        self.__paws[symbol] = paw

    def get_valence(self, structure: Structure):
        potcar = Potcar([self.get_paw(el.symbol) for el in structure.elements], functional='PBE_64')
        valences = []
        for element in potcar:
            valence = 0.0
            for level in element.get_electron_configuration():
                valence += level[2]
            valences.append(valence)
        return list(map(int, valences))

    def get_lattice(self, structure: Structure):
        template = '''LATTICE
        SYSDIM       3D
        SYSTYPE      BULK
        BRAVAIS      {bravais_code}
        ALAT         {alat}
        A(1)         {a}
        A(2)         {b}
        A(3)         {c}'''
        structure_prim = SgA(structure).get_primitive_standard_structure()
        primitive_vectros = structure_prim.lattice.matrix / structure.lattice.a
        places = ['{bravais_code}', '{alat}', '{a}', '{b}', '{c}']
        replaces = [f"{' '.join([str(i) for i in Pearson.from_symbol(SgA(structure).get_pearson_symbol()[:2])][1:])}",
                    f'{structure.lattice.a / 0.529177:.4f}',
                    ''.join([f'{i:>10.4f}' for i in primitive_vectros[0]]),
                    ''.join([f'{i:>10.4f}' for i in primitive_vectros[1]]),
                    ''.join([f'{i:>10.4f}' for i in primitive_vectros[2]])
                    ]
        for src, dst in zip(places, replaces):
            template = template.replace(src, dst)
        # print(template)
        return template

    def get_sites(self, structure_prim: Structure):
        filestring = f'\t\tIQ{" "*10}QX{" "*10}QY{" "*10}QZ\n'
        structure = SgA(structure_prim).get_conventional_standard_structure()
        for i, site in enumerate(structure_prim.sites):
            pos = site.coords / structure.lattice.a
            filestring += f"\t\t{i + 1}\t\t{pos[0]:.4f}\t\t{pos[1]:.4f}\t\t{pos[2]:.4f}\n"
        return filestring[:-1]

    def get_occupation(self, structure_prim: Structure):
        filestring = '\t\tIQ\t\tIREFQ\t\tIMQ\t\tNOQ\t\tITOQ\t\tCONC\n'
        occupation = {'IQ': [],     # Индекс атома
                      'IREFQ': [],  # Какому референс-потенциалу принадлежит
                      'IMQ': [],    # Магнитная подрешетка
                      'NOQ': [],    # Число состояний в данной позиции
                      'ITOQ': [],   # Тип атома (ссылается на блок TYPES)
                      'CONC': []}   # Концентрация (1.0 — чистый элемент)

        unique_elements = [el.symbol for el in structure_prim.composition.elements]

        element_to_it = {el: i + 1 for i, el in enumerate(unique_elements)}
        for i, site in enumerate(structure_prim.sites):
            occupation['IQ'].append(i+1)
            occupation['IREFQ'].append(element_to_it[site.specie.symbol])
            occupation['IMQ'].append(element_to_it[site.specie.symbol])
            occupation['NOQ'].append(1)
            occupation['ITOQ'].append(element_to_it[site.specie.symbol])
            occupation['CONC'].append(1.0)
            filestring += '\t\t' + '\t\t'.join([str(val[i]) for val in occupation.values()])
            filestring += '\n'
        return filestring[:-1]

    def get_ref_system(self, structure_prim: Structure):
        reference_system = {'IREF': [],
                            'VREF': [],  # Усредненный потенциал (для старта)
                            'RMTREF': []}  # Радиус сферы МТ (может быть 0 на старте)
        filestring = f'NREF\t{len(structure_prim.elements)}\n'
        filestring += '\t\tIREF\t\tVREF\t\tRMTREF\n'
        unique_elements = [el.symbol for el in structure_prim.composition.elements]
        for i, _ in enumerate(unique_elements):
            filestring += f'\t\t{i+1}\t\t4.0\t\t0.0\n'
        return filestring[:-1]

    def get_mag_direction(self, structure_prim: Structure):
        mag_direction = {'IQ': [],
                         'QMTET': [],  # Угол θ (в радианах или градусах — зависит от настроек)
                         'QMPHI': []}  # Угол φ (азимут)
        filestring = '\t\tIQ\t\tQMTET\t\tQMPHI\n'
        for i, _ in enumerate(structure_prim.sites):
            filestring += f'\t\t{i+1}\t\t0.0\t\t0.0\n'
        return filestring[:-1]

    def get_mesh_info(self, structure_prim: Structure):
        structure = SgA(structure_prim).get_conventional_standard_structure()
        filestring = '\t\tIM\t\tR(1)\t\tDX\t\tJRMT\t\tRMT\t\tJRWS\t\tRWS\n'
        mesh_info = {'IM': [],  # Индекс сетки (по атомам)
                     'R(1)': [],  # Начальный радиус сетки
                     'DX': [],  # Шаг сетки
                     'JRMT': [],  # Индекс точки, соответствующей RMT
                     'RMT': [],  # Радиус МТ-сферы
                     'JRWS': [],  # Индекс для границы Wigner-Seitz
                     'RWS': []}  # Радиус WS-сферы
        base_radii = {element.symbol: element.atomic_radius_calculated for element in structure.composition}
        rws = get_rws(structure, base_radii=base_radii)
        for i, el in enumerate(structure_prim.composition.elements):
            mesh_info['IM'].append(i+1)
            mesh_info['R(1)'].append(0.000001)
            mesh_info['DX'].append(f'{np.log(rws[el.symbol]/0.000001)/720:.12f}')
            mesh_info['JRMT'].append(0)
            mesh_info['RMT'].append(f'{rws[el.symbol] * 0.85:.6f}')
            mesh_info['JRWS'].append(721)
            mesh_info['RWS'].append(f'{rws[el.symbol]:.6f}')
            filestring += '\t\t' + '\t\t'.join([str(val[i]) for val in mesh_info.values()])
            filestring += '\n'
        return filestring[:-1]


    def get_types(self, structure_prim: Structure):
        types = {'IT': [],  # Номер типа атома
                 'TXTT': [],  # Химическое обозначение
                 'ZT': [],  # Заряд ядра
                 'NCORT': [],  # Кол-во электронов в core
                 'NVALT': [],  # Кол-во валентных
                 'NSEMCORSHLT': []}  # Полусердечные состояния (0 — игнорируются)
        filestring = '\t\tIT\t\tTXTT\t\tZT\t\tNCORT\t\tNVALT\t\tNSEMCORSHLT\n'
        valences = get_valence(structure_prim)
        for i, el in enumerate(structure_prim.composition.elements):
            types['IT'].append(i + 1)
            types['TXTT'].append(el.symbol)
            types['ZT'].append(el.Z)
            types['NCORT'].append(el.Z - valences[i])
            types['NVALT'].append(valences[i])
            types['NSEMCORSHLT'].append(0)
            filestring += '\t\t' + '\t\t'.join([str(val[i]) for val in types.values()])
            filestring += '\n'
        return filestring[:-1]


    def generate_pot(self, structure: Structure):
        structure = SgA(structure).get_conventional_standard_structure()
        prim_s = SgA(structure).get_primitive_standard_structure()

        replaces = {'{generation_date}': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    '{system}': prim_s.composition.to_pretty_string().replace('1', ''),
                    '{nq}': f'{prim_s.num_sites}',  # Кол-во неэквивалентных позиций в кристалле
                    '{nt}': f'{len(structure.elements)}',  # Кол-во уникальных типов атомов
                    '{nm}': f'{len(structure.elements)}',  # Кол-во магнитных подрешёток (или групп атомов с разными магнитными свойствами)
                    '{irel}': '2',
                    '{ne}': str(sum(get_valence(prim_s))),
                    '{ibzint}': '2',
                    '{nktab}': '4000',
                    '{lattice}': get_lattice(structure),
                    '{sites}': get_sites(prim_s),
                    '{occupations}': get_occupation(prim_s),
                    '{ref_system}': get_ref_system(prim_s),
                    '{mag_directions}': get_mag_direction(prim_s),
                    '{mesh_info}': get_mesh_info(prim_s),
                    '{types}': get_types(prim_s)}

        with open('template_pot.txt', 'r') as f:
            template = f.read()
        for key, val in replaces.items():
            template = template.replace(key, val)
        return template


    def write_pot(self, filestring: str, output_path: Path):
        with open(output_path, 'w') as f:
            f.write(filestring)


if __name__ == '__main__':
    pass

