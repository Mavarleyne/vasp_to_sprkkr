"""
sprkkr_potgen.py

Генератор .pot-файлов для SPR-KKR/xband.
Использует общие данные из sprkkr_sysgen.generate_sys_data().
"""

import datetime
import math
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar, Potcar

from sys_from_gpt import generate_sys_data


paws = {
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


def get_valence(structure: Structure):
    potcar = Potcar([paws[el.symbol] for el in structure.elements], functional='PBE_64')
    valences = []
    for element in potcar:
        valence = 0.0
        for level in element.get_electron_configuration():
            valence += level[2]
        valences.append(valence)
    # print({el.symbol.split('_')[0]: val for el, val in zip(potcar, valences)})
    return {el.symbol.split('_')[0]: val for el, val in zip(potcar, valences)}


def generate_pot_file(
    data: dict,
    xc_potential="PBE",
    relativity=2,
    mesh_type="EXPONENTIAL",
    scf_tol=1e-5,
    scf_mix=0.2,
):

    bravais = data["bravais"]
    system_name = data["system_name"]
    prim = data["prim_structure"]
    atom_types = data["atom_types"]
    site2type = data["site2type"]
    rws_dict = data["rws_dict"]
    prim_matrix = data["prim_matrix"]
    a_au = data["a_au"]
    a_ang = data['a_ang']

    valence = get_valence(prim)

    now = datetime.datetime.now().strftime("%a %d %b %H:%M:%S %Z %Y")

    NQ = len(prim.sites)
    NT = len(atom_types)
    NM = NT

    filestring = []
    filestring.append("*******************************************************************************\n")
    filestring.append(f"HEADER    'SCF-start data created by sprkkr_potgen  {now}'\n")
    filestring.append("*******************************************************************************\n")
    filestring.append(f"TITLE     'SPR-KKR calculation for {system_name}'\n")
    filestring.append(f"SYSTEM    {system_name}\n")
    filestring.append("PACKAGE   SPRKKR\n")
    filestring.append("FORMAT    6  (21.05.2007)\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("GLOBAL SYSTEM PARAMETER\n")
    filestring.append(f"NQ{NQ:>19}\n")
    filestring.append(f"NT{NT:>19}\n")
    filestring.append(f"NM{NM:>19}\n")
    filestring.append(f"IREL{relativity:>17}\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("SCF-INFO\n")
    filestring.append("INFO      NONE\n")
    filestring.append("SCFSTATUS START\n")
    filestring.append("FULLPOT   F\n")
    filestring.append("BREITINT  F\n")
    filestring.append("NONMAG    F\n")
    filestring.append("ORBPOL    NONE\n")
    filestring.append("EXTFIELD  F\n")
    filestring.append("BLCOUPL   F\n")
    filestring.append("BEXT          0.0000000000\n")
    filestring.append("SEMICORE  F\n")
    filestring.append("LLOYD     F\n")
    filestring.append("NE               30\n")
    filestring.append("IBZINT            2\n")
    filestring.append("NKTAB             0\n")
    filestring.append(f"XC-POT    {xc_potential}\n")
    filestring.append("SCF-ALG   BROYDEN2\n")
    filestring.append("SCF-ITER           0\n")
    filestring.append(f"SCF-MIX       {scf_mix:.10}\n")
    filestring.append(f"SCF-TOL       {scf_tol:.10}\n")
    filestring.append("RMSAVV    999999.0000000000\n")
    filestring.append("RMSAVB    999999.0000000000\n")
    filestring.append("EF            0.0000000000\n")
    filestring.append("VMTZ          0.0000000000\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("LATTICE\n")
    filestring.append("SYSDIM       3D\n")
    filestring.append("SYSTYPE      BULK\n")
    filestring.append(f"BRAVAIS{bravais}\n")
    filestring.append(f"ALAT{a_au:>22.10}\n")
    for i, v in enumerate(prim_matrix, start=1):
        filestring.append(f"A({i}){v[0]:>22.10f}    {v[1]:>16.10f}    {v[2]:>16.10f}\n")

    filestring.append("*******************************************************************************\n")

    filestring.append("SITES\n")
    filestring.append("CARTESIAN T\n")
    filestring.append("BASSCALE      1.0000000000    1.0000000000    1.0000000000\n")
    filestring.append("        IQ      QX              QY              QZ\n")
    for i, site in enumerate(prim.sites, start=1):
        x, y, z = site.coords / a_ang
        filestring.append(f"{i:>10d}{x:>16.10f}{y:>16.10f}{z:>16.10f}\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("OCCUPATION\n")
    filestring.append("        IQ     IREFQ       IMQ       NOQ  ITOQ  CONC\n")
    for i, tindex in enumerate(site2type, start=1):
        filestring.append(f"{i:10d}{tindex:10d}{tindex:10d}{1:10d}{tindex:6d} 1.000\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("REFERENCE SYSTEM\n")
    filestring.append(f"NREF{NT:>16}\n")
    filestring.append("      IREF      VREF            RMTREF\n")
    for i in range(1, NT + 1):
        filestring.append(f"{i:10d}    4.0000000000    0.0000000000\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("MAGNETISATION DIRECTION\n")
    filestring.append("KMROT              0\n")
    filestring.append("QMVEC         0.0000000000    0.0000000000    0.0000000000\n")
    filestring.append("        IQ      QMTET           QMPHI \n")
    for i in range(1, NQ + 1):
        filestring.append(f"{i:10d}    0.0000000000    0.0000000000\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("MESH INFORMATION\n")
    filestring.append(f"MESH-TYPE {mesh_type} \n")
    filestring.append("   IM      R(1)            DX         JRMT      RMT        JRWS      RWS\n")
    for i, (label, Z, avg_mag, inds) in enumerate(atom_types, start=1):
        rws = rws_dict.get(label.split("_")[0], 2.6)
        R1 = 1e-6
        JRWS = 721
        DX = math.log(rws / R1) / (JRWS - 1)
        RMT = rws * 0.85
        filestring.append(f"{i:5d}    {R1:.10f}    {DX:.10f}    0   {RMT: .10f}  {JRWS}   {rws: .10f}\n")
    filestring.append("*******************************************************************************\n")

    filestring.append("TYPES\n")
    filestring.append("   IT     TXTT        ZT     NCORT     NVALT    NSEMCORSHLT\n")
    for i, (label, Z, avg_mag, inds) in enumerate(atom_types, start=1):
        # val = {26: 8, 29: 11, 22: 4, 25: 7, 24: 6, 78: 10}.get(Z, 8)
        filestring.append(f"{i:5d}     {label:<4s}{Z:>14d}{int(Z - valence[label.split('_')[0]]):10d}{int(valence[label.split('_')[0]]):10d}{0:15d}\n")
    filestring.append("*******************************************************************************\n")

    return ''.join(filestring)
    # print(f"✅ POT-файл успешно записан: {output_path}")


if __name__ == '__main__':
    # pos = Poscar.from_file(Path("for_spr/Ti4Fe8Cu4/119/FiM/POSCAR"))
    # magmoms = [-0.534, 1.755, 2.137, -0.035]  # при необходимости
    pos = Poscar.from_file(Path("for_spr/Mn8Cr4Pt4/139/FiM/POSCAR"))
    magmoms = [3.178, 3.178, -2.593, 0.157]
    structure = pos.structure
    pot_text = generate_pot_file(structure, magmoms)

    # print(sys_text)
    print('#' * 100)
    print(pot_text)
