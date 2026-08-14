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

# from Wigner_Seitz_radius import get_rws_physical
from System import generate_sys_data
from pymatgen.io.vasp.inputs import PotcarSingle
PotcarSingle.functional_dir['PBE_64'] = "/home/buche/VaspPot_64/potpaw_PBE.64"



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
    potcars_dir = Path("/home/buche/VaspPot_64/potpaw_PBE.64")

    # Множество доступных имён потенциалов (имена подпапок = символы POTCAR)
    available = {p.name for p in potcars_dir.iterdir() if p.is_dir()}

    potcars = []

    for el in structure.composition.elements:
        symbol = el.symbol

        candidates = [f'{symbol}_pv']
        if symbol in paws:
            candidates.append(paws[symbol])
        candidates.append(f'{symbol}_sv')
        candidates.append(symbol)

        chosen = None
        for candidate in candidates:
            if candidate in available:
                chosen = candidate
                break

        if chosen is None:
            raise OSError(
                f"Не найден POTCAR для элемента '{symbol}' в {potcars_dir}. "
                f"Проверенные варианты: {candidates}"
            )

        potcars.append(chosen)

    potcar = Potcar(potcars, functional='PBE_64')
    valences = []
    for element in potcar:
        # print(dir(element))
        # print(element.electron_configuration)
        # exit()
        valence = 0.0
        for level in element.electron_configuration:
            # print()
            # exit()
            valence += level[2]
        valences.append(valence)
    return {el.symbol.split('_')[0]: val for el, val in zip(potcar, valences)}


def pot_header_section(system_name: str) -> str:
    now = datetime.datetime.now().strftime("%a %d %b %H:%M:%S %Z %Y")

    s = "*******************************************************************************\n"
    s += f"HEADER    'SCF-start data created by sprkkr_potgen  {now}'\n"
    s += "*******************************************************************************\n"
    s += f"TITLE     'SPR-KKR calculation for {system_name}'\n"
    s += f"SYSTEM    {system_name}\n"
    s += "PACKAGE   SPRKKR\n"
    s += "FORMAT    6  (21.05.2007)\n"
    s += "*******************************************************************************\n"
    return s


def pot_global_section(data, relativity: int) -> str:
    prim = data["prim_structure"]
    # atom_types = data["atom_types"]
    symmetrized = data['symmetrized']

    NQ = len(prim.sites)
    NT = len(symmetrized.equivalent_indices)
    NM = NT

    s = "GLOBAL SYSTEM PARAMETER\n"
    s += f"NQ{NQ:>19}\n"
    s += f"NT{NT:>19}\n"
    s += f"NM{NM:>19}\n"
    s += f"IREL{relativity:>17}\n"
    s += "*******************************************************************************\n"
    return s


def pot_scf_section(xc_potential, scf_mix, scf_tol) -> str:
    s = "SCF-INFO\n"
    s += "INFO      NONE\n"
    s += "SCFSTATUS START\n"
    s += "FULLPOT   F\n"
    s += "BREITINT  F\n"
    s += "NONMAG    F\n"
    s += "ORBPOL    NONE\n"
    s += "EXTFIELD  F\n"
    s += "BLCOUPL   F\n"
    s += "BEXT          0.0000000000\n"
    s += "SEMICORE  F\n"
    s += "LLOYD     F\n"
    s += "NE               30\n"
    s += "IBZINT            2\n"
    s += "NKTAB             0\n"
    s += f"XC-POT    {xc_potential}\n"
    s += "SCF-ALG   BROYDEN2\n"
    s += "SCF-ITER           0\n"
    s += f"SCF-MIX       {scf_mix:.10}\n"
    s += f"SCF-TOL       {scf_tol:.10}\n"
    s += "RMSAVV    999999.0000000000\n"
    s += "RMSAVB    999999.0000000000\n"
    s += "EF            0.0000000000\n"
    s += "VMTZ          0.0000000000\n"
    s += "*******************************************************************************\n"
    return s


def pot_lattice_section(data) -> str:
    prim_matrix = data["prim_matrix"]
    a_au = data["a_au"]
    bravais = data["bravais"]

    s = "LATTICE\n"
    s += "SYSDIM       3D\n"
    s += "SYSTYPE      BULK\n"
    s += f"BRAVAIS{bravais}\n"
    s += f"ALAT{a_au:>22.10}\n"

    for i, v in enumerate(prim_matrix, start=1):
        s += f"A({i}){v[0]:>22.10f}    {v[1]:>16.10f}    {v[2]:>16.10f}\n"

    s += "*******************************************************************************\n"
    return s


def pot_sites_section(data) -> str:
    prim = data["prim_structure"]
    a_ang = data["a_ang"]

    s = "SITES\n"
    s += "CARTESIAN T\n"
    s += "BASSCALE      1.0000000000    1.0000000000    1.0000000000\n"
    s += "        IQ      QX              QY              QZ\n"

    for i, site in enumerate(prim.sites, start=1):
        x, y, z = site.coords / a_ang
        s += f"{i:>10d}{x:>16.10f}{y:>16.10f}{z:>16.10f}\n"

    s += "*******************************************************************************\n"
    return s


def pot_occupation_section(data) -> str:
    symmetrized = data["symmetrized"]
    prim = data["prim_structure"]
    s = "OCCUPATION\n"
    s += "        IQ     IREFQ       IMQ       NOQ  ITOQ  CONC\n"

    for IQ, site in enumerate(prim.sites, start=1):
        s += f"{IQ:10d}{site.type_idx:10d}{site.type_idx:10d}{1:10d}{site.type_idx:6d} 1.000\n"

    s += "*******************************************************************************\n"
    return s


def pot_reference_section(data) -> str:
    NT = len(data["symmetrized"].equivalent_indices)

    s = "REFERENCE SYSTEM\n"
    s += f"NREF{NT:>16}\n"
    s += "      IREF      VREF            RMTREF\n"

    for i in range(1, NT + 1):
        s += f"{i:10d}    4.0000000000    0.0000000000\n"

    s += "*******************************************************************************\n"
    return s


def pot_magnetisation_section(data) -> str:
    NQ = len(data["prim_structure"].sites)

    s = "MAGNETISATION DIRECTION\n"
    s += "KMROT              0\n"
    s += "QMVEC         0.0000000000    0.0000000000    0.0000000000\n"
    s += "        IQ      QMTET           QMPHI \n"

    for i in range(1, NQ + 1):
        s += f"{i:10d}    0.0000000000    0.0000000000\n"

    s += "*******************************************************************************\n"
    return s


def pot_mesh_section(data, mesh_type="EXPONENTIAL") -> str:
    symmetrized = data["symmetrized"].equivalent_sites
    # symmetrized_idx = data["symmetrized"].equivalent_indices
    rws_dict = data["rws_dict"]

    s = "MESH INFORMATION\n"
    s += f"MESH-TYPE {mesh_type} \n"
    s += "   IM      R(1)            DX         JRMT      RMT        JRWS      RWS\n"

    for i, equiv_set in enumerate(symmetrized, start=1):
        rws = rws_dict.get(equiv_set[0].specie.symbol, 2.6)
        R1 = 1e-6
        JRWS = 721
        DX = math.log(rws / R1) / (JRWS - 1)
        RMT = rws * 0.85
        if hasattr(equiv_set[0], 'type_idx'):
            IM = equiv_set[0].type_idx
        else:
            raise AttributeError('Site doesn\'t have attribute "type_idx"')
        s += f"{IM:5d}    {R1:.10f}    {DX:.10f}    0   {RMT: .10f}  {JRWS}   {rws: .10f}\n"

    s += "*******************************************************************************\n"
    return s


def pot_types_section(data, valence) -> str:
    # symmetrized = data["symmetrized"].equivalent_sites
    symmetrized_idx = data["symmetrized"].equivalent_indices
    prim = data['prim_structure']
    sites = prim.sites
    s = "TYPES\n"
    s += "   IT     TXTT        ZT     NCORT     NVALT    NSEMCORSHLT\n"

    for i, sym_sites in enumerate(symmetrized_idx, start=1):
        type = sites[sym_sites[0]].type
        base = sites[sym_sites[0]].specie.symbol
        nval = int(valence[base])
        Z = sites[sym_sites[0]].specie.Z
        ncore = int(Z - nval)

        s += f"{i:5d}     {type:<4s}{Z:>14d}{ncore:10d}{nval:10d}{0:15d}\n"
    s += "*******************************************************************************\n"
    return s


def generate_pot_from_data(
    data: dict,
    xc_potential="PBE",
    relativity=2,
    mesh_type="EXPONENTIAL",
    scf_tol=1e-5,
    scf_mix=0.2,
):

    valence = get_valence(data["prim_structure"])

    pot = ""
    pot += pot_header_section(data["system_name"])
    pot += pot_global_section(data, relativity)
    pot += pot_scf_section(xc_potential, scf_mix, scf_tol)
    pot += pot_lattice_section(data)
    pot += pot_sites_section(data)
    pot += pot_occupation_section(data)
    pot += pot_reference_section(data)
    pot += pot_magnetisation_section(data)
    pot += pot_mesh_section(data, mesh_type)
    # pot += pot_mesh_section_physical(data)
    pot += pot_types_section(data, valence)

    return pot


if __name__ == '__main__':
    # pos = Poscar.from_file(Path("for_spr/Ti4Fe8Cu4/119/FiM/POSCAR"))
    # magmoms = [-0.534, 1.755, 2.137, -0.035]  # при необходимости
    pos = Poscar.from_file(Path("SPR_KKR/Al/Tsharp/CONTCAR"))
    magmoms = [3.178, 3.178, -2.593, 0.157]

    structure = pos.structure
    data = generate_sys_data(structure)
    # pot_text = generate_pot_file(structure, magmoms)
    pot_text = generate_pot_from_data(data)

    # print(sys_text)
    print('#' * 100)
    print(pot_text)
