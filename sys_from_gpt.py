"""
sprkkr_sysgen.py

Генератор .sys (xband-style) на основе данных из generate_sys_data().
"""

from pathlib import Path
import datetime
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from typing import Dict, List, Tuple

from sys_data_gpt import generate_sys_data  # импортируем общую функцию
import numpy as np


# ============================================================
# Форматирование
# ============================================================
def _fmt(x: float, width: float = 18, prec: int = 12) -> str:
    return f"{x:{width}.{prec}f}"


def _fmt_short(x: float, width: int = 12, prec: int = 6) -> str:
    return f"{x:{width}.{prec}f}"


# ============================================================
# Формирование разделов SYS
# ============================================================
def header_section(system_name: str, xband_version: str = "6.3") -> str:
    now = datetime.datetime.now().strftime("%a %d %b %H:%M:%S %Z %Y")
    return (
        f"system data-file created by sprkkr_sysgen on {now}\n"
        f"{system_name}\n"
        f"xband-version\n{xband_version}\n"
        f"dimension\n3D\n"
    )


def bravais_section(spacegroup_number: Tuple[int], bravais: str) -> str:
    """Формирует описание браavais и номер группы симметрии."""
    return (
        f"Bravais lattice\n"
        f"{bravais}\n"
        f"space group number (ITXC and AP)\n"
        f"{spacegroup_number[0]:4d}  {spacegroup_number[1]}\n"
        f"structure type\nUNKNOWN\n"
    )


def lattice_section(data: Dict) -> str:
    """Блок lattice, формируется напрямую из данных generate_sys_data."""
    prim_matrix = data["prim_matrix"]
    a_au = data["a_au"]
    a_ang = data["a_ang"]
    prim = data["prim_structure"]
    conv = data["conv_structure"]
    c_over_a = conv.lattice.c / conv.lattice.a
    angles = conv.lattice.angles

    s = "lattice parameter A  [a.u.]\n"
    s += f"    {_fmt(a_au)}\n"
    s += "ratio of lattice parameters  b/a  c/a\n"
    s += f"    {_fmt(1.0)}    {_fmt(c_over_a)}\n"
    s += "lattice parameters  a b c  [a.u.]\n"
    s += f"    {_fmt(a_au)}    {_fmt(a_au)}    {_fmt(a_au * c_over_a)}\n"
    s += "lattice angles  alpha beta gamma  [deg]\n"
    s += f"   {_fmt(angles[0])}   {_fmt(angles[1])}   {_fmt(angles[2])}\n"
    s += "primitive vectors     (cart. coord.) [A]\n"
    for v in prim_matrix:
        s += f"    {_fmt(v[0])}   {_fmt(v[1])}   {_fmt(v[2])}\n"
    return s


def sites_section(data: Dict) -> str:
    """Создаёт секцию NQ (atomic sites)."""
    prim = data["prim_structure"]
    rws_dict = data["rws_dict"]
    a_ang = data['a_ang']

    n_sites = len(prim.sites)
    s = f"number of sites NQ\n  {n_sites}\n"
    s += " IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"
    for i, site in enumerate(prim.sites, start=1):
        frac = site.coords / a_ang
        elem = site.specie.symbol
        rws = rws_dict.get(elem, 2.5)
        s += (
            f"{i:3d}   {site.type_idx:1d}   {_fmt_short(frac[0])}   {_fmt_short(frac[1])}   {_fmt_short(frac[2])}"
            f"       {_fmt(rws, width=12, prec=12)}   3    1   {site.type_idx}\n"
        )

    s += "number of sites classes NCL\n"
    s += f"  {n_sites}\n"
    s += "ICL WYCK NQCL IQECL (equivalent sites)\n"

    symmetrized = data["symmetrized"].equivalent_sites
    symmetrized_idx = data["symmetrized"].equivalent_indices
    for i, _ in enumerate(symmetrized, start=1):
        sites_str = " ".join([str(j + 1) for j in symmetrized_idx[i - 1]])
        s += f"  {i:1d}   -    {len(sites_str.split())}  {sites_str}\n"
    return s


def atom_types_section(data: Dict) -> str:
    """Секция ATOM TYPES."""
    symmetrized = data["symmetrized"].equivalent_sites
    symmetrized_idx = data["symmetrized"].equivalent_indices
    s = "number of atom types NT\n"
    s += f"  {len(symmetrized_idx)}\n"
    s += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"

    for i, equiv in enumerate(symmetrized, start=1):
        label = equiv[0].type
        Z = equiv[0].specie.Z
        sites_str = " ".join([str(j + 1) for j in symmetrized_idx[i - 1]])
        s += f"  {i:2d}  {Z:3d}  {label:<8} {len(equiv):3d} 1.000  {sites_str}\n"

    return s


# ============================================================
# Основная функция
# ============================================================
def generate_sys_from_data(data: Dict) -> str:
    """Генерирует полный текст .sys-файла на основе словаря из generate_sys_data."""
    sys_txt = ""
    sys_txt += header_section(data["system_name"])
    sys_txt += bravais_section(data["spacegroup"], data["bravais"])
    sys_txt += lattice_section(data)
    sys_txt += sites_section(data)
    sys_txt += atom_types_section(data)
    return sys_txt


# ============================================================
# Пример использования
# ============================================================
if __name__ == "__main__":
    from pymatgen.io.vasp import Poscar, Outcar


    path = Path("SPR_KKR/Al/Tsharp/CONTCAR")
    pos = Poscar.from_file(path)
    out = Outcar(path.parent / 'OUTCAR')
    magmoms = [mag['tot'] for mag in out.magnetization]
    structure = pos.structure

    data = generate_sys_data(structure)
    print(generate_sys_from_data(data))

    exit()
    # for i, _ in enumerate(structure.sites):
    #     structure.sites[i].spin = magmoms[i]

    sga = SpacegroupAnalyzer(structure)
    prim = sga.get_primitive_standard_structure(keep_site_properties=True)
    for d in prim.sites[0].__dict__:
        print(d)

    print(prim.sites[0].__dict__)

    if hasattr(prim.sites[0], 'spin'):
        print('Spin saved')
        print(prim.sites[0].spin)
    exit()


    data = generate_sys_data(structure, magmoms)
    sys_text = generate_sys_from_data(data)

    print(sys_text)
    # with open(f"{data['system_name']}.sys", "w") as f:
    #     f.write(sys_text)



