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


def bravais_section(spacegroup_number: int, bravais: str) -> str:
    """Формирует описание браavais и номер группы симметрии."""
    return (
        f"Bravais lattice\n"
        f"{bravais}\n"
        f"space group number (ITXC and AP)\n"
        f"{spacegroup_number:4d}  361\n"
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
    site2type = data["site2type"]
    a_ang = data['a_ang']

    n_sites = len(prim.sites)
    s = f"number of sites NQ\n  {n_sites}\n"
    s += " IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"
    for i, site in enumerate(prim.sites, start=1):
        # frac = np.mod(site.frac_coords, 1.0)
        # frac = site.frac_coords
        frac = site.coords / a_ang
        elem = site.specie.symbol
        rws = rws_dict.get(elem, 2.5)
        s += (
            f"{i:3d}   {i:1d}   {_fmt_short(frac[0])}   {_fmt_short(frac[1])}   {_fmt_short(frac[2])}"
            f"       {_fmt(rws, width=12, prec=12)}   3    1   {site2type[i-1]}\n"
        )

    s += "number of sites classes NCL\n"
    s += f"  {n_sites}\n"
    s += "ICL WYCK NQCL IQECL (equivalent sites)\n"
    # for i in range(1, n_sites + 1):
    #     s += f"  {i:1d}   -    1  {i}\n"

    atom_types = data["atom_types"]
    for i, (label, Z, avg_mag, site_inds) in enumerate(atom_types, start=1):
        sites_str = " ".join(str(idx + 1) for idx in site_inds)
        s += f"  {i:1d}   -    {len(sites_str.split())}  {sites_str}\n"
    return s


def atom_types_section(data: Dict) -> str:
    """Секция ATOM TYPES."""
    atom_types = data["atom_types"]
    s = "number of atom types NT\n"
    s += f"  {len(atom_types)}\n"
    s += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
    for i, (label, Z, avg_mag, site_inds) in enumerate(atom_types, start=1):
        sites_str = " ".join(str(idx + 1) for idx in site_inds)
        s += f"  {i:2d}  {Z:3d}  {label:<8} {len(site_inds):3d} 1.000  {sites_str}\n"
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
    from pymatgen.io.vasp import Poscar

    pos = Poscar.from_file(Path("for_spr/Ti4Fe8Cu4/119/FiM/POSCAR"))
    # pos = Poscar.from_file(Path("for_spr/Mn8Cr4Pt4/139/FiM/POSCAR"))
    structure = pos.structure
    magmoms = [-0.534, 1.755, 2.137, -0.035]  # при необходимости
    # magmoms = [3.178, 3.178, -2.593, 0.157]

    # pos = Poscar.from_file(Path("for_spr/Mn8Cr4Pt4/139/FiM/POSCAR"))
    # structure = pos.structure
    # magmoms = [3.178, 3.178, -2.593, 0.157]

    data = generate_sys_data(structure, magmoms)
    sys_text = generate_sys_from_data(data)

    print(sys_text)
    # with open(f"{data['system_name']}.sys", "w") as f:
    #     f.write(sys_text)



