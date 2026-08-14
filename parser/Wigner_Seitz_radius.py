import os
import numpy as np
from pathlib import Path
from pymatgen.core import Structure


base_radii = {
    "H":  0.53,   "He": 0.31,
    "Li": 1.67,   "Be": 1.12, "B":  0.87,   "C":  0.67,   "N":  0.56,
    "O":  0.48,   "F":  0.42, "Ne": 0.38,
    "Na": 1.90,   "Mg": 1.45, "Al": 1.18,  "Si": 1.11,  "P":  0.98,
    "S":  0.87,   "Cl": 0.79, "Ar": 0.71,
    "K":  2.43,   "Ca": 1.94, "Sc": 1.84,  "Ti": 1.76,  "V":  1.71,
    "Cr": 1.66,   "Mn": 1.61, "Fe": 1.56,  "Co": 1.52,  "Ni": 1.49,
    "Cu": 1.45,   "Zn": 1.42, "Ga": 1.36,  "Ge": 1.25,  "As": 1.14,
    "Se": 1.03,   "Br": 0.94, "Kr": 0.87,
    "Rb": 2.65,   "Sr": 2.19, "Y":  2.12,  "Zr": 2.06,  "Nb": 1.98,
    "Mo": 1.90,   "Tc": 1.83, "Ru": 1.78,  "Rh": 1.73,  "Pd": 1.69,
    "Ag": 1.65,   "Cd": 1.61, "In": 1.56,  "Sn": 1.45,  "Sb": 1.33,
    "Te": 1.23,   "I":  1.15, "Xe": 1.08,
    "Cs": 2.98,   "Ba": 2.53,

    # Лантаноиды
    "La": None,   "Ce": None, "Pr": 2.47,  "Nd": 2.06,  "Pm": 2.05,
    "Sm": 2.38,   "Eu": 2.31, "Gd": 2.33,  "Tb": 2.25,  "Dy": 2.28,
    "Ho": 2.26,   "Er": 2.26, "Tm": 2.22,  "Yb": 2.22,  "Lu": 2.17,

    # 6-й период (после лантаноидов)
    "Hf": 2.08,   "Ta": 2.00, "W":  1.93,  "Re": 1.88,  "Os": 1.85,
    "Ir": 1.80,   "Pt": 1.77, "Au": 1.74,  "Hg": 1.71,  "Tl": 1.56,
    "Pb": 1.54,   "Bi": 1.43, "Po": 1.35,  "At": 1.27,  "Rn": 1.20,

    # Актиноиды и сверхтяжёлые элементы
    "Fr": None,   "Ra": None, "Ac": None,
    "Th": None,   "Pa": None, "U":  None, "Np": None, "Pu": None,
    "Am": None,   "Cm": None, "Bk": None, "Cf": None, "Es": None,
    "Fm": None,   "Md": None, "No": None, "Lr": None,
    "Rf": None,   "Db": None, "Sg": None, "Bh": None, "Hs": None,
    "Mt": None,   "Ds": None, "Rg": None, "Cn": None, "Nh": None,
    "Fl": None,   "Mc": None, "Lv": None, "Ts": None, "Og": None,
}


def load_rws_vst(rws_file_path="rws.vst"):
    """
    Загружает таблицу базовых радиусов Вигнера-Зейтца из файла rws.vst.

    Args:
        rws_file_path (str): Путь к файлу rws.vst

    Returns:
        dict: Словарь {элемент: радиус в Å}
    """
    rws_vst = {}
    try:
        with open(Path(rws_file_path), 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        element, radius = parts[-1], float(parts[1])
                        rws_vst[element] = radius
    except FileNotFoundError:
        print(f"Файл {rws_file_path} не найден. Используются дефолтные радиусы.")
        # Дефолтные радиусы для Mn, Cr, Pt (на основе целевых значений)
        # rws_vst = {'Mn': 1.35, 'Cr': 1.34, 'Pt': 1.45}
    # print(rws_vst)
    return rws_vst


def compute_rmt_from_packing(structure, safety_scale=0.97):
    """
    Вычисляет RMT через сферический packing.
    safety_scale < 1 гарантирует отсутствие overlap.
    """

    import numpy as np

    sites = structure.sites
    n = len(sites)

    # матрица расстояний
    min_dist = [1e9] * n

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = structure.get_distance(i, j)
            if d < min_dist[i]:
                min_dist[i] = d

    # Базовый RMT = половина минимального расстояния
    rmt = []
    for i, site in enumerate(sites):
        elem = site.specie.symbol
        base = base_radii.get(elem, 1.3)
        rmt_i = min(min_dist[i] / 2, base)
        rmt.append(rmt_i * safety_scale)

    return rmt


def get_rws_physical(structure, symmetrized, volume=None):
    import math

    if volume is None:
        volume = structure.volume

    # --- 1. RMT packing ---
    rmt_sites = compute_rmt_from_packing(structure)

    # --- 2. Усредняем по symm class ---
    rmt_class = []
    for eq_sites in symmetrized.equivalent_sites:
        vals = [rmt_sites[s.index] for s in eq_sites]
        rmt_class.append(sum(vals) / len(vals))

    # --- 3. Масштабируем чтобы заполнить объём ---
    total = 0.0
    total_sites = 0
    for r, eq_sites in zip(rmt_class, symmetrized.equivalent_sites):
        total += len(eq_sites) * r**3
        total_sites += len(eq_sites)

    sphere_vol = total * (4 * math.pi / 3)
    scale = (volume / sphere_vol) ** (1/3)

    rmt_class_scaled = [r * scale for r in rmt_class]

    # --- 4. Получаем RWS из RMT ---
    # типично RWS ≈ 1.12–1.18 * RMT
    rws_class = [r * 1.15 for r in rmt_class_scaled]

    return rmt_class_scaled, rws_class


def get_rws(structure: Structure, RBAS: np.ndarray, ALAT: float, rws_file_path='rws.vst', base_radii=None, volume_override=None):
    """
    Масштабирует радиусы Вигнера-Зейтца для 3D-систем, чтобы суммарный объём сфер соответствовал объёму ячейки.

    Args:
        structure (Structure): Объект структуры из pymatgen
        RBAS: Базис как в .sys файле
        ALAT: Постоянная решётки, заданная в .sys
        rws_file_path (str, optional): Путь к файлу rws.vst
        base_radii (dict, optional): Словарь базовых радиусов {элемент: радиус в Å}
        volume_override (float, optional): Переопределение объёма ячейки в Å³

    Returns:
        dict: Словарь радиусов Вигнера-Зейтца {элемент: радиус в Å}

    Raises:
        ValueError: Если сумма атомов не соответствует числу сайтов в структуре
    """
    # print(structure)
    # Константы
    PI = 3.141592653589793238462643

    # Получение объёма ячейки
    if volume_override is not None:
        volume = volume_override
    elif RBAS and ALAT:
        volume = abs(np.linalg.det(RBAS)) * ALAT**3
    else:
        # Вычисление объёма через векторы решётки: |a · (b × c)|
        volume = structure.volume


    # Получение состава и числа атомов
    composition = structure.composition
    natoms = {element.symbol: int(composition[element]) for element in composition}
    num_atoms = sum(natoms.values())
    if num_atoms != structure.num_sites:
        raise ValueError(f"Сумма атомов ({num_atoms}) не равна num_atoms ({structure.num_sites})")

    # Загрузка базовых радиусов
    if base_radii is None:
        base_radii = load_rws_vst(rws_file_path) if rws_file_path else {
            element.symbol: element.atomic_radius_calculated if element.atomic_radius_calculated else 1.3
            for element in composition
        }

    # Формирование словаря радиусов и числа атомов
    atomic_radii = {element: (base_radii.get(element, 1.3), natoms[element])
                    for element in natoms}


    # Сумма кубов базовых радиусов, умноженных на количество атомов
    total_r3 = sum(count * radius ** 3 for _, (radius, count) in atomic_radii.items())

    # Суммарный объём сфер до масштабирования
    v_rws = total_r3 * (4.0 * PI / 3.0)

    # Масштабный коэффициент k
    scale_rws = (volume / v_rws) ** (1.0 / 3.0)

    # Масштабированные радиусы
    ws_radii = {element: radius * scale_rws / 0.529177 for element, (radius, _) in atomic_radii.items()}

    return ws_radii


import math
from pymatgen.core import Structure


# Экспериментальные значения Skriver из rws.vst (колонка 1, по Z)
# 0 означает "нет эксп. значения" → используется rws.vst колонка 3
TABRWS_EXPERIMENTAL = {
    19: 4.862, 20: 4.122, 21: 3.427, 22: 3.052, 23: 2.818,
    24: 2.684, 25: 2.699, 26: 2.662, 27: 2.621, 28: 2.602, 29: 2.669,
    37: 5.197, 38: 4.494, 39: 3.761, 40: 3.347, 41: 3.071,
    42: 2.922, 43: 2.840, 44: 2.791, 45: 2.809, 46: 2.873, 47: 3.005,
    55: 5.656, 56: 4.652,
    57: 3.920, 58: 3.800, 59: 3.818, 60: 3.804, 61: 3.783,
    62: 3.768, 63: 4.263, 64: 3.764, 65: 3.720, 66: 3.704,
    67: 3.687, 68: 3.668, 69: 3.649, 70: 4.052, 71: 3.624,
    72: 3.301, 73: 3.069, 74: 2.945, 75: 2.872, 76: 2.825,
    77: 2.835, 78: 2.897, 79: 3.002,
    87: 5.900, 88: 4.790, 89: 3.900,
    90: 3.756, 91: 3.430, 92: 3.221, 93: 3.140, 94: 3.181,
    95: 3.614, 96: 3.641, 97: 3.550,
}


def get_rws_xband(
    structure: Structure,
    rws_vst_path: str = "rws.vst",
):
    """
    Вычисляет RWS по алгоритму xband "adjust R_WS by scaling":

    1. Начальный RWS каждого сайта = TABRWS[Z]:
       - сначала ищем экспериментальное значение (колонка 1 rws.vst / TABRWS_EXPERIMENTAL)
       - если 0 или отсутствует — берём колонка 3 (rws.vst fallback)
    2. scale = (V_cell / Σ(4π/3·rws³))^(1/3)  — итерация по всем сайтам
    3. RWS_final = rws_initial · scale

    Возвращает dict {element: rws_au} — одинаковый RWS для всех сайтов
    одного элемента, что гарантирует единый mesh в SPR-KKR.

    Args:
        structure:    примитивная структура pymatgen
        rws_vst_path: путь к файлу rws.vst

    Returns:
        dict {element_symbol: rws_scaled_au}
    """
    ANG_TO_AU = 1.889726125
    PI = math.pi

    # --- Шаг 1: загрузить rws.vst (колонки 1 и 3) ---
    rws_vst_exp = {}    # колонка 1: экспериментальные (по Z)
    rws_vst_calc = {}   # колонка 3: расчётные (по Z)

    try:
        with open(rws_vst_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    z = int(parts[0])
                except ValueError:
                    continue

                # колонка 1: экспериментальная (может отсутствовать)
                # колонка 3 (последняя перед символом элемента): расчётная
                # формат: Z  [exp_val]  calc_val  Symbol
                if len(parts) == 3:
                    # нет экспериментального значения: Z  calc_val  Symbol
                    try:
                        rws_vst_calc[z] = float(parts[1])
                    except ValueError:
                        pass
                elif len(parts) >= 4:
                    # есть экспериментальное: Z  exp_val  calc_val  Symbol
                    try:
                        exp_val = float(parts[1])
                        if exp_val > 0:
                            rws_vst_exp[z] = exp_val
                        rws_vst_calc[z] = float(parts[2])
                    except ValueError:
                        pass
    except FileNotFoundError:
        print(f"Файл {rws_vst_path} не найден. Используются встроенные значения.")

    # --- Шаг 2: назначить начальный RWS каждому элементу ---
    # Приоритет: rws.vst колонка 1 > TABRWS_EXPERIMENTAL > rws.vst колонка 3
    def get_initial_rws(z: int) -> float:
        # 1. экспериментальное из rws.vst
        if z in rws_vst_exp and rws_vst_exp[z] > 0:
            return rws_vst_exp[z]
        # 2. встроенное экспериментальное
        if z in TABRWS_EXPERIMENTAL and TABRWS_EXPERIMENTAL[z] > 0:
            return TABRWS_EXPERIMENTAL[z]
        # 3. расчётное из rws.vst (колонка 3)
        if z in rws_vst_calc and rws_vst_calc[z] > 0:
            return rws_vst_calc[z]
        # 4. запасное значение
        return 2.5

    elem_initial_rws = {
        el.symbol: get_initial_rws(el.Z)
        for el in structure.composition.elements
    }

    # --- Шаг 3: масштабирование (итерация по сайтам, как в xband) ---
    # V_RWS = Σ по всем сайтам rws_initial(site)³ · 4π/3
    v_rws = sum(
        elem_initial_rws[site.specie.symbol] ** 3
        for site in structure.sites
    ) * (4.0 * PI / 3.0)

    # Объём ячейки в а.е.³
    v_cell_au = structure.volume * ANG_TO_AU ** 3

    scale = (v_cell_au / v_rws) ** (1.0 / 3.0)

    # --- Шаг 4: итоговые значения ---
    rws_scaled = {
        elem: rws_init * scale
        for elem, rws_init in elem_initial_rws.items()
    }

    return rws_scaled
