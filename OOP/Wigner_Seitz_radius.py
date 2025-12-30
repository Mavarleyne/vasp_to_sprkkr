import os
from pymatgen.core import Structure


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
        with open(os.path.expanduser(rws_file_path), 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        element, radius = parts[0], float(parts[1])
                        rws_vst[element] = radius
    except FileNotFoundError:
        print(f"Файл {rws_file_path} не найден. Используются дефолтные радиусы.")
        # Дефолтные радиусы для Mn, Cr, Pt (на основе целевых значений)
        # rws_vst = {'Mn': 1.35, 'Cr': 1.34, 'Pt': 1.45}
    return rws_vst


def get_rws(structure: Structure, rws_file_path=None, base_radii=None, volume_override=None):
    """
    Масштабирует радиусы Вигнера-Зейтца для 3D-систем, чтобы суммарный объём сфер соответствовал объёму ячейки.

    Args:
        structure (Structure): Объект структуры из pymatgen
        rws_file_path (str, optional): Путь к файлу rws.vst
        base_radii (dict, optional): Словарь базовых радиусов {элемент: радиус в Å}
        volume_override (float, optional): Переопределение объёма ячейки в Å³

    Returns:
        dict: Словарь радиусов Вигнера-Зейтца {элемент: радиус в Å}

    Raises:
        ValueError: Если сумма атомов не соответствует числу сайтов в структуре
    """
    # Константы
    PI = 3.141592653589793238462643

    # Получение объёма ячейки
    if volume_override is not None:
        volume = volume_override
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
            element.symbol: element.atomic_radius if element.atomic_radius else 1.3
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