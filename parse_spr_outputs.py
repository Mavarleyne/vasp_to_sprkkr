# -*- coding: utf-8 -*-
"""
sprkkr_parser.py
================
Модуль для чтения выходных файлов пакета SPR-KKR (версии 6.x, 7.x, 8.x).

Поддерживаемые версии
---------------------
* 6.3.1  — формат Jij: ITAUIJ/ITAUJI + J_ij [Ry] + J_ij [eV], таблица типов без VAL/COR
* 7.7.x  — тот же формат, что 6.3.x
* 8.6.x  — формат Jij: IT IQ JT JQ + J_ij [meV] + J_ij [eV], таблица типов с VAL/COR

Функции
-------
detect_version(text_lines)          → str
parse_type_table_jxc(text_lines)    → list[dict]
parse_type_table_scf(text_lines)    → list[dict]
parse_lattice(text_lines)           → dict
parse_jij(text_lines, da_max, dr_max_count) → np.ndarray  shape (N,6)
parse_magmoms_scf(text_lines)       → list[dict]
parse_tc_mfa(text_lines)            → float | None

Высокоуровневые функции-фасады
------------------------------
read_jxc(path, da_max, dr_max_count) → dict
read_scf(path)                       → dict
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional

import numpy as np
import scipy.constants as constants

# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные паттерны (re.compile для скорости)
# ─────────────────────────────────────────────────────────────────────────────

# Строка с номером версии в шапке файла (строка ~14):
#   *       KKRGEN  VERSION  8.6.0          (C) 2021  H. Ebert   *
_RE_VERSION = re.compile(r'KKR\w+\s+VERSION\s+([\d.]+)', re.IGNORECASE)

# Заголовок таблицы типов атомов — два варианта:
#   v6/v7:  type TXTT     NL mesh    RMT     RWS   NAT  CONC   on sites
#   v8.6:   type TXTT   NL VAL COR mesh    RMT     RWS   NAT  CONC  on sites
_RE_TYPE_TABLE_HDR_NEW = re.compile(
    r'type\s+TXTT\s+NL\s+VAL\s+COR\s+mesh\s+RMT\s+RWS\s+NAT\s+CONC\s+on\s+sites'
)
_RE_TYPE_TABLE_HDR_OLD = re.compile(
    r'type\s+TXTT\s+NL\s+mesh\s+RMT\s+RWS\s+NAT\s+CONC\s+on\s+sites'
)

# Заголовок блока обменных интегралов — два варианта:
#   v8.6:  IT   IQ   JT   JQ    N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [meV]  J_ij [eV]
#   v6/v7: ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]
_RE_JIJ_HDR_NEW = re.compile(
    r'IT\s+IQ\s+JT\s+JQ\s+N1\s+N2\s+N3\s+DRX\s+DRY\s+DRZ\s+DR\s+J_ij \[meV\]'
)
_RE_JIJ_HDR_OLD = re.compile(
    r'ITAUIJ\s+ITAUJI\s+N1\s+N2\s+N3\s+DRX\s+DRY\s+DRZ\s+DR\s+J_ij \[Ry\]'
)

# Строка с IQ/IT/JQ/JT заголовком пары атомов (общая для всех версий):
#   IQ =  1 IT =  1                     JQ =  3 JT =  2
_RE_IQ_IT_PAIR = re.compile(
    r'IQ\s*=\s*(\d+)\s+IT\s*=\s*(\d+)\s+JQ\s*=\s*(\d+)\s+JT\s*=\s*(\d+)'
)

# Заголовок INIT_MOD_LATTICE
_RE_LATTICE_START = re.compile(r'<INIT_MOD_LATTICE>')
# Строка с примитивными векторами: "(  -0.50000,   0.50000,   0.50000 )"
_RE_PRIM_VEC = re.compile(
    r'\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)\s*\)'
)
# Конец секции векторов: "primitive vectors of reciprocal space Bravais lattice"
_RE_LATTICE_END = re.compile(r'2\*pi/a')

# Строка с постоянной решётки:
#   lattice constant  ALAT               5.35941
_RE_ALAT = re.compile(r'lattice\s+constant\s+ALAT\s+([\d.]+)')

# Секция магнитных моментов в SCF.out:
#   results extrapolated to corrected Fermi energy:
_RE_MAGMOM_BLOCK_START = re.compile(r'results extrapolated to corrected Fermi energy:')
# Строка с атомом над блоком DOS: "  31 E= 0.8023 0.0000  IT= 1  Fe"
_RE_ATOM_LINE = re.compile(
    r'\d+\s+E=\s*[\d.]+\s+[\d.]+\s+IT=\s*(\d+)\s+(\S+)'
)
# Строка суммы в блоке DOS:  "sum  10.8060   8.0000    4.6723   2.3035 ..."
_RE_SUM_LINE = re.compile(
    r'^[\s\w]*sum\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
)

# Температура Кюри в MFA:
#   Curie temperature within mean field approximation  T_C =   1458.2 K
_RE_TC_MFA = re.compile(
    r'Curie temperature within mean field approximation\s+T_C\s*=\s*([-\d.]+)\s*K'
)

# Константа для перевода из Ry в Дж: 1 Ry = constants.physical_constants['Rydberg constant times hc in J'][0]
# Но в исходном коде используется прямой eV→J: float(eV) * constants.e
# Для старых версий: J_Ry [Ry] → J_eV [eV] (уже есть в файле), затем умножаем на e*2
# Для новой версии: J_meV [meV] → J_eV [eV] (уже есть в файле)


# ─────────────────────────────────────────────────────────────────────────────
# Низкоуровневые парсеры
# ─────────────────────────────────────────────────────────────────────────────

def detect_version(lines: List[str]) -> str:
    """
    Определить версию SPR-KKR по первым строкам файла (JXC или SCF).

    Ищет строку вида:
        *       KKRGEN  VERSION  8.6.0 ...

    Parameters
    ----------
    lines : list[str]
        Строки файла (уже разбитые по ``\\n``).

    Returns
    -------
    str
        Строка версии, например ``'8.6.0'``, ``'7.7.2'``, ``'6.3.1'``.

    Raises
    ------
    RuntimeError
        Если версия не найдена в первых 50 строках.
    """
    for line in lines[:50]:
        m = _RE_VERSION.search(line)
        if m:
            return m.group(1)
    raise RuntimeError(
        'Не удалось определить версию SPR-KKR. '
        'Убедитесь, что файл содержит строку вида "KKRGEN  VERSION  X.Y.Z".'
    )


def _version_family(version: str) -> str:
    """
    Вернуть семейство версии: ``'new'`` (≥8) или ``'old'`` (6/7).
    """
    major = int(version.split('.')[0])
    return 'new' if major >= 8 else 'old'


# ─── Таблица типов атомов ────────────────────────────────────────────────────

def parse_type_table_jxc(lines: List[str], version: Optional[str] = None) -> List[dict]:
    """
    Прочитать таблицу типов атомов из JXC.out.

    Формат v6/v7::

        type TXTT     NL mesh    RMT     RWS   NAT  CONC   on sites
           1  Fe        3   1   2.3207  2.6388   1  1.000     1

    Формат v8.6::

        type TXTT   NL VAL COR mesh    RMT     RWS   NAT  CONC  on sites
           1  Al      3   3  10   1   2.3559  2.6869   1  1.000   1

    Parameters
    ----------
    lines : list[str]
    version : str, optional
        Если не передана — определяется автоматически.

    Returns
    -------
    list[dict]
        Каждый словарь содержит ключи:

        * ``idx``   — int, номер типа (1-based)
        * ``name``  — str, название типа («Fe», «Fe_1» и т.д.)
        * ``nat``   — int, количество атомов этого типа в ячейке
        * ``conc``  — float, концентрация
        * ``sites`` — list[int], номера позиций (1-based)
        * ``rmt``   — float, muffin-tin радиус
        * ``rws``   — float, Wigner-Seitz радиус
    """
    if version is None:
        version = detect_version(lines)
    family = _version_family(version)

    # Ищем заголовок таблицы
    hdr_re = _RE_TYPE_TABLE_HDR_NEW if family == 'new' else _RE_TYPE_TABLE_HDR_OLD
    table_start = None
    for i, line in enumerate(lines):
        if hdr_re.search(line):
            table_start = i + 1
            break

    if table_start is None:
        raise ValueError('Не найден заголовок таблицы типов атомов в JXC файле.')

    # Читаем строки таблицы до первой пустой / короткой строки
    result = []
    for line in lines[table_start:]:
        parts = line.split()
        # Строка данных: первый элемент — целое число (номер типа)
        if not parts or not parts[0].isdigit():
            if result:  # уже что-то прочитали → конец таблицы
                break
            continue

        idx = int(parts[0])
        name = parts[1]

        if family == 'new':
            # 1 Al 3 3 10 1 2.3559 2.6869 1 1.000 1 [2 3...]
            # idx name NL VAL COR mesh RMT RWS NAT CONC sites...
            rmt = float(parts[6])
            rws = float(parts[7])
            nat = int(parts[8])
            conc = float(parts[9])
            sites = [int(s) for s in parts[10:]]
        else:
            # 1 Fe 3 1 2.3207 2.6388 1 1.000 1 [2 ...]
            # idx name NL mesh RMT RWS NAT CONC sites...
            rmt = float(parts[4])
            rws = float(parts[5])
            nat = int(parts[6])
            conc = float(parts[7])
            sites = [int(s) for s in parts[8:]]

        result.append({
            'idx': idx,
            'name': name,
            'nat': nat,
            'conc': conc,
            'sites': sites,
            'rmt': rmt,
            'rws': rws,
        })

    return result


def parse_type_table_scf(lines: List[str], version: Optional[str] = None) -> List[dict]:
    """
    Прочитать таблицу типов атомов из SCF.out.

    Интерфейс и возвращаемый формат идентичны :func:`parse_type_table_jxc`.
    """
    # Формат таблиц в SCF.out идентичен JXC.out
    return parse_type_table_jxc(lines, version)


# ─── Параметры решётки ───────────────────────────────────────────────────────

def parse_lattice(lines: List[str]) -> dict:
    """
    Прочитать параметры решётки из JXC.out или SCF.out.

    Ищет секцию ``<INIT_MOD_LATTICE>``, извлекает примитивные векторы
    (в единицах *a*) и базисные векторы, а также постоянную решётки ALAT
    (в боровских радиусах).

    Parameters
    ----------
    lines : list[str]

    Returns
    -------
    dict с ключами:

    * ``alat``            — float, постоянная решётки (Bohr)
    * ``primitive_vecs``  — np.ndarray shape (3,3), вектора трансляций (в ед. *a*)
    * ``basis_vecs``      — np.ndarray shape (N,3), базис (в ед. *a*)
    * ``lens_angstrom``   — np.ndarray shape (3,), длины ребёр ячейки (Å)

    Notes
    -----
    Перевод из атомных единиц в Å: 1 Bohr = 0.52917721090 Å.
    """
    # Находим секцию INIT_MOD_LATTICE
    start = None
    for i, line in enumerate(lines):
        if _RE_LATTICE_START.search(line):
            start = i
            break
    if start is None:
        raise ValueError('<INIT_MOD_LATTICE> секция не найдена.')

    # Собираем координаты: сначала 3 примитивных вектора, потом базис
    # Прекращаем на строке "2*pi/a" (начало обратного пространства)
    all_vecs: List[List[float]] = []
    for line in lines[start:]:
        if _RE_LATTICE_END.search(line):
            break
        m = _RE_PRIM_VEC.search(line)
        if m:
            all_vecs.append([float(m.group(1)), float(m.group(2)), float(m.group(3))])

    if len(all_vecs) < 3:
        raise ValueError('Не удалось найти 3 примитивных вектора в секции INIT_MOD_LATTICE.')

    prim = np.array(all_vecs[:3])
    basis = np.array(all_vecs[3:]) if len(all_vecs) > 3 else np.zeros((1, 3))

    # Ищем ALAT (первое вхождение после начала файла, вне секции SCF)
    alat = None
    for line in lines:
        m = _RE_ALAT.search(line)
        if m:
            alat = float(m.group(1))
            break
    if alat is None:
        raise ValueError('Постоянная решётки ALAT не найдена.')

    # Длины ребёр примитивной ячейки в Å
    bohr_to_ang = 0.52917721090
    lens = np.array([
        alat * np.linalg.norm(prim[i]) * bohr_to_ang
        for i in range(3)
    ])

    return {
        'alat': alat,
        'primitive_vecs': np.round(prim, 5),
        'basis_vecs': np.round(basis, 5),
        'lens_angstrom': lens,
    }


# ─── Обменные интегралы Jij ──────────────────────────────────────────────────

def _make_jij_schema(version: str) -> dict:
    """
    Вернуть словарь с параметрами парсинга Jij-строки для данной версии.
    """
    family = _version_family(version)
    if family == 'new':
        # IT IQ JT JQ  N1 N2 N3  DRX DRY DRZ  DR  J_meV  J_eV
        return {
            'header_re': _RE_JIJ_HDR_NEW,
            'line_length': 13,
            'it_pos': 0,   # позиция IT в data-строке
            'iq_pos': 1,
            'jt_pos': 2,
            'jq_pos': 3,
            'n1_pos': 4,
            'dr_pos': 10,
            'jij_pos': 12,   # J_ij [eV]
            'jij_unit': 'eV',
        }
    else:
        # ITAUIJ ITAUJI  N1 N2 N3  DRX DRY DRZ  DR  J_Ry  J_eV
        return {
            'header_re': _RE_JIJ_HDR_OLD,
            'line_length': 11,
            'it_pos': None,   # IT/JT берём из строки заголовка пары
            'iq_pos': None,
            'jt_pos': None,
            'jq_pos': None,
            'n1_pos': 2,
            'dr_pos': 8,
            'jij_pos': 10,   # J_ij [eV]
            'jij_unit': 'eV',
        }


def parse_jij(
    lines: List[str],
    da_max: float = 10.0,
    dr_max_count: int = -1,
    version: Optional[str] = None,
) -> np.ndarray:
    """
    Прочитать обменные интегралы J_ij из JXC.out.

    Parameters
    ----------
    lines : list[str]
        Строки JXC.out.
    da_max : float
        Максимальное расстояние пары (в единицах *a*). Пары с DR > da_max
        пропускаются.
    dr_max_count : int
        Максимальное число координационных сфер. -1 = без ограничения.
        Сферы нумеруются глобально по возрастанию DR по всему файлу
        (не по порядку встречи строк, т.к. блоки разных пар IT/JT
        могут чередоваться в произвольном порядке DR).
    version : str, optional
        Версия SPR-KKR. Определяется автоматически, если не передана.

    Returns
    -------
    np.ndarray, shape (N, 6)
        Столбцы: ``[IT-1, JT-1, N1, N2, N3, J_ij_meV]``

        * IT-1, JT-1 — индексы типов (0-based, как в VAMPIRE)
        * N1, N2, N3 — целочисленные трансляции
        * J_ij_meV   — обменный интеграл в меВ

        Дублированные строки удалены, массив отсортирован.

    Raises
    ------
    RuntimeError
        Если версия не поддерживается.
    """
    if version is None:
        version = detect_version(lines)

    schema = _make_jij_schema(version)
    hdr_re = schema['header_re']
    ll = schema['line_length']
    n1_pos = schema['n1_pos']
    dr_pos = schema['dr_pos']
    jij_pos = schema['jij_pos']
    family = _version_family(version)

    # Текущая пара IT/JT (нужна для старых версий)
    curr_IT: float = 0.0
    curr_JT: float = 0.0

    # Сырые записи: [IT-1, JT-1, N1, N2, N3, J_meV, DR]
    # DR хранится отдельно для последующего глобального ранжирования
    raw: List[List[float]] = []

    i = 0
    while i < len(lines):
        line = lines[i]
        parts = line.split()

        if not parts:
            i += 1
            continue

        # Заголовок пары IQ=.. IT=.. JQ=.. JT=..  (все версии)
        m_pair = _RE_IQ_IT_PAIR.search(line)
        if m_pair:
            curr_IT = float(m_pair.group(2))
            curr_JT = float(m_pair.group(4))
            i += 1
            continue

        # Заголовок столбцов Jij-таблицы
        if hdr_re.search(line):
            i += 1
            if i >= len(lines):
                break
            data_parts = lines[i].split()
            if len(data_parts) < ll:
                i += 1
                continue

            try:
                dr    = float(data_parts[dr_pos])
                j_ev  = float(data_parts[jij_pos])
            except (ValueError, IndexError):
                i += 1
                continue

            # Фильтр по da_max применяем сразу — очевидно избыточные пары
            if dr > da_max:
                i += 1
                continue

            # Фильтр шумовых J
            # if abs(j_ev * 1000) < 1e-2:   # |J| < 0.01 meV
            #     i += 1
            #     continue

            n1 = float(data_parts[n1_pos])
            n2 = float(data_parts[n1_pos + 1])
            n3 = float(data_parts[n1_pos + 2])

            if family == 'new':
                IT = float(data_parts[schema['it_pos']])
                JT = float(data_parts[schema['jt_pos']])
            else:
                IT = curr_IT
                JT = curr_JT

            j_mev = j_ev * 1000
            j_si = j_ev * 2 * constants.e

            raw.append([IT - 1, JT - 1,  n1,  n2,  n3, j_mev, dr])
            raw.append([JT - 1, IT - 1, -n1, -n2, -n3, j_mev, dr])

        i += 1

    if not raw:
        return np.empty((0, 6))

    raw_arr = np.array(raw)  # shape (M, 7): cols 0-5 = данные, col 6 = DR

    # ── Глобальное ранжирование координационных сфер ──────────────────────────
    # DR может не возрастать монотонно по ходу файла (блоки разных пар IT/JT
    # чередуются в произвольном порядке), поэтому нумеруем сферы глобально:
    # собираем все уникальные DR, сортируем, присваиваем номера 1, 2, 3, ...
    if dr_max_count != -1:
        unique_drs = np.unique(np.round(raw_arr[:, 6], 4))
        dr_to_sphere = {dr: idx + 1 for idx, dr in enumerate(unique_drs)}
        sphere_nums = np.array([dr_to_sphere[round(dr, 4)] for dr in raw_arr[:, 6]])
        mask = sphere_nums < dr_max_count
        raw_arr = raw_arr[mask]

    if raw_arr.size == 0:
        return np.empty((0, 6))

    # Возвращаем только первые 6 столбцов (без DR)
    J_arr = np.array(raw_arr[:, :])
    # Сортировка и удаление дублей
    # idx = np.lexsort(J_arr.T)
    # J_arr = J_arr[idx]
    # mask = np.concatenate(([True], np.any(np.diff(J_arr, axis=0) != 0, axis=1)))
    # J_arr = J_arr[mask]
    out = sorted(J_arr, key=lambda x: (x[0], x[1]))
    # out = sorted(J_arr, key=lambda x: (x[-1]))
    J_arr = np.array(out)
    print(J_arr)
    return J_arr[:, :6]


# ─── Магнитные моменты из SCF.out ───────────────────────────────────────────
def parse_magmoms_scf(lines: List[str]) -> List[dict]:
    """
    Прочитать магнитные моменты из последней итерации SCF.out.

    Ищет **последний** блок ``results extrapolated to corrected Fermi energy:``
    и извлекает из него строки ``sum`` для каждого типа атомов.

    Parameters
    ----------
    lines : list[str]

    Returns
    -------
    list[dict]
        Каждый словарь:

        * ``idx``    — int, номер типа (1-based, из строки IT=)
        * ``name``   — str, название («Fe», «Co» и т.д.)
        * ``m_spin`` — float, спиновый момент (μB)
        * ``m_orb``  — float, орбитальный момент (μB)
        * ``m_total``— float, полный момент = m_spin + m_orb

    Notes
    -----
    Строка ``sum`` в блоке DOS SCF.out имеет вид::

        sum  10.8060   8.0000    4.6723   2.3035   0.00000  0.00000  ...

    Столбцы: DOS  NOS  P_spin  m_spin  P_orb  m_orb  ...
    """
    # Находим последний блок "results extrapolated..."
    last_block_start = None
    for i, line in enumerate(lines):
        if _RE_MAGMOM_BLOCK_START.search(line):
            last_block_start = i

    if last_block_start is None:
        raise ValueError('Блок "results extrapolated to corrected Fermi energy" не найден.')

    result: List[dict] = []
    current_it: Optional[int] = None
    current_name: Optional[str] = None
    in_dos_block = False

    for line in lines[last_block_start:]:
        parts = line.split()

        # Строка с IT= и именем атома: "  31 E= 0.8023 0.0000  IT= 1  Fe"
        m_atom = _RE_ATOM_LINE.search(line)
        if m_atom:
            current_it = int(m_atom.group(1))
            current_name = m_atom.group(2)
            in_dos_block = True
            continue

        # Строка sum (заканчивает блок DOS для одного атома)
        if in_dos_block and parts and parts[0] == 'sum':
            try:
                # sum DOS NOS P_spin m_spin P_orb m_orb ...
                # [0] [1]  [2]  [3]    [4]   [5]   [6]
                m_spin = float(parts[4])
                m_orb = float(parts[6]) if len(parts) > 6 else 0.0
            except (ValueError, IndexError):
                in_dos_block = False
                continue

            result.append({
                'idx': current_it,
                'name': current_name,
                'm_spin': m_spin,
                'm_orb': m_orb,
                'm_total': m_spin + m_orb,
            })
            in_dos_block = False
            current_it = None
            current_name = None

        # Блок итогов "TOT ..." — конец последней итерации
        if parts and parts[0] == 'TOT' and result:
            break

    return result


# ─── Температура Кюри MFA ───────────────────────────────────────────────────

def parse_tc_mfa(lines: List[str]) -> Optional[float]:
    """
    Прочитать температуру Кюри в приближении среднего поля из JXC.out.

    Parameters
    ----------
    lines : list[str]

    Returns
    -------
    float | None
        Значение T_C (K), или None если строка не найдена.
    """
    for line in reversed(lines):
        m = _RE_TC_MFA.search(line)
        if m:
            return float(m.group(1))
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Высокоуровневые фасады
# ─────────────────────────────────────────────────────────────────────────────

def read_jxc(
    path,
    da_max: float = 10.0,
    dr_max_count: int = -1,
) -> dict:
    """
    Прочитать JXC.out и вернуть все данные, нужные для конвертации в VAMPIRE.

    Parameters
    ----------
    path : str | Path
        Путь к JXC.out файлу (поддерживаются glob-паттерны типа ``'*JXC.out'``).
    da_max : float
        Максимальное расстояние Jij-пар (в единицах *a*).
    dr_max_count : int
        Число координационных сфер (-1 = все).

    Returns
    -------
    dict с ключами:

    * ``version``    — str
    * ``type_table`` — list[dict]  (см. :func:`parse_type_table_jxc`)
    * ``lattice``    — dict        (см. :func:`parse_lattice`)
    * ``jij``        — np.ndarray  (см. :func:`parse_jij`)
    * ``tc_mfa``     — float | None
    * ``path``       — Path
    """
    path = _resolve_path(path)
    lines = path.read_text(encoding='utf-8', errors='replace').split('\n')

    version = detect_version(lines)
    return {
        'version': version,
        'type_table': parse_type_table_jxc(lines, version),
        'lattice': parse_lattice(lines),
        'jij': parse_jij(lines, da_max, dr_max_count, version),
        'tc_mfa': parse_tc_mfa(lines),
        'path': path,
    }


def read_scf(path) -> dict:
    """
    Прочитать SCF.out и вернуть все данные, нужные для конвертации в VAMPIRE.

    Parameters
    ----------
    path : str | Path
        Путь к SCF.out (поддерживаются glob-паттерны).

    Returns
    -------
    dict с ключами:

    * ``version``    — str
    * ``type_table`` — list[dict]  (см. :func:`parse_type_table_scf`)
    * ``lattice``    — dict        (см. :func:`parse_lattice`)
    * ``magmoms``    — list[dict]  (см. :func:`parse_magmoms_scf`)
    * ``path``       — Path
    """
    path = _resolve_path(path)
    lines = path.read_text(encoding='utf-8', errors='replace').split('\n')

    version = detect_version(lines)
    return {
        'version': version,
        'type_table': parse_type_table_scf(lines, version),
        'lattice': parse_lattice(lines),
        'magmoms': parse_magmoms_scf(lines),
        'path': path,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные утилиты
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_path(path) -> Path:
    """
    Преобразовать строку/Path к Path, раскрывая glob-паттерн если нужно.
    """
    if isinstance(path, str):
        path = Path(path)
    if not path.exists():
        # Попробуем glob относительно CWD или родителя
        candidates = list(Path('.').glob(str(path)))
        if not candidates:
            raise FileNotFoundError(f'Файл не найден: {path}')
        path = candidates[0]
    return path


# ─────────────────────────────────────────────────────────────────────────────
# __all__
# ─────────────────────────────────────────────────────────────────────────────

__all__ = [
    'detect_version',
    'parse_type_table_jxc',
    'parse_type_table_scf',
    'parse_lattice',
    'parse_jij',
    'parse_magmoms_scf',
    'parse_tc_mfa',
    'read_jxc',
    'read_scf',
]

if __name__ == '__main__':
    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/Al/L21/*SCF_auto.out')
    lines = wd.read_text().split('\n')
    lattice = parse_lattice(lines)
    for key, value in lattice.items():
        print(f'{key}:\n{value}')