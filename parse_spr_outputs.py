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


# ─── Паттерны для секции OPTIMIZE_BASIS ─────────────────────────────────────

# Заголовок секции оптимизации базиса (только в SCF.out, только для NQ > 1)
_RE_OPT_BASIS_HDR = re.compile(r'<OPTIMIZE_BASIS_3D>')
# Заголовок строки-шапки таблицы OLD/NEW
_RE_OPT_BASIS_COLS = re.compile(r'OLD\s+basis\s+vectors\s+NEW\s+basis\s+vectors')


# ─── Параметры решётки ───────────────────────────────────────────────────────

def parse_lattice(lines: List[str]) -> dict:
    """
    Прочитать параметры решётки из JXC.out (примитивные векторы и NEW-базис).

    Читает секцию ``<INIT_MOD_LATTICE>`` и первое вхождение
    ``basis vectors in units of a`` (= NEW-базис после оптимизации).
    Постоянная решётки берётся из строки ``lattice constant  ALAT``.

    Parameters
    ----------
    lines : list[str]

    Returns
    -------
    dict с ключами:

    * ``alat``           — float, постоянная решётки (Bohr)
    * ``primitive_vecs`` — np.ndarray shape (3,3), прим. вектора (ед. *a*)
    * ``basis_vecs``     — np.ndarray shape (N,3), NEW-базис (ед. *a*)
    * ``lens_angstrom``  — np.ndarray shape (3,), длины ребёр ячейки (Å)

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

    prim  = np.array(all_vecs[:3])
    basis = np.array(all_vecs[3:]) if len(all_vecs) > 3 else np.zeros((1, 3))

    # Ищем ALAT (первое вхождение)
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
        alat * np.linalg.norm(prim[k]) * bohr_to_ang
        for k in range(3)
    ])

    return {
        'alat': alat,
        'primitive_vecs': np.round(prim, 5),
        'basis_vecs': np.round(basis, 5),
        'lens_angstrom': lens,
    }


def parse_basis_scf(lines: List[str]) -> dict:
    """
    Прочитать OLD- и NEW-координаты атомов из SCF.out и вычислить
    трансляции OLD→NEW в единицах примитивных векторов.

    Данные берутся из двух мест файла:

    1. ``<OPTIMIZE_BASIS_3D>`` — таблица ``OLD basis vectors / NEW basis vectors``
       (присутствует только при NQ > 1; если секция отсутствует, OLD = NEW = basis
       из ``<INIT_MOD_LATTICE>`` и трансляции нулевые).
    2. Примитивные векторы из первого блока ``primitive vectors for Bravais lattice``
       после ``<INIT_MOD_LATTICE>`` (уже содержатся в :func:`parse_lattice`, но
       здесь они нужны для вычисления трансляций).

    Трансляция для атома *i* определяется как вектор, переводящий его OLD-позицию
    в NEW-позицию, выраженный в единицах примитивных векторов::

        t_i = (new_i − old_i) @ inv(prim)   →  целые числа ∈ ℤ³

    Эти трансляции используются в :func:`parse_jij` для коррекции N1, N2, N3:
    при переходе от расчётного (NEW) базиса к исходному (OLD) базису вектор
    трансляции между атомами IQ и JQ корректируется как::

        N_corr = N_calc + t_JQ − t_IQ

    Parameters
    ----------
    lines : list[str]
        Строки SCF.out.

    Returns
    -------
    dict с ключами:

    * ``old_basis``     — np.ndarray shape (NQ, 3), исходные координаты (ед. *a*)
    * ``new_basis``     — np.ndarray shape (NQ, 3), оптимизированные координаты
    * ``translations``  — np.ndarray shape (NQ, 3), трансляции в прим. ед. (целые)
    * ``primitive_vecs``— np.ndarray shape (3, 3), прим. векторы (ед. *a*)

    Notes
    -----
    Если в SCF.out нет секции ``<OPTIMIZE_BASIS_3D>`` (однобазисные кристаллы),
    возвращаются нулевые трансляции и old_basis = new_basis.
    """
    # ── 1. Примитивные векторы из первого STRINIT-блока ──────────────────────
    # Используем уже готовый parse_lattice (он читает из INIT_MOD_LATTICE)
    lattice = parse_lattice(lines)
    prim = lattice['primitive_vecs']           # shape (3, 3)
    new_basis_init = lattice['basis_vecs']     # shape (NQ, 3) — NEW базис

    # ── 2. Секция OPTIMIZE_BASIS_3D ──────────────────────────────────────────
    opt_start = None
    for i, line in enumerate(lines):
        if _RE_OPT_BASIS_HDR.search(line):
            opt_start = i
            break

    if opt_start is None:
        # Секция отсутствует: OLD = NEW, трансляции = 0
        n = new_basis_init.shape[0]
        return {
            'old_basis': new_basis_init.copy(),
            'new_basis': new_basis_init.copy(),
            'translations': np.zeros((n, 3), dtype=float),
            'primitive_vecs': prim,
        }

    # Ищем строку-шапку "OLD basis vectors   NEW basis vectors"
    col_hdr = None
    for i in range(opt_start, min(opt_start + 20, len(lines))):
        if _RE_OPT_BASIS_COLS.search(lines[i]):
            col_hdr = i
            break
    if col_hdr is None:
        raise ValueError(
            'Заголовок "OLD basis vectors  NEW basis vectors" не найден '
            'в секции <OPTIMIZE_BASIS_3D>.'
        )

    # Читаем пары (old_vec, new_vec) — каждая строка содержит два вектора
    # Формат:  "(  x1, y1, z1 )   (  x2, y2, z2 )"
    _RE_TWO_VECS = re.compile(
        r'\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)\s*\)'
        r'.*?'
        r'\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)\s*\)'
    )
    old_list: List[List[float]] = []
    new_list: List[List[float]] = []

    for line in lines[col_hdr + 1:]:
        m = _RE_TWO_VECS.search(line)
        if m:
            old_list.append([float(m.group(1)), float(m.group(2)), float(m.group(3))])
            new_list.append([float(m.group(4)), float(m.group(5)), float(m.group(6))])
        elif old_list:
            # Первая непустая строка без двух векторов → конец таблицы
            stripped = line.strip()
            if stripped and not stripped.startswith('('):
                break

    if not old_list:
        raise ValueError('Строки OLD/NEW basis vectors не найдены.')

    old_basis = np.array(old_list)   # shape (NQ, 3)
    new_basis = np.array(new_list)   # shape (NQ, 3)

    # ── 3. Трансляции: t_i = (new_i − old_i) @ inv(prim) ────────────────────
    prim_inv = np.linalg.inv(prim)
    translations = (new_basis - old_basis) @ prim_inv   # shape (NQ, 3), должны быть ≈ целыми

    # Округляем до ближайшего целого (должны быть точно целыми ± числ. погрешность)
    translations = np.round(translations).astype(float)

    return {
        'old_basis': np.round(old_basis, 5),
        'new_basis': np.round(new_basis, 5),
        'translations': translations,
        'primitive_vecs': prim,
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
    translations: Optional[np.ndarray] = None,
    jij_Joul = True
) -> np.ndarray:
    """
    Прочитать обменные интегралы J_ij из JXC.out.

    Если переданы ``translations`` (полученные из :func:`parse_basis_scf`),
    то трансляционные индексы N1, N2, N3 корректируются для перехода от
    NEW-базиса (используемого в расчёте SPR-KKR) к OLD-базису (исходным
    физическим координатам атомов)::

        N_corr = N_calc + translations[JQ-1] − translations[IQ-1]

    где IQ, JQ — индексы сайтов (1-based) пары обменного взаимодействия.

    Parameters
    ----------
    lines : list[str]
        Строки JXC.out.
    da_max : float
        Жёсткий отсев по расстоянию: пары с DR > da_max пропускаются
        ещё при чтении файла.
    dr_max_count : int
        Максимальное число координационных сфер. -1 = без ограничения.

        Фильтрация выполняется **после** сбора всех данных и сортировки
        по DR: уникальные значения DR нумеруются начиная с 1 (1-я сфера —
        ближайшие соседи), и из итогового массива удаляются все строки,
        относящиеся к сферам с номером > dr_max_count.

        Такой подход корректен даже если в JXC.out записи идут не по
        возрастанию DR (что происходит при нескольких типах атомов).
    version : str, optional
        Версия SPR-KKR. Определяется автоматически, если не передана.
    translations : np.ndarray, shape (NQ, 3), optional
        Трансляции сайтов в единицах примитивных векторов — результат
        ``parse_basis_scf(scf_lines)['translations']``.
        Если None — коррекция не применяется (N остаются как в JXC.out).

    Returns
    -------
    np.ndarray, shape (N, 8)
        Столбцы: ``[IT-1, IQ-1, JT-1, JQ-1, N1, N2, N3, J_ij_SI]``

        * IT-1, JT-1 — индексы типов (0-based)
        * IQ-1, JQ-1 — индексы сайтов (0-based)
        * N1, N2, N3 — целочисленные трансляции (скорректированные, если
          передан параметр ``translations``); массив отсортирован по DR
        * J_ij_SI    — обменный интеграл в Дж: ``J_eV * e * 2``

        Дублированные строки удалены.

    Raises
    ------
    RuntimeError
        Если версия не поддерживается.
    """
    if version is None:
        version = detect_version(lines)

    schema = _make_jij_schema(version)
    hdr_re  = schema['header_re']
    ll      = schema['line_length']
    n1_pos  = schema['n1_pos']
    dr_pos  = schema['dr_pos']
    jij_pos = schema['jij_pos']
    family  = _version_family(version)

    use_translations = translations is not None

    J: List[List[float]] = []

    # Текущие индексы пары (обновляются при встрече заголовка "IQ= IT= JQ= JT=")
    curr_IT:  float = 0.0
    curr_JT:  float = 0.0
    curr_IQ:  float = 0.0
    curr_JQ:  float = 0.0

    i = 0
    while i < len(lines):
        line  = lines[i]
        parts = line.split()

        if not parts:
            i += 1
            continue

        # ── Заголовок пары: "IQ =  1 IT =  1   JQ =  2 JT =  1" ─────────────
        m_pair = _RE_IQ_IT_PAIR.search(line)
        if m_pair:
            curr_IQ = float(m_pair.group(1))
            curr_IT = float(m_pair.group(2))
            curr_JQ = float(m_pair.group(3))
            curr_JT = float(m_pair.group(4))
            i += 1
            continue

        # ── Заголовок столбцов Jij-блока ─────────────────────────────────────
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

            if dr > da_max:
                i += 1
                continue
            if abs(j_ev * 1000) < 1e-4:   # |J| < 0.0001 meV → пропускаем
                i += 1
                continue

            # N1, N2, N3 из строки данных
            n1 = float(data_parts[n1_pos])
            n2 = float(data_parts[n1_pos + 1])
            n3 = float(data_parts[n1_pos + 2])

            # IT, JT, IQ, JQ
            if family == 'new':
                IT = float(data_parts[schema['it_pos']])
                IQ = float(data_parts[schema['iq_pos']])
                JT = float(data_parts[schema['jt_pos']])
                JQ = float(data_parts[schema['jq_pos']])
            else:
                IT = curr_IT
                JT = curr_JT
                IQ = curr_IQ
                JQ = curr_JQ

            # ── Коррекция N на трансляции OLD→NEW ────────────────────────────
            # В JXC используется NEW-базис: N*prim + new[JQ] - new[IQ] = r
            # В VAMPIRE нужен OLD-базис:    N'*prim + old[JQ] - old[IQ] = r
            # Откуда:  N' = N + t_JQ - t_IQ,  t_k = (new_k - old_k) @ inv(prim)
            if use_translations:
                iq_idx = int(IQ) - 1
                jq_idx = int(JQ) - 1
                t_IQ = translations[iq_idx]
                t_JQ = translations[jq_idx]
                n1 += t_JQ[0] - t_IQ[0]
                n2 += t_JQ[1] - t_IQ[1]
                n3 += t_JQ[2] - t_IQ[2]

            if jij_Joul:
                j_si = j_ev * constants.e * 2  # Дж
            else:
                j_si = j_ev
            # dr сохраняем последним столбцом — нужен для фильтрации
            # по координационным сферам после сборки всего массива
            J.append([IT - 1, IQ - 1, JT - 1, JQ - 1,  n1,  n2,  n3, j_si, dr])
            J.append([JT - 1, JQ - 1, IT - 1, IQ - 1, -n1, -n2, -n3, j_si, dr])

        i += 1

    if not J:
        return np.empty((0, 8))

    J_arr = np.array(J)   # shape (N, 9): [..., dr]

    # ── Сортировка по dr, затем дедупликация ─────────────────────────────────
    # Сортируем по (dr, остальным столбцам) — dr последний, но lexsort читает
    # с конца, поэтому dr (столбец 8) будет первичным ключом сортировки.
    idx   = np.lexsort(J_arr.T)
    J_arr = J_arr[idx]
    mask  = np.concatenate(([True], np.any(np.diff(J_arr, axis=0) != 0, axis=1)))
    J_arr = J_arr[mask]

    # ── Фильтрация по числу координационных сфер ─────────────────────────────
    # Теперь массив отсортирован по dr → уникальные значения dr идут по возрастанию.
    # Координационная сфера = группа строк с одинаковым dr (с допуском 1e-4).
    if dr_max_count != -1:
        dr_vals = J_arr[:, -1]            # столбец dr
        sphere_count = 0
        prev_dr      = -1.0
        cutoff_idx   = len(J_arr)         # по умолчанию — брать всё

        for k, dr in enumerate(dr_vals):
            if abs(dr - prev_dr) > 1e-4:
                sphere_count += 1
                prev_dr = dr
                if sphere_count > dr_max_count:
                    cutoff_idx = k
                    break

        J_arr = J_arr[:cutoff_idx]

    # Убираем вспомогательный столбец dr и возвращаем итоговый массив
    return J_arr[:, :-1]


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
    scf_path=None,
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
    scf_path : str | Path, optional
        Путь к соответствующему SCF.out. Если передан, из него читаются
        OLD/NEW базисы и трансляции, которые применяются при извлечении
        N1, N2, N3 в :func:`parse_jij`. Поддерживаются glob-паттерны.

    Returns
    -------
    dict с ключами:

    * ``version``      — str
    * ``type_table``   — list[dict]  (см. :func:`parse_type_table_jxc`)
    * ``lattice``      — dict        (см. :func:`parse_lattice`)
    * ``basis_info``   — dict | None (см. :func:`parse_basis_scf`; None если scf_path не передан)
    * ``jij``          — np.ndarray  (см. :func:`parse_jij`)
    * ``tc_mfa``       — float | None
    * ``path``         — Path
    """
    path = _resolve_path(path)
    lines = path.read_text(encoding='utf-8', errors='replace').split('\n')

    version    = detect_version(lines)
    lattice    = parse_lattice(lines)
    type_table = parse_type_table_jxc(lines, version)

    # Трансляции из SCF — если передан путь
    basis_info = None
    translations = None
    if scf_path is not None:
        scf_p = _resolve_path(scf_path)
        scf_lines = scf_p.read_text(encoding='utf-8', errors='replace').split('\n')
        basis_info   = parse_basis_scf(scf_lines)
        translations = basis_info['translations']

    jij = parse_jij(lines, da_max, dr_max_count, version, translations)

    return {
        'version':    version,
        'type_table': type_table,
        'lattice':    lattice,
        'basis_info': basis_info,
        'jij':        jij,
        'tc_mfa':     parse_tc_mfa(lines),
        'path':       path,
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

    * ``version``      — str
    * ``type_table``   — list[dict]  (см. :func:`parse_type_table_scf`)
    * ``lattice``      — dict        (см. :func:`parse_lattice`)
    * ``basis_info``   — dict        (см. :func:`parse_basis_scf`)
    * ``magmoms``      — list[dict]  (см. :func:`parse_magmoms_scf`)
    * ``path``         — Path
    """
    path = _resolve_path(path)
    lines = path.read_text(encoding='utf-8', errors='replace').split('\n')

    version = detect_version(lines)
    return {
        'version':    version,
        'type_table': parse_type_table_scf(lines, version),
        'lattice':    parse_lattice(lines),
        'basis_info': parse_basis_scf(lines),
        'magmoms':    parse_magmoms_scf(lines),
        'path':       path,
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
    'parse_basis_scf',
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

    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/Al/L21/*JXC_auto.out')
    lines = wd.read_text().split('\n')
    lattice = parse_lattice(lines)
    for key, value in lattice.items():
        print(f'{key}:\n{value}')