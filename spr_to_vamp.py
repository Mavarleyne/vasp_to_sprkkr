# -*- coding: utf-8 -*-
"""
sprkkr_to_vampire.py
====================
Генератор входных файлов VAMPIRE из выходных файлов SPR-KKR.

Использует sprkkr_parser для чтения данных расчёта (JXC.out, SCF.out)
и записывает vampire.UCF, vampire.mat, input.

Поддерживаемые версии SPR-KKR: 6.x, 7.x, 8.x.

Использование
-------------
Из командной строки::

    python sprkkr_to_vampire.py /path/to/calc_dir  --dr-max 10 --da-max 10.0

Из кода::

    from sprkkr_to_vampire import generate_vampire_inputs
    generate_vampire_inputs(Path('/path/to/calc_dir'), dr_max=10)

Структура каталога расчёта
--------------------------
В директории должны находиться файлы *JXC.out и *SCF.out.
Выходные файлы VAMPIRE пишутся в ту же директорию.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np

from parse_spr_outputs import read_jxc, read_scf

# ─────────────────────────────────────────────────────────────────────────────
# Глобальные параметры симуляции VAMPIRE (можно переопределять)
# ─────────────────────────────────────────────────────────────────────────────

DA_MAX          = 10.0    # максимальное расстояние Jij (в единицах a)
MACROCELL_ATOMS = 20_000  # целевое число атомов в суперячейке
T_MIN           = 0       # K
T_MAX           = 1800    # K
T_STEP          = 10      # K
MC_STEP         = 1000    # число шагов Монте-Карло на температурную точку


# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции
# ─────────────────────────────────────────────────────────────────────────────

def _get_macrosize(n_basis_atoms: int,
                   required_n_atoms: int = MACROCELL_ATOMS) -> tuple[int, int, int]:
    """
    Вычислить размер суперячейки (Nx, Ny, Nz) так, чтобы суммарное
    число атомов было близко к ``required_n_atoms``.
    """
    size = round((required_n_atoms / n_basis_atoms) ** (1 / 3))
    return size, size, size


def _find_file(directory: Path, pattern: str) -> Path:
    """
    Найти единственный файл по glob-паттерну в директории.

    Raises
    ------
    FileNotFoundError : если файл не найден.
    """
    hits = list(directory.glob(pattern))
    if not hits:
        raise FileNotFoundError(
            f'Файл "{pattern}" не найден в {directory}'
        )
    return hits[0]


# ─────────────────────────────────────────────────────────────────────────────
# Классы данных материалов
# ─────────────────────────────────────────────────────────────────────────────

class Material:
    """
    Один тип атома (материал) для VAMPIRE.

    Параметры
    ---------
    idx : int
        Номер типа, 1-based.
    label : str
        Метка («Fe», «Fe_1» и т.д.).
    spin : float
        Магнитный момент, μB. Знак сохраняется; VAMPIRE берёт абсолютное значение.
    element : str
        Химический символ элемента.
    sites : list[int]
        Номера атомных позиций (1-based), занятых этим типом.
    """

    def __init__(self, idx: int, label: str, spin: float,
                 element: str, sites: list[int]):
        self.idx = idx
        self.label = label
        self.spin = spin
        self.element = element
        self.sites = sites


class Materials:
    """
    Коллекция материалов, собранная из результатов парсера SPR-KKR.

    Строится из словарей, возвращённых :func:`sprkkr_parser.read_scf`.

    Parameters
    ----------
    scf_data : dict
        Результат :func:`sprkkr_parser.read_scf`.
    """

    def __init__(self, scf_data: dict):
        self._materials: list[Material] = []

        # Создаём lookup: idx_типа → m_spin (из последней итерации SCF)
        magmom_by_idx: dict[int, float] = {
            m['idx']: abs(m['m_spin'])
            for m in scf_data['magmoms']
        }

        for entry in scf_data['type_table']:
            idx     = entry['idx']
            name    = entry['name']
            element = name.split('_')[0]
            sites   = entry['sites']
            spin    = magmom_by_idx.get(idx, 0.0)

            self._materials.append(
                Material(idx=idx, label=name, spin=spin,
                         element=element, sites=sites)
            )

    def __iter__(self):
        return iter(self._materials)

    def __len__(self):
        return len(self._materials)

    def __repr__(self):
        lines = ['Materials:']
        for m in self._materials:
            lines.append(
                f'  [{m.idx}] {m.label:12s}  spin={m.spin:.4f} µB  '
                f'element={m.element}  sites={m.sites}'
            )
        return '\n'.join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Генераторы файлов VAMPIRE
# ─────────────────────────────────────────────────────────────────────────────

def write_mat_file(materials: Materials) -> str:
    """
    Сформировать содержимое файла ``vampire.mat``.

    Parameters
    ----------
    materials : Materials

    Returns
    -------
    str
        Текст файла.
    """
    lines = [f'material:num-materials = {len(materials)}']
    for mat in materials:
        lines += [
            '#---------------------------------------------------',
            f'# Material {mat.idx}  ({mat.label})',
            '#---------------------------------------------------',
            f'material[{mat.idx}]:material-name={mat.label}',
            f'material[{mat.idx}]:damping-constant=1.0',
        ]
        if abs(mat.spin) < 0.15:
            lines.append(f'material[{mat.idx}]:non-magnetic=keep')
        else:
            lines.append(f'material[{mat.idx}]:atomic-spin-moment={mat.spin:.6f} !muB')
        lines += [
            f'material[{mat.idx}]:uniaxial-anisotropy-constant=0.0',
            f'material[{mat.idx}]:material-element={mat.element}',
            f'material[{mat.idx}]:initial-spin-direction = 0.0,0.0,1.0',
            f'material[{mat.idx}]:uniaxial-anisotropy-direction = 0.0,0.0,1.0',
        ]
    return '\n'.join(lines)


def write_ucf_file(lattice: dict, materials: Materials, jij: np.ndarray) -> str:
    """
    Сформировать содержимое файла ``vampire.UCF``.

    Parameters
    ----------
    lattice : dict
        Результат :func:`sprkkr_parser.parse_lattice`.
    materials : Materials
    jij : np.ndarray
        Массив обменных интегралов, shape (N, 6):
        [IT-1, JT-1, N1, N2, N3, J_SI].

    Returns
    -------
    str
        Текст UCF-файла.
    """
    lens = lattice['lens_angstrom']          # (x, y, z) в Å
    prim = lattice['primitive_vecs']         # (3, 3) в ед. a
    basis = lattice['basis_vecs']            # (M, 3)

    n_atoms = basis.shape[0]
    n_mat   = len(materials)

    lines = ['# Unit cell size (Angstrom):']
    lines.append(f'{lens[0]:.6f} {lens[1]:.6f} {lens[2]:.6f}')
    lines.append('# Unit cell lattice vectors:')
    for row in prim:
        lines.append(f'{row[0]:10.6f}  {row[1]:10.6f}  {row[2]:10.6f}')

    lines.append('# Atoms')
    lines.append(f'{n_atoms} {n_mat}')

    # Для каждой позиции находим номер материала (0-based для VAMPIRE)
    sites_to_mat: dict[int, int] = {}
    for mat in materials:
        for s in mat.sites:
            sites_to_mat[s] = mat.idx - 1   # 0-based

    for i in range(n_atoms):
        mat_idx = sites_to_mat.get(i + 1, 0)
        bx, by, bz = basis[i]
        lines.append(f'{i} {bx:.5f} {by:.5f} {bz:.5f} {mat_idx}')

    lines.append('# Interactions')
    n_int = jij.shape[0]
    lines.append(f'{n_int} isotropic')

    for k, row in enumerate(jij):
        it, jt, n1, n2, n3, j = row
        lines.append(
            f'{k:d} {int(it):d} {int(jt):d} '
            f'{int(n1):d} {int(n2):d} {int(n3):d} {j:.6e}'
        )

    return '\n'.join(lines)


def write_input_file(lattice: dict, materials: Materials,
                     t_min: int = T_MIN, t_max: int = T_MAX,
                     t_step: int = T_STEP, mc_step: int = MC_STEP) -> str:
    """
    Сформировать содержимое файла ``input`` для VAMPIRE.

    Parameters
    ----------
    lattice : dict
        Результат :func:`sprkkr_parser.parse_lattice`.
    materials : Materials
    t_min, t_max, t_step : int
        Диапазон температур для кривой намагниченности (K).
    mc_step : int
        Число шагов Монте-Карло.

    Returns
    -------
    str
        Текст input-файла.
    """
    lens = lattice['lens_angstrom']
    basis = lattice['basis_vecs']
    n_basis = basis.shape[0]
    macro = _get_macrosize(n_basis)

    lines = [
        '#------------------------------------------',
        '# Creation attributes:',
        '#------------------------------------------',
        'create:full',
        'create:periodic-boundaries-x',
        'create:periodic-boundaries-y',
        'create:periodic-boundaries-z',
        '#------------------------------------------',
        'material:file=vampire.mat',
        'material:unit-cell-file = "vampire.UCF"',
        '#------------------------------------------',
        '# System Dimensions:',
        '#------------------------------------------',
        f'dimensions:unit-cell-size-x = {lens[0]:.6f} !A',
        f'dimensions:unit-cell-size-y = {lens[1]:.6f} !A',
        f'dimensions:unit-cell-size-z = {lens[2]:.6f} !A',
        f'dimensions:system-size-x = {0.1 * lens[0] * macro[0]:.4f} !nm',
        f'dimensions:system-size-y = {0.1 * lens[1] * macro[1]:.4f} !nm',
        f'dimensions:system-size-z = {0.1 * lens[2] * macro[2]:.4f} !nm',
        '#------------------------------------------',
        '# Simulation attributes:',
        '#------------------------------------------',
        f'sim:minimum-temperature = {t_min}',
        f'sim:maximum-temperature = {t_max}',
        f'sim:temperature-increment = {t_step}',
        f'sim:equilibration-time-steps = {mc_step // 10}',
        f'sim:loop-time-steps = {mc_step}',
        'sim:time-steps-increment = 1',
        '#------------------------------------------',
        '# Program and integrator details',
        '#------------------------------------------',
        'sim:program=curie-temperature',
        'sim:integrator = monte-carlo',
        '#------------------------------------------',
        '# Data output',
        '#------------------------------------------',
        '#output:real-time',
        'output:temperature',
        '#output:material-magnetisation',
        'output:magnetisation-length',
        'output:mean-total-energy',
    ]
    return '\n'.join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Основная функция генерации
# ─────────────────────────────────────────────────────────────────────────────

def generate_vampire_inputs(
    calc_dir: Path,
    *,
    da_max: float = DA_MAX,
    dr_max: int = -1,
    output_dir: Path | None = None,
    t_min: int = T_MIN,
    t_max: int = T_MAX,
    t_step: int = T_STEP,
    mc_step: int = MC_STEP,
    verbose: bool = True,
) -> Path:
    """
    Прочитать результаты расчёта SPR-KKR и записать входные файлы VAMPIRE.

    Parameters
    ----------
    calc_dir : Path
        Директория с файлами \*JXC.out и \*SCF.out.
    da_max : float
        Максимальное расстояние для Jij (в единицах *a*).
    dr_max : int
        Число координационных сфер (-1 = все).
    output_dir : Path, optional
        Куда записать файлы VAMPIRE. По умолчанию — ``calc_dir``.
    t_min, t_max, t_step : int
        Температурный диапазон (K).
    mc_step : int
        Число МК-шагов.
    verbose : bool
        Выводить информацию о процессе.

    Returns
    -------
    Path
        Директория, в которую записаны файлы.

    Raises
    ------
    FileNotFoundError
        Если *JXC.out или *SCF.out не найдены в ``calc_dir``.

    Examples
    --------
    # >>> from pathlib import Path
    # >>> from sprkkr_to_vampire import generate_vampire_inputs
    # >>> out = generate_vampire_inputs(Path('/data/Fe_L21'), dr_max=10)
    # >>> print(out)
    /data/Fe_L21
    """
    if isinstance(calc_dir, str):
        calc_dir = Path(calc_dir)

    if output_dir is None:
        output_dir = calc_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Чтение SPR-KKR файлов ─────────────────────────────────────────────
    jxc_path = _find_file(calc_dir, '*JXC.out')
    scf_path = _find_file(calc_dir, '*SCF.out')

    if verbose:
        print(f'[SPR-KKR→VAMPIRE] Директория: {calc_dir}')
        print(f'  JXC: {jxc_path.name}')
        print(f'  SCF: {scf_path.name}')

    jxc_data = read_jxc(jxc_path, da_max=da_max, dr_max_count=dr_max)
    scf_data = read_scf(scf_path)

    version  = jxc_data['version']
    lattice  = jxc_data['lattice']   # берём решётку из JXC (там всегда есть)
    jij      = jxc_data['jij']
    tc_mfa   = jxc_data['tc_mfa']

    if verbose:
        print(f'  SPR-KKR версия : {version}')
        print(f'  Ячейка (Å)     : {np.round(lattice["lens_angstrom"], 4)}')
        print(f'  Атомов в базисе: {lattice["basis_vecs"].shape[0]}')
        print(f'  Jij пар        : {jij.shape[0]}')
        if tc_mfa is not None:
            print(f'  T_C (MFA)      : {tc_mfa} K')

    # ── 2. Сборка материалов ─────────────────────────────────────────────────
    materials = Materials(scf_data)

    if verbose:
        print(f'\n{materials}')

    # ── 3. Генерация и запись файлов ─────────────────────────────────────────
    mat_text = write_mat_file(materials)
    ucf_text = write_ucf_file(lattice, materials, jij)
    inp_text = write_input_file(
        lattice, materials,
        t_min=t_min, t_max=t_max, t_step=t_step, mc_step=mc_step,
    )

    (output_dir / 'vampire.mat').write_text(mat_text + '\n', encoding='utf-8')
    (output_dir / 'vampire.UCF').write_text(ucf_text + '\n', encoding='utf-8')
    (output_dir / 'input').write_text(inp_text + '\n', encoding='utf-8')

    if verbose:
        print(f'\nЗаписаны файлы:')
        print(f'  {output_dir / "vampire.mat"}')
        print(f'  {output_dir / "vampire.UCF"}')
        print(f'  {output_dir / "input"}')

    return output_dir


def generate_vampire_inputs_recursive(
    root_path: Path,
    depth: int,
    dr_max: int,
    *,
    da_max: float = DA_MAX,
    subdir: str = 'vampire',
    **kwargs,
) -> None:
    """
    Рекурсивно обойти дерево каталогов и создать входные файлы VAMPIRE
    для каждого найденного расчёта SPR-KKR.

    Parameters
    ----------
    root_path : Path
        Корневая директория для обхода.
    depth : int
        Глубина вложенности *JXC.out относительно ``root_path``
        (число частей пути ``relative_to(root_path)`` минус 1).
    dr_max : int
        Число координационных сфер (-1 = все).
    da_max : float
        Максимальное расстояние Jij.
    subdir : str
        Имя поддиректории внутри каждого расчёта, куда копируются
        SPR-KKR файлы и пишутся выходные файлы VAMPIRE.
    **kwargs
        Дополнительные параметры для :func:`generate_vampire_inputs`
        (``t_min``, ``t_max``, ``t_step``, ``mc_step``, ``verbose``).

    Examples
    --------
    >>> generate_vampire_inputs_recursive(
    ...     Path('/data/SPR_KKR_Fe2CoZ'), depth=2, dr_max=-1
    ... )
    """
    for jxc_path in root_path.rglob('*JXC.out'):
        rel_parts = jxc_path.relative_to(root_path).parts
        if len(rel_parts) != depth + 1:
            continue

        calc_dir = jxc_path.parent
        out_dir = calc_dir / subdir

        # Создаём поддиректорию и копируем файлы расчёта
        out_dir.mkdir(exist_ok=True)
        for pattern in ('*JXC.out', '*SCF.out'):
            src_list = list(calc_dir.glob(pattern))
            if src_list:
                shutil.copy2(src_list[0], out_dir / src_list[0].name)

        try:
            generate_vampire_inputs(
                out_dir,
                da_max=da_max,
                dr_max=dr_max,
                **kwargs,
            )
        except Exception as exc:
            print(f'[ОШИБКА] {calc_dir}: {exc}')


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description='Генерация входных файлов VAMPIRE из выходных файлов SPR-KKR.'
    )
    p.add_argument('calc_dir', type=Path,
                   help='Директория с *JXC.out и *SCF.out')
    p.add_argument('--da-max', type=float, default=DA_MAX,
                   help=f'Макс. расстояние Jij в ед. a (по умолч. {DA_MAX})')
    p.add_argument('--dr-max', type=int, default=-1,
                   help='Число координационных сфер (-1 = все, по умолч.)')
    p.add_argument('--output-dir', type=Path, default=None,
                   help='Директория для записи файлов VAMPIRE (по умолч. = calc_dir)')
    p.add_argument('--t-min', type=int, default=T_MIN)
    p.add_argument('--t-max', type=int, default=T_MAX)
    p.add_argument('--t-step', type=int, default=T_STEP)
    p.add_argument('--mc-step', type=int, default=MC_STEP)
    return p.parse_args()


if __name__ == '__main__':
    # args = _parse_args()
    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/Al/L21/vampire_manual')
    generate_vampire_inputs(
        wd,
        dr_max=2
    )