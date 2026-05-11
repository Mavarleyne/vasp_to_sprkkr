from pathlib import Path
from typing import Dict, List
import numpy as np
from matplotlib import pyplot as plt


def main():
    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/Al/L21/vampire_new/*JXC.out')
    fig, ax = plot_exchanges_from_jxc(wd)
    plt.show()

    # wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ/Al/L21/vampire_manual/JXC.out')
    # fig, ax = plot_exchanges_from_jxc(wd)
    # plt.show()
    return

    wd = Path('/home/buche/VaspTesting/Danil/magnetocaloric_nn/SPR_KKR_Fe2CoZ')
    (wd / 'pics').mkdir(exist_ok=True)
    for p in wd.rglob('output'):
        path = p.parent
        fig, ax = plot_exchanges_from_jxc(path / '*JXC.out')

        name = f'{path.parts[-3]}_{path.parts[-2]}'
        fig.savefig(wd / 'pics' / name, dpi=300)
        plt.show()


class AtomType:
    idx: int
    label: str
    sites: List[int]

    def __init__(self, idx, label, sites):
        self.idx = idx
        self.label = label
        self.sites = sites

    def __repr__(self):
        return f'{{idx: {self.idx}, label: {self.label}, sites: {np.array(self.sites)}}}'

    def __str__(self):
        return f'Type_idx: {self.idx}, Label: {self.label}, on sites: {np.array(self.sites)}'


def parse_sprkkr_jxc(file_path: Path) -> Dict[str, np.ndarray]:
    """Парсер обменных констант J_ij из файлов SPR-KKR
    Параметры:
        file_path: путь к файлу вывода SPR-KKR.
    Возвращает:
        Словарь вида {'Fe_1-Fe_2': array([[DR_1, J_1_avg], ...])}
    """
    raw_data = {}
    types = {}
    jxc_out = file_path.read_text(encoding="utf-8").split("\n")

    if 'VERSION' in jxc_out[13]:
        kkr_ver = jxc_out[13].split()[3]

    if kkr_ver == '8.6.0':
        line_length = 13
        dr_pos = 10
        jij_pos = 12
        flag_line = 'IT   IQ   JT   JQ    N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [meV]  J_ij [eV]'
        table_head = 'type TXTT   NL VAL COR mesh    RMT     RWS   NAT  CONC  on sites'
        sites_start = 10
        table_len = 11
    elif kkr_ver == '7.7.2':
        line_length = 11
        dr_pos = 8
        jij_pos = 10
        flag_line = 'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]'
        table_head = 'type TXTT     NL mesh    RMT     RWS   NAT  CONC   on sites'
        sites_start = 8
        table_len = 9
    elif kkr_ver == '6.3.1':
        line_length = 11
        dr_pos = 8
        jij_pos = 10
        flag_line = 'ITAUIJ ITAUJI   N1 N2 N3    DRX    DRY    DRZ     DR     J_ij [Ry]  J_ij [eV]'
    else:
        raise RuntimeError('Unsupportable version of SPR-KKR, try to set other version')


    # --- ШАГ 1: Поиск таблицы с типами атомов ---
    flag = False
    for i, line in enumerate(jxc_out):
        aim = table_head

        if flag and len(line.split()) < table_len:
            break

        if flag:
            l = line.split()
            type_idx = int(l[0])
            label = l[1]
            sites = [int(i) for i in l[sites_start:]]
            types[type_idx] = AtomType(type_idx, label, sites)

        if aim in line:
            flag = True

    print(types)
    # exit()

    # --- ШАГ 2: Сбор данных J_ij (ваша исходная структура) ---
    for i in range(len(jxc_out)):
        # Ищем строку заголовка по ее началу
        if flag_line in jxc_out[i]:
            # Проверяем, есть ли следующая строка с данными
            data_line = jxc_out[i + 1].strip().split()
        # Проверяем, что строка содержит данные (минимум 13 колонок)
            try:
                it = int(jxc_out[i - 2].split()[5])  # IT
                jt = int(jxc_out[i - 2].split()[11])  # JT
                dr = float(data_line[dr_pos])  # DR
                jij = float(data_line[jij_pos]) * 1000  # J_ij [meV]

                # Переводим индексы в имена (например, '1' -> 'Fe_1')
                # print(types)
                atom1 = types[it].label
                atom2 = types[jt].label
                key = f"{atom1}-{atom2}"
                reverse_key = f'{atom2}-{atom1}'

                if reverse_key in raw_data.keys() and key != reverse_key:
                    continue

                if key not in raw_data:
                    raw_data[key] = []

                raw_data[key].append([dr, jij])
            except ValueError:
                continue

    # --- ШАГ 3: Усреднение данных ---
    result = {}
    for pair, values in raw_data.items():
        arr = np.array(values, dtype=np.float64)

        # Округляем расстояния до 4 знака
        # rounded_dr = np.round(arr[:, 0], decimals=4)
        dr = arr[:, 0]

        # Находим уникальные расстояния
        unique_dr, inverse_indices = np.unique(
            dr, return_inverse=True
        )

        # Суммируем и усредняем
        sum_jij = np.zeros_like(unique_dr)
        np.add.at(sum_jij, inverse_indices, arr[:, 1])
        counts = np.bincount(inverse_indices)
        avg_jij = sum_jij / counts

        result[pair] = np.column_stack((unique_dr, avg_jij))

    return result


def plot_exchanges_from_jxc(path_to_jxc: Path):
    '''

    :param path_to_jxc: Full path to *JXC.out
    :return: Figure
    '''
    jij = parse_sprkkr_jxc(path_to_jxc)
    fig, ax = plt.subplots(figsize=(9, 6))

    for pair, curve in jij.items():
        x = curve[:, 0]
        y = curve[:, 1]
        if np.allclose(y, np.zeros(y.shape), atol=0.2):
            continue

        ax.plot(x, y, label=pair, linestyle='-', marker='s')
    ax.set_xlim(0.4, 1.5)
    ax.set_xlabel('d/a')
    ax.set_ylabel(r'$J_{ij}$, meV')
    ax.grid(True)
    handles, labels = plt.gca().get_legend_handles_labels()

    sorted_labels_handles = sorted(zip(labels, handles), key=lambda t: t[0])
    labels, handles = zip(*sorted_labels_handles)

    ax.legend(handles, labels, ncol=3, fontsize=12)
    ax.legend(ncol=3, fontsize=12)
    plt.show()
    return fig, ax

if __name__ == '__main__':
    main()
