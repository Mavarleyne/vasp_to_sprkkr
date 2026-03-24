import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings

from scipy.signal import savgol_filter
from scipy.optimize import curve_fit


# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции (общие для обоих вариантов)
# ─────────────────────────────────────────────────────────────────────────────

def critical_law(T, A, Tc, beta):
    """M = A · (1 − T/Tc)^beta,  0 при T >= Tc."""
    ratio = np.asarray(1.0 - T / Tc, dtype=float)
    return np.where(ratio > 0, A * ratio ** beta, 0.0)


def critical_law_fixed_beta(T, A, Tc, beta=0.365):
    """Критический закон с фиксированным beta (по умолчанию 0.365)."""
    ratio = np.asarray(1.0 - T / Tc, dtype=float)
    return np.where(ratio > 0, A * ratio ** beta, 0.0)


def _prepare_data(data):
    """
    Общий препроцессинг: дедупликация, сортировка, |M|, сглаживание.
    Возвращает (T, M, M_smooth) или (None, None, None) при ошибке.
    """
    if data is None or len(data) < 10:
        warnings.warn(
            "curie_critical_fit: слишком мало точек данных (< 10).",
            UserWarning,
        )
        return None, None, None

    T_unique, idx = np.unique(data[:, 0], return_index=True)
    M_unique = np.abs(data[:, 1][idx])
    order = np.argsort(T_unique)
    T = T_unique[order]
    M = M_unique[order]

    n = len(T)
    window = min(21, n if n % 2 == 1 else n - 1)
    if window < 5:
        warnings.warn(
            "curie_critical_fit: слишком мало уникальных точек для сглаживания.",
            UserWarning,
        )
        return None, None, None

    poly = min(3, window - 2)
    M_smooth = np.clip(savgol_filter(M, window_length=window, polyorder=poly), 0, None)
    return T, M, M_smooth


def _estimate_Tc(T, M_smooth):
    """
    Первичная оценка Tc по максимуму |dM/dT|.
    Защита от краевых эффектов — запасной вариант через медиану.
    """
    dT = np.gradient(T)
    dT = np.where(np.abs(dT) < 1e-12, 1e-12, dT)
    dM_dT = np.gradient(M_smooth) / dT

    Tc_est = float(T[np.argmax(np.abs(dM_dT))])
    T_range = T[-1] - T[0]
    edge = 0.05 * T_range

    if Tc_est < T[0] + edge or Tc_est > T[-1] - edge:
        mid = (T > T[0] + 0.1 * T_range) & (T < T[-1] - 0.1 * T_range)
        if mid.sum() >= 5:
            Tc_est = float(T[mid][np.argmax(np.abs(dM_dT[mid]))])
        else:
            Tc_est = float(np.median(T))

    return Tc_est


def _fit_bounds_and_starts(T, M_smooth, Tc_est, beta_fixed=None):
    """
    Возвращает: T_fit, M_fit, lower_bounds, upper_bounds, trials
    trials — список стартовых (A, Tc, beta).
    beta_fixed=None → beta свободна; иначе beta зафиксирована.
    """
    T_range = T[-1] - T[0]
    M_threshold = 0.03 * float(M_smooth.max())

    # Сужаем диапазон фита: [0.75·Tc_est, Tc_est]  — критическая область
    lower_T_cutoff = Tc_est - 0.25 * (Tc_est - T[0])
    mask = (T >= lower_T_cutoff) & (T < Tc_est) & (M_smooth > M_threshold)
    T_fit = T[mask]
    M_fit = M_smooth[mask]

    if len(T_fit) < 8:
        # Запасной вариант: весь диапазон ниже Tc_est
        mask = (T < Tc_est) & (M_smooth > M_threshold)
        T_fit = T[mask]
        M_fit = M_smooth[mask]

    T_step    = float(np.median(np.diff(T)))
    T_fit_max = float(T_fit.max()) if len(T_fit) else Tc_est
    A0        = float(M_fit.max()) if len(M_fit) else 1.0

    Tc_lo = T_fit_max + T_step
    Tc_hi = max(T[-1] + 0.20 * T_range, Tc_lo + 1.0)
    Tc0   = float(np.clip(Tc_est, Tc_lo + 0.1, Tc_hi - 0.1))

    if beta_fixed is not None:
        # beta не оптимизируется — bounds/trials только для A и Tc
        lower_bounds = [0.0,          Tc_lo]
        upper_bounds = [max(A0*2, 1), Tc_hi]
        rng = np.random.default_rng(42)
        trials = [
            [A0 * rng.uniform(0.9, 1.1), float(np.clip(Tc0 + rng.uniform(-T_step, T_step), Tc_lo+0.1, Tc_hi-0.1))]
            for _ in range(20)
        ]
    else:
        lower_bounds = [0.0,          Tc_lo,  0.10]
        upper_bounds = [max(A0*2, 1), Tc_hi,  0.60]
        rng = np.random.default_rng(42)
        beta_trials = [0.33, 0.36, 0.25, 0.50]
        trials = [
            [
                A0  * rng.uniform(0.9, 1.1),
                float(np.clip(Tc0 + rng.uniform(-T_step, T_step), Tc_lo+0.1, Tc_hi-0.1)),
                float(np.clip(b * rng.uniform(0.9, 1.1), 0.11, 0.59)),
            ]
            for b in beta_trials for _ in range(5)
        ]

    return T_fit, M_fit, lower_bounds, upper_bounds, trials


def _save_plot(T, M, M_smooth, T_fit, M_fit, T_model, M_model,
               Tc, tc_mfa, beta, beta_label, save_path, note=""):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(T, M,        color="steelblue", alpha=0.4, lw=1,   label="Исходные данные")
    ax.plot(T, M_smooth, color="steelblue", lw=2,  ls="--",    label="Сглаженные данные")
    ax.scatter(T_fit, M_fit, s=15, color="steelblue", alpha=0.7, zorder=3,
               label=f"Точки фита ({len(T_fit)} шт.)")
    ax.plot(T_model, M_model, color="crimson", lw=2.5,
            label=rf"Критический фит: $\beta={beta_label}$")
    ax.axvline(Tc, color="crimson", ls=":", lw=1.5)
    ax.scatter([Tc], [0], color="crimson", zorder=6, s=80)
    ax.annotate(
        f"$T_C = {Tc:.1f}$ K",
        xy=(Tc, 0), xytext=(12, 16), textcoords="offset points",
        fontsize=11, color="crimson",
        arrowprops=dict(arrowstyle="->", color="crimson", lw=1.2),
    )
    fig.text(0.8, 0.8, r'$T_C^{MFA}$ = ' + f'{tc_mfa:.2f}' + ' K', backgroundcolor='#C8C8C8')
    if note:
        ax.text(0.02, 0.05, note, transform=ax.transAxes,
                fontsize=9, color="darkorange",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", edgecolor="orange"))
    ax.set_xlabel("Температура (K)", fontsize=12)
    ax.set_ylabel("Намагниченность (норм.)", fontsize=12)
    ax.set_title(rf"Температура Кюри — критическая аппроксимация"
                 "\n" rf"$M = A\,(1 - T/T_C)^{{\beta}}$", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.4)
    ax.set_ylim(bottom=-0.02)
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"График сохранён: {save_path}")


# ═════════════════════════════════════════════════════════════════════════════
# ВАРИАНТ 1 — beta фиксирован = 0.365 (теоретическое значение Гейзенберга 3D)
# ═════════════════════════════════════════════════════════════════════════════

BETA_HEISENBERG = 0.365   # теоретическое значение для 3D-модели Гейзенберга


def curie_fit_fixed_beta(data: np.ndarray, Tc_mfa: float, save_path: str, beta: float = BETA_HEISENBERG):
    """
    Определяет температуру Кюри T_C при ФИКСИРОВАННОМ критическом показателе beta.

    Beta задаётся извне (по умолчанию 0.365 — 3D-Гейзенберг).
    Оптимизируются только два параметра: амплитуда A и T_C.

    Преимущества фиксированного beta:
      - Меньше степеней свободы → более устойчивый фит при малом числе точек.
      - Физически обоснован: при уверенности в модели не нужно подбирать beta.
      - Исключает ложные решения с beta > 0.5.

    Параметры
    ----------
    data      : np.ndarray (N, 2) — [T (K), M (норм.)]
    save_path : str
    beta      : float — фиксированное значение критического показателя

    Возвращает
    ----------
    (Tc, beta) | (None, None)
    """

    T, M, M_smooth = _prepare_data(data)
    if T is None:
        return None, None

    Tc_est = _estimate_Tc(T, M_smooth)
    T_fit, M_fit, lower_bounds, upper_bounds, trials = \
        _fit_bounds_and_starts(T, M_smooth, Tc_est, beta_fixed=beta)

    if len(T_fit) < 8:
        warnings.warn(
            f"curie_fit_fixed_beta: недостаточно точек для фита ({len(T_fit)}). "
            f"Возвращаем оценку по производной T_C ≈ {Tc_est:.1f} K.",
            UserWarning,
        )
        return Tc_est, beta

    # Обёртка: подставляем фиксированный beta внутрь двухпараметрической функции
    def model(T_, A_, Tc_):
        return critical_law_fixed_beta(T_, A_, Tc_, beta)

    popt = None
    best_rss = np.inf

    for p0 in trials:
        try:
            p, cov = curve_fit(
                model, T_fit, M_fit,
                p0=p0,
                bounds=(lower_bounds, upper_bounds),
                maxfev=50_000,
                method="trf",
            )
            rss = float(np.sum((model(T_fit, *p) - M_fit) ** 2))
            if rss < best_rss:
                best_rss = rss
                popt = p
                best_cov = cov
        except RuntimeError:
            continue

    if popt is None:
        warnings.warn(
            f"curie_fit_fixed_beta: curve_fit не сошёлся. "
            f"Возвращаем оценку по производной T_C ≈ {Tc_est:.1f} K.",
            UserWarning,
        )
        return Tc_est, beta

    A, Tc = popt
    try:
        Tc_err = float(np.sqrt(best_cov[1, 1]))
    except Exception:
        Tc_err = float("nan")

    print(f"[fixed beta] T_C    = {Tc:.2f} ± {Tc_err:.2f} K")
    print(f"[fixed beta] Tc_MFA = {Tc_mfa:.2f} K")
    print(f"[fixed beta] beta   = {beta:.4f}  (зафиксирован, 3D-Гейзенберг)")
    print(f"[fixed beta] A      = {A:.4f}")

    T_model = np.linspace(T_fit.min(), Tc * 0.9999, 500)
    M_model = critical_law_fixed_beta(T_model, A, Tc, beta)
    _save_plot(T, M, M_smooth, T_fit, M_fit, T_model, M_model,
               Tc, Tc_mfa, beta, f"{beta:.3f} (фикс.)", save_path)

    return Tc, beta


# ═════════════════════════════════════════════════════════════════════════════
# ВАРИАНТ 2 — beta свободна, но с обработкой случая beta > 0.5
# ═════════════════════════════════════════════════════════════════════════════

BETA_MAX_PHYSICAL = 0.50   # граница физической достоверности
BETA_HEISENBERG_3D = 0.365


def curie_fit_free_beta(data: np.ndarray, Tc_mfa: float, save_path: str):
    """
    Определяет T_C и beta методом нелинейного МНК.
    Beta оптимизируется свободно в диапазоне [0.10, 0.60].

    Если оптимальный beta > 0.50 (выход за теорию среднего поля):
      1. Выдаётся предупреждение с диагнозом причины.
      2. Автоматически выполняется повторный фит в суженном диапазоне
         температур [0.85·T_C, T_C] — только критическая область.
      3. Если beta снова > 0.50 — результат помечается как ненадёжный,
         но возвращается вместе с флагом и предупреждением (не исключение).
      4. На графике появляется аннотация-предупреждение.

    Возвращает
    ----------
    (Tc, beta, reliable: bool)
        reliable=False означает, что beta > 0.5 даже после сужения диапазона.
    """

    T, M, M_smooth = _prepare_data(data)
    if T is None:
        return None, None, False

    Tc_est = _estimate_Tc(T, M_smooth)

    def _run_fit(T_fit, M_fit):
        """Внутренний запуск curve_fit, возвращает (popt, pcov) или None."""
        if len(T_fit) < 8:
            return None, None

        T_step    = float(np.median(np.diff(T)))
        T_range   = T[-1] - T[0]
        A0        = float(M_fit.max())
        T_fit_max = float(T_fit.max())

        Tc_lo = T_fit_max + T_step
        Tc_hi = max(T[-1] + 0.20 * T_range, Tc_lo + 1.0)
        Tc0   = float(np.clip(Tc_est, Tc_lo + 0.1, Tc_hi - 0.1))

        lower_bounds = [0.0,          Tc_lo,  0.10]
        upper_bounds = [max(A0*2, 1), Tc_hi,  0.60]

        rng = np.random.default_rng(42)
        beta_trials = [0.33, 0.36, 0.25, 0.50]
        trials = [
            [
                A0 * rng.uniform(0.9, 1.1),
                float(np.clip(Tc0 + rng.uniform(-T_step, T_step), Tc_lo+0.1, Tc_hi-0.1)),
                float(np.clip(b * rng.uniform(0.9, 1.1), 0.11, 0.59)),
            ]
            for b in beta_trials for _ in range(5)
        ]

        best_popt, best_pcov, best_rss = None, None, np.inf
        for p0 in trials:
            try:
                p, cov = curve_fit(
                    critical_law, T_fit, M_fit,
                    p0=p0,
                    bounds=(lower_bounds, upper_bounds),
                    maxfev=50_000,
                    method="trf",
                )
                rss = float(np.sum((critical_law(T_fit, *p) - M_fit) ** 2))
                if rss < best_rss:
                    best_rss = rss
                    best_popt = p
                    best_pcov = cov
            except RuntimeError:
                continue
        return best_popt, best_pcov

    # ── Попытка 1: стандартный диапазон фита ─────────────────────────────────
    T_range   = T[-1] - T[0]
    M_threshold = 0.03 * float(M_smooth.max())
    lower_T1  = Tc_est - 0.25 * (Tc_est - T[0])
    mask1 = (T >= lower_T1) & (T < Tc_est) & (M_smooth > M_threshold)
    if mask1.sum() < 8:
        mask1 = (T < Tc_est) & (M_smooth > M_threshold)

    T_fit1 = T[mask1]
    M_fit1 = M_smooth[mask1]

    popt, pcov = _run_fit(T_fit1, M_fit1)

    if popt is None:
        warnings.warn(
            f"curie_fit_free_beta: curve_fit не сошёлся. "
            f"Возвращаем оценку по производной T_C ≈ {Tc_est:.1f} K.",
            UserWarning,
        )
        return Tc_est, None, False

    A, Tc, beta = popt
    reliable = True
    note = ""

    # ── Проверка beta > 0.5 ───────────────────────────────────────────────────
    if beta > BETA_MAX_PHYSICAL:
        # Диагностика причины
        n_fit = len(T_fit1)
        frac  = (Tc_est - T_fit1.min()) / (Tc_est - T[0]) if Tc_est > T[0] else 1.0

        reason_parts = []
        if frac > 0.5:
            reason_parts.append("слишком широкий диапазон фита (включены точки далеко от T_C)")
        if n_fit < 20:
            reason_parts.append(f"мало точек вблизи T_C ({n_fit} шт.)")
        if not reason_parts:
            reason_parts.append("возможны конечно-размерные эффекты или шум данных")

        reason_str = "; ".join(reason_parts)
        warnings.warn(
            f"curie_fit_free_beta: beta = {beta:.4f} > {BETA_MAX_PHYSICAL} "
            f"(выход за предел теории среднего поля). Причина: {reason_str}. "
            f"Выполняется повторный фит в суженном диапазоне [0.85·T_C, T_C].",
            UserWarning,
        )

        # ── Попытка 2: только критическая область [0.85·Tc_est, Tc_est] ──────
        lower_T2 = Tc_est - 0.15 * (Tc_est - T[0])
        mask2 = (T >= lower_T2) & (T < Tc_est) & (M_smooth > M_threshold)
        T_fit2 = T[mask2]
        M_fit2 = M_smooth[mask2]

        popt2, pcov2 = _run_fit(T_fit2, M_fit2)

        if popt2 is not None:
            A2, Tc2, beta2 = popt2

            if beta2 <= BETA_MAX_PHYSICAL:
                # Сужение помогло — принимаем новый результат
                print(
                    f"[free beta] Повторный фит (суженный диапазон) успешен: "
                    f"beta = {beta2:.4f}  (было {beta:.4f})"
                )
                A, Tc, beta, pcov = A2, Tc2, beta2, pcov2
                T_fit1, M_fit1 = T_fit2, M_fit2
                note = ""
            else:
                # beta по-прежнему > 0.5 — результат помечаем как ненадёжный
                reliable = False
                note = (
                    f"⚠ beta = {beta2:.3f} > 0.5 даже в суженном диапазоне.\n"
                    f"Результат ненадёжен. Возможны конечно-размерные эффекты\n"
                    f"или недостаточно точек вблизи T_C."
                )
                warnings.warn(
                    f"curie_fit_free_beta: beta = {beta2:.4f} > {BETA_MAX_PHYSICAL} "
                    f"даже после сужения диапазона. "
                    f"Результат возвращается с флагом reliable=False.",
                    UserWarning,
                )
                A, Tc, beta, pcov = A2, Tc2, beta2, pcov2
                T_fit1, M_fit1 = T_fit2, M_fit2
        else:
            # Повторный фит не сошёлся — оставляем первый результат, ненадёжный
            reliable = False
            note = (
                f"⚠ beta = {beta:.3f} > 0.5, повторный фит не сошёлся.\n"
                f"Результат ненадёжен."
            )
            warnings.warn(
                "curie_fit_free_beta: повторный фит в суженном диапазоне "
                "не сошёлся. Используется первый результат (ненадёжный).",
                UserWarning,
            )

    # ── Погрешности ──────────────────────────────────────────────────────────
    try:
        perr = np.sqrt(np.diag(pcov))
        Tc_err, beta_err = perr[1], perr[2]
    except Exception:
        Tc_err = beta_err = float("nan")

    status = "OK" if reliable else "НЕНАДЁЖНО"
    print(f"[free beta | {status}] T_C  = {Tc:.2f} ± {Tc_err:.2f} K")
    print(f"[fixed beta] Tc_MFA = {Tc_mfa:.2f} K")
    print(f"[free beta | {status}] beta = {beta:.4f} ± {beta_err:.4f}")
    print(f"[free beta | {status}] A    = {A:.4f}")
    if not reliable:
        print(f"  → beta > {BETA_MAX_PHYSICAL}: результат за пределами физически "
              f"допустимого диапазона. Рекомендуется finite-size scaling.")

    T_model = np.linspace(T_fit1.min(), Tc * 0.9999, 500)
    M_model = critical_law(T_model, A, Tc, beta)
    beta_label = f"{beta:.3f}"
    if not reliable:
        beta_label += " (!)"

    _save_plot(T, M, M_smooth, T_fit1, M_fit1, T_model, M_model,
               Tc, Tc_mfa, beta, beta_label, save_path, note=note)

    return Tc, beta, reliable