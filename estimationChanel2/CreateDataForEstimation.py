"""
ГЕНЕРАТОР ТЕСТОВЫХ ДАННЫХ ДЛЯ ОЦЕНКИ ПАРАМЕТРОВ КАНАЛА

ФОРМУЛЫ ИЗ СТАТЬИ:
- (1) Тестовый сигнал: S(t) = sum(a_m * cos(omega_m * t - phi_m))
- (3) Отклик канала: s_m(t) = sum(A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k))
- (4) Полезный сигнал: s(t) = sum(s_m(t))
- (2) Наблюдаемый сигнал: xi(t) = s(t) + n(t)
"""

import numpy as np


def generate_test_signal_and_channel_response(M, a_m, omega_m, phi_m, N, A_k, tau_k, psi_k, fs, T1, T2, SNR_dB=None):
    """
    Генерация тестового сигнала и отклика многолучевого канала

    МАТЕМАТИЧЕСКАЯ МОДЕЛЬ (из раздела 1 статьи):

    ФОРМУЛА (1): Тестовый сигнал
    S(t) = sum_{m=1}^{M} a_m * cos(omega_m * t - phi_m)

    ФОРМУЛА (3): Отклик канала на гармонику с номером m
    s_m(t) = sum_{k=1}^{N} A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k)

    ФОРМУЛА (4): Полезный сигнал на выходе канала
    s(t) = sum_{m=1}^{M} s_m(t) = sum_{m=1}^{M} sum_{k=1}^{N} A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k)

    ФОРМУЛА (2): Наблюдаемый сигнал
    xi(t) = s(t) + n(t)

    ПАРАМЕТРЫ:
    ----------
    M : int - количество синусоид в тестовом сигнале
    a_m : array (M,) - амплитуды синусоид
    omega_m : array (M,) - угловые частоты [рад/с]
    phi_m : array (M,) - начальные фазы [рад]
    N : int - количество лучей
    A_k : array (N,) - амплитуды лучей
    tau_k : array (N,) - задержки лучей [с]
    psi_k : array (N,) - фазовые сдвиги лучей [рад]
    fs : float - частота дискретизации [Гц]
    T1, T2 : float - интервал наблюдения [с]
    SNR_dB : float or None - ОСШ в дБ (None - без шума)

    ВОЗВРАЩАЕТ:
    -----------
    dict : с полями:
        - 't': временная сетка
        - 'S_test': тестовый сигнал S(t) (формула 1)
        - 's_useful': полезный сигнал s(t) (формула 4)
        - 'xi_observed': наблюдаемый сигнал xi(t) (формула 2)
        - 'params': параметры модели
        - 'SNR_dB': заданный SNR
    """

    # Временная сетка
    t = np.arange(T1, T2, 1 / fs)
    n_samples = len(t)

    # ========================================================================
    # ФОРМУЛА (1): Тестовый сигнал S(t) = sum(a_m * cos(omega_m * t - phi_m))
    # ========================================================================
    S_test = np.zeros(n_samples)
    for m in range(M):
        S_test += a_m[m] * np.cos(omega_m[m] * t - phi_m[m])

    # ========================================================================
    # ФОРМУЛА (4): Полезный сигнал s(t) = sum_m sum_k A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k)
    # ========================================================================
    s_useful = np.zeros(n_samples)
    for m in range(M):
        for k in range(N):
            # ФОРМУЛА (3): s_m(t) = sum_k A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k)
            # Аргумент косинуса: omega_m*(t - tau_k) - phi_m - psi_k
            s_useful += A_k[k] * a_m[m] * np.cos(omega_m[m] * (t - tau_k[k]) - phi_m[m] - psi_k[k])

    # ========================================================================
    # ФОРМУЛА (2): xi(t) = s(t) + n(t) - добавление шума
    # ========================================================================
    if SNR_dB is not None:
        Ps = np.mean(s_useful ** 2)
        Pn = Ps / (10 ** (SNR_dB / 10))
        noise = np.sqrt(Pn) * np.random.randn(n_samples)
        xi_observed = s_useful + noise
    else:
        xi_observed = s_useful.copy()

    # Результат
    result = {
        't': t,
        'S_test': S_test,  # Формула (1)
        's_useful': s_useful,  # Формула (4)
        'xi_observed': xi_observed,  # Формула (2)
        'params': {
            'M': M, 'N': N,
            'a_m': a_m, 'omega_m': omega_m, 'phi_m': phi_m,
            'A_k': A_k, 'tau_k': tau_k, 'psi_k': psi_k,
            'fs': fs, 'T1': T1, 'T2': T2
        },
        'SNR_dB': SNR_dB
    }

    return result


# ============================================================================
# ПРИМЕР ИСПОЛЬЗОВАНИЯ
# ============================================================================

if __name__ == "__main__":
    # ----- ПАРАМЕТРЫ -----

    # Тестовый сигнал (формула 1)
    M = 3
    f_m = np.array([100, 200, 300])  # Гц
    omega_m = 2 * np.pi * f_m  # рад/с
    a_m = np.array([1.0, 0.8, 0.6])
    phi_m = np.array([0, np.pi / 4, np.pi / 3])

    # Канал (лучи)
    N = 2
    A_k = np.array([0.7, 0.5])
    tau_k = np.array([500e-6, 1200e-6])  # 500 и 1200 мкс
    psi_k = np.array([0.52, -0.53])

    # Дискретизация
    fs = 20 * np.max(f_m)  # 6000 Гц
    T1 = 0.0
    T2 = 0.01  # 10 мс

    # Шум
    SNR_dB = 2

    # ----- ГЕНЕРАЦИЯ -----

    data = generate_test_signal_and_channel_response(
        M=M, a_m=a_m, omega_m=omega_m, phi_m=phi_m,
        N=N, A_k=A_k, tau_k=tau_k, psi_k=psi_k,
        fs=fs, T1=T1, T2=T2,
        SNR_dB=SNR_dB
    )

    print(f"\nРезультат:")
    print(f"  - t: {len(data['t'])} отсчетов")
    print(f"  - S_test: {data['S_test'].shape}")
    print(f"  - s_useful: {data['s_useful'].shape}")
    print(f"  - xi_observed: {data['xi_observed'].shape}")
    print(f"  - SNR: {data['SNR_dB']} дБ")