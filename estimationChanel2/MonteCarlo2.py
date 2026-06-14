"""
Monte Carlo моделирование для оценки параметров канала (ДВА ЛУЧА)
Сравнение с границей Крамера-Рао (CRB) - ИСПРАВЛЕННАЯ ВЕРСИЯ
"""

import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')


def compute_matrices_C_S(omega_m, phi_m, tau):
    M = len(omega_m)
    N = len(tau)
    C = np.zeros((M, N))
    S = np.zeros((M, N))
    for m in range(M):
        for k in range(N):
            phase = omega_m[m] * tau[k] + phi_m[m]
            C[m, k] = np.cos(phase)
            S[m, k] = np.sin(phase)
    return C, S


def compute_matrix_H(C, S):
    M, N = C.shape
    H = np.zeros((2 * M, 2 * N))
    H[:M, :N] = C
    H[:M, N:] = -S
    H[M:, :N] = S
    H[M:, N:] = C
    return H


# ИСПРАВЛЕНО: Добавлен аргумент T_dur и множитель coef = T_dur / 2
def compute_matrix_Q(omega_m, a_m, tau, N0, T_dur):
    M = len(omega_m)
    N = len(tau)
    Qc = np.zeros((N, N))
    Qcs = np.zeros((N, N))
    coef = T_dur / 2  # Множитель из формул (80)-(83)

    for i in range(N):
        for k in range(N):
            sum_cos = 0.0
            sum_sin = 0.0
            for m in range(M):
                weight = 2 * a_m[m] ** 2 / N0
                delta = tau[i] - tau[k]
                sum_cos += weight * np.cos(omega_m[m] * delta)
                sum_sin += weight * np.sin(omega_m[m] * (tau[k] - tau[i]))
            Qc[i, k] = coef * sum_cos
            Qcs[i, k] = coef * sum_sin

    Q = np.vstack([np.hstack([Qc, Qcs]), np.hstack([Qcs.T, Qc])])
    return Q


def compute_sufficient_statistics(xi, t, omega_m, a_m, N0, fs):
    M = len(omega_m)
    dt = 1 / fs
    x_m = np.zeros(M)
    y_m = np.zeros(M)
    for m in range(M):
        cos_sum = np.sum(xi * np.cos(omega_m[m] * t)) * dt
        sin_sum = np.sum(xi * np.sin(omega_m[m] * t)) * dt
        x_m[m] = (2.0 / N0) * cos_sum
        y_m[m] = (2.0 / N0) * sin_sum
    Z_prime = np.concatenate([a_m * x_m, a_m * y_m])
    return Z_prime, x_m, y_m


def compute_likelihood(Z_prime, H, Q):
    try:
        Q_inv = np.linalg.inv(Q)
        return 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime
    except np.linalg.LinAlgError:
        return -np.inf


# ИСПРАВЛЕНО: Убрана искусственная нормировка
def estimate_amplitudes_and_phases(Z_prime, H, Q):
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    N = len(A_cs) // 2
    A_c = A_cs[:N]
    A_s = A_cs[N:]
    A_est = np.sqrt(A_c ** 2 + A_s ** 2)
    psi_est = np.arctan2(A_s, A_c)
    return A_est, psi_est


def estimate_delays_2d(Z_prime, omega_m, a_m, phi_m, tau1_range, tau2_range, T_dur, N0):
    """Двумерный поиск задержек для двух лучей"""
    n1, n2 = len(tau1_range), len(tau2_range)
    L_grid = np.zeros((n1, n2))

    for i, t1 in enumerate(tau1_range):
        for j, t2 in enumerate(tau2_range):
            tau = np.array([t1, t2])
            C, S = compute_matrices_C_S(omega_m, phi_m, tau)
            H = compute_matrix_H(C, S)
            Q = compute_matrix_Q(omega_m, a_m, tau, N0, T_dur)
            L_grid[i, j] = compute_likelihood(Z_prime, H, Q)

    idx = np.unravel_index(np.argmax(L_grid), L_grid.shape)
    tau_est = np.array([tau1_range[idx[0]], tau2_range[idx[1]]])
    return tau_est, L_grid


# ИСПРАВЛЕНО: Обновленная функция CRB для двух лучей
def compute_crb_parameters_two_paths(f_m, a_m, SNR_power, T_dur, fs, A_true):
    """
    Вычисляет средние CRB для двух лучей.
    """
    N_samples = T_dur * fs
    omega_m = 2 * np.pi * f_m

    sum_a2 = np.sum(a_m**2)
    mean_omega2 = np.sum(omega_m**2 * a_m**2) / sum_a2
    mean_omega = np.sum(omega_m * a_m**2) / sum_a2
    sigma_omega2 = mean_omega2 - mean_omega**2

    base_factor = 1.0 / (N_samples * SNR_power)

    # Средняя CRB по амплитуде
    crb_A_avg = np.mean([A / np.sqrt(sum_a2) for A in A_true])
    crb_A_percent = crb_A_avg / np.mean(A_true) * 100

    # Средняя CRB по задержке
    crb_tau_sec = 1 / np.sqrt(sum_a2 * sigma_omega2)
    crb_tau_us = crb_tau_sec * 1e6

    # Средняя CRB по фазе
    crb_psi_rad = np.sqrt(mean_omega2) / np.sqrt(sum_a2 * sigma_omega2)

    # Применяем общий множитель
    crb_A_percent *= np.sqrt(base_factor)
    crb_tau_us *= np.sqrt(base_factor)
    crb_psi_rad *= np.sqrt(base_factor)

    return crb_A_percent, crb_tau_us, crb_psi_rad


if __name__ == "__main__":

    print("=" * 80)
    print("MONTE CARLO МОДЕЛИРОВАНИЕ - ДВА ЛУЧА (ИСПРАВЛЕНО)")
    print("=" * 80)

    # ========== ПАРАМЕТРЫ ДЛЯ ДВУХ ЛУЧЕЙ ==========
    # Увеличиваем полосу частот для разрешения лучей!
    M = 6
    f_m = np.array([100, 2000, 4000, 6000, 8000, 10000])  # Полоса ~9900 Гц
    omega_m = 2 * np.pi * f_m
    a_m = np.ones(M)
    phi_m = np.zeros(M)

    N_true = 2
    A_true = np.array([0.7, 0.5])
    tau_true = np.array([400e-6, 1500e-6])  # Разность = 1100 мкс >> 400 мкс
    psi_true = np.array([0.52, -0.53])

    fs = 10000  # 10 кГц достаточно для 3000 Гц
    T1 = 0.0
    T2 = 1    # 100 мс
    T_dur = T2 - T1
    N0 = 1.0

    SNR_dB_list = [0, 5, 10, 15, 20, 25, 30]
    n_iterations = 100  # Меньше итераций из-за 2D поиска

    # Массивы для результатов (по каждому лучу отдельно)
    tau1_rmse_us = np.zeros(len(SNR_dB_list))
    tau2_rmse_us = np.zeros(len(SNR_dB_list))
    A1_rmse_percent = np.zeros(len(SNR_dB_list))
    A2_rmse_percent = np.zeros(len(SNR_dB_list))
    psi1_rmse_rad = np.zeros(len(SNR_dB_list))
    psi2_rmse_rad = np.zeros(len(SNR_dB_list))

    crb_tau_us = np.zeros(len(SNR_dB_list))
    crb_A_percent = np.zeros(len(SNR_dB_list))
    crb_psi_rad = np.zeros(len(SNR_dB_list))

    print(f"\nЗапуск Monte Carlo симуляции...")
    print(f"  Частоты: {f_m} Гц")
    print(f"  T_dur = {T_dur} с, fs = {fs} Гц")
    print(f"  N_samples = {T_dur * fs}")
    print(f"  Истинные τ₁ = {tau_true[0]*1e6:.0f} мкс, τ₂ = {tau_true[1]*1e6:.0f} мкс")
    print(f"  SNR уровни: {SNR_dB_list}")
    print(f"  Итераций на SNR: {n_iterations}\n")

    for idx_snr, SNR_dB in enumerate(tqdm(SNR_dB_list, desc="Обработка SNR")):
        tau1_errors_us = []
        tau2_errors_us = []
        A1_errors_percent = []
        A2_errors_percent = []
        psi1_errors_rad = []
        psi2_errors_rad = []

        # Сетка поиска задержек (вокруг истинных значений)
        tau1_range = np.linspace(tau_true[0] - 300e-6, tau_true[0] + 300e-6, 31)
        tau2_range = np.linspace(tau_true[1] - 300e-6, tau_true[1] + 300e-6, 31)

        for iter_num in range(n_iterations):
            try:
                data = generate_test_signal_and_channel_response(
                    M, a_m, omega_m, phi_m,
                    N_true, A_true, tau_true, psi_true,
                    fs, T1, T2, SNR_dB
                )

                Z_prime, x_m, y_m = compute_sufficient_statistics(
                    data['xi_observed'], data['t'], omega_m, a_m, N0, fs
                )

                tau_est, L_grid = estimate_delays_2d(
                    Z_prime, omega_m, a_m, phi_m, tau1_range, tau2_range, T_dur, N0
                )

                # Оценка амплитуд и фаз при найденных задержках
                C, S = compute_matrices_C_S(omega_m, phi_m, tau_est)
                H_opt = compute_matrix_H(C, S)
                Q_opt = compute_matrix_Q(omega_m, a_m, tau_est, N0, T_dur)
                A_est, psi_est = estimate_amplitudes_and_phases(Z_prime, H_opt, Q_opt)

                # Сортируем лучи по задержкам для корректного сравнения
                if tau_est[0] > tau_est[1]:
                    tau_est = tau_est[::-1]
                    A_est = A_est[::-1]
                    psi_est = psi_est[::-1]

                tau1_err_us = (tau_est[0] - tau_true[0]) * 1e6
                tau2_err_us = (tau_est[1] - tau_true[1]) * 1e6
                A1_err_percent = (A_est[0] - A_true[0]) / A_true[0] * 100
                A2_err_percent = (A_est[1] - A_true[1]) / A_true[1] * 100
                psi1_err_rad = psi_est[0] - psi_true[0]
                psi2_err_rad = psi_est[1] - psi_true[1]

                # Фильтрация грубых ошибок (лучи не перепутались)
                if abs(tau1_err_us) < 500 and abs(tau2_err_us) < 500:
                    tau1_errors_us.append(tau1_err_us)
                    tau2_errors_us.append(tau2_err_us)
                    A1_errors_percent.append(A1_err_percent)
                    A2_errors_percent.append(A2_err_percent)
                    psi1_errors_rad.append(psi1_err_rad)
                    psi2_errors_rad.append(psi2_err_rad)

            except Exception as e:
                pass

        # Сохранение RMSE
        if len(tau1_errors_us) > 0:
            tau1_rmse_us[idx_snr] = np.sqrt(np.mean(np.array(tau1_errors_us)**2))
            tau2_rmse_us[idx_snr] = np.sqrt(np.mean(np.array(tau2_errors_us)**2))
            A1_rmse_percent[idx_snr] = np.sqrt(np.mean(np.array(A1_errors_percent)**2))
            A2_rmse_percent[idx_snr] = np.sqrt(np.mean(np.array(A2_errors_percent)**2))
            psi1_rmse_rad[idx_snr] = np.sqrt(np.mean(np.array(psi1_errors_rad)**2))
            psi2_rmse_rad[idx_snr] = np.sqrt(np.mean(np.array(psi2_errors_rad)**2))

        # CRB для двух лучей
        SNR_power = 10**(SNR_dB/10)
        crb_A_percent[idx_snr], crb_tau_us[idx_snr], crb_psi_rad[idx_snr] = \
            compute_crb_parameters_two_paths(f_m, a_m, SNR_power, T_dur, fs, A_true)

    # ========== ВЫВОД РЕЗУЛЬТАТОВ ==========
    print("\n" + "=" * 80)
    print("РЕЗУЛЬТАТЫ MONTE CARLO (ДВА ЛУЧА)")
    print("=" * 80)

    print(f"\n{'SNR (дБ)':>10} | {'τ₁ RMSE (мкс)':>15} | {'τ₂ RMSE (мкс)':>15} | {'A₁ RMSE (%)':>15} | {'A₂ RMSE (%)':>15}")
    print("-" * 85)
    for idx, SNR_dB in enumerate(SNR_dB_list):
        print(f"{SNR_dB:10.0f} | {tau1_rmse_us[idx]:15.4f} | {tau2_rmse_us[idx]:15.4f} | {A1_rmse_percent[idx]:15.4f} | {A2_rmse_percent[idx]:15.4f}")

    print("\n" + "=" * 80)
    print("ГРАНИЦА КРАМЕРА-РАО (CRB) ДЛЯ ДВУХ ЛУЧЕЙ")
    print("=" * 80)

    print(f"\n{'SNR (дБ)':>10} | {'CRB τ (мкс)':>18} | {'CRB A (%)':>15} | {'CRB ψ (рад)':>15}")
    print("-" * 65)
    for idx, SNR_dB in enumerate(SNR_dB_list):
        print(f"{SNR_dB:10.0f} | {crb_tau_us[idx]:18.6f} | {crb_A_percent[idx]:15.6f} | {crb_psi_rad[idx]:15.6f}")

    # ========== ГРАФИКИ ==========
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Задержка τ₁ и τ₂
    ax = axes[0, 0]
    ax.plot(SNR_dB_list, tau1_rmse_us, 'bo-', linewidth=2, markersize=8, label='RMSE τ₁')
    ax.plot(SNR_dB_list, tau2_rmse_us, 'gs-', linewidth=2, markersize=8, label='RMSE τ₂')
    ax.plot(SNR_dB_list, crb_tau_us, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('RMSE τ, мкс')
    ax.set_title('Ошибка оценки задержек')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. Амплитуда A₁ и A₂
    ax = axes[0, 1]
    ax.plot(SNR_dB_list, A1_rmse_percent, 'bo-', linewidth=2, markersize=8, label='RMSE A₁')
    ax.plot(SNR_dB_list, A2_rmse_percent, 'gs-', linewidth=2, markersize=8, label='RMSE A₂')
    ax.plot(SNR_dB_list, crb_A_percent, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('RMSE A, %')
    ax.set_title('Ошибка оценки амплитуд')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. Фаза ψ₁ и ψ₂
    ax = axes[1, 0]
    ax.plot(SNR_dB_list, psi1_rmse_rad, 'bo-', linewidth=2, markersize=8, label='RMSE ψ₁')
    ax.plot(SNR_dB_list, psi2_rmse_rad, 'gs-', linewidth=2, markersize=8, label='RMSE ψ₂')
    ax.plot(SNR_dB_list, crb_psi_rad, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('RMSE ψ, рад')
    ax.set_title('Ошибка оценки фаз')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. Отношение RMSE/CRB (усредненное)
    ax = axes[1, 1]
    ratio1 = tau1_rmse_us / crb_tau_us
    ratio2 = tau2_rmse_us / crb_tau_us
    ax.plot(SNR_dB_list, ratio1, 'bo-', linewidth=2, markersize=8, label='τ₁/CRB')
    ax.plot(SNR_dB_list, ratio2, 'gs-', linewidth=2, markersize=8, label='τ₂/CRB')
    ax.axhline(1, color='r', linestyle='--', alpha=0.7, label='CRB (оптимум)')
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('RMSE / CRB')
    ax.set_title('Относительная эффективность оценки τ')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle('Monte Carlo анализ оценки параметров канала (ДВА ЛУЧА)', fontsize=14)
    plt.tight_layout()
    plt.show()

    print("\n" + "=" * 80)
    print("МОДЕЛИРОВАНИЕ ЗАВЕРШЕНО")
    print("=" * 80)