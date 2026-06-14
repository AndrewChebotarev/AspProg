"""
Monte Carlo моделирование для оценки параметров канала
Сравнение с границей Крамера-Рао (CRB) - ЛИНЕЙНЫЙ МАСШТАБ
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


def compute_matrix_Q(omega_m, a_m, tau, N0):
    M = len(omega_m)
    N = len(tau)
    Qc = np.zeros((N, N))
    Qcs = np.zeros((N, N))
    for i in range(N):
        for k in range(N):
            sum_cos = 0.0
            sum_sin = 0.0
            for m in range(M):
                weight = 2 * a_m[m] ** 2 / N0
                delta = tau[i] - tau[k]
                sum_cos += weight * np.cos(omega_m[m] * delta)
                sum_sin += weight * np.sin(omega_m[m] * (tau[k] - tau[i]))
            Qc[i, k] = sum_cos
            Qcs[i, k] = sum_sin
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
        return 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    except:
        return -np.inf


def estimate_amplitudes_and_phases(Z_prime, H, Q, T_dur):
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    N = len(A_cs) // 2
    A_c = A_cs[:N]
    A_s = A_cs[N:]
    A_est = np.sqrt(A_c ** 2 + A_s ** 2)
    psi_est = np.arctan2(A_s, A_c)
    return A_est, psi_est


def estimate_delay(Z_prime, omega_m, a_m, phi_m, tau_true, T_dur, N0):
    tau_range = np.linspace(tau_true - 100e-6, tau_true + 100e-6, 201)
    L_vals = []
    for tau_test in tau_range:
        tau = np.array([tau_test])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L_vals.append(compute_likelihood(Z_prime, H, Q))
    L_vals = np.array(L_vals)

    idx_max = np.argmax(L_vals)
    tau_est = tau_range[idx_max]

    if 0 < idx_max < len(L_vals) - 1:
        x = tau_range[idx_max-1:idx_max+2]
        y = L_vals[idx_max-1:idx_max+2]
        a = (y[0] - 2*y[1] + y[2]) / (2 * (x[0] - x[1])**2)
        b = (y[2] - y[0]) / (2 * (x[2] - x[0]))
        if a != 0:
            tau_est = x[1] - b / (2 * a)

    return tau_est


def compute_crb_parameters(f_m, a_m, SNR_power, T_dur, fs, A):
    N_samples = T_dur * fs
    M = len(f_m)
    omega_m = 2 * np.pi * f_m

    sum_a2 = np.sum(a_m**2)

    mean_omega = np.sum(omega_m * a_m**2) / sum_a2
    mean_omega2 = np.sum(omega_m**2 * a_m**2) / sum_a2
    sigma_omega2 = mean_omega2 - mean_omega**2

    crb_A_abs = A / np.sqrt(N_samples * SNR_power * sum_a2)
    crb_A_percent = crb_A_abs / A * 100

    crb_tau_sec = 1 / np.sqrt(N_samples * SNR_power * sum_a2 * sigma_omega2)
    crb_tau_us = crb_tau_sec * 1e6

    crb_psi_rad = np.sqrt(mean_omega2) / np.sqrt(N_samples * SNR_power * sum_a2 * sigma_omega2)

    return crb_A_percent, crb_tau_us, crb_psi_rad


if __name__ == "__main__":

    print("=" * 80)
    print("MONTE CARLO МОДЕЛИРОВАНИЕ - ЛИНЕЙНЫЙ МАСШТАБ")
    print("=" * 80)

    M = 2
    f_m = np.array([900, 1100])
    omega_m = 2 * np.pi * f_m
    a_m = np.array([1.0, 1.0])
    phi_m = np.array([0, 0])

    N_true = 1
    A_true = 0.7
    tau_true = 600e-6
    psi_true = 0.52

    fs = 20000
    T1 = 0.0
    T2 = 1.0
    T_dur = T2 - T1
    N0 = 1.0

    SNR_dB_list = [-30, -20, -10, -5, 0, 5, 10, 15, 20, 30, 40]
    n_iterations = 200

    tau_rmse_us = np.zeros(len(SNR_dB_list))
    A_rmse_percent = np.zeros(len(SNR_dB_list))
    psi_rmse_rad = np.zeros(len(SNR_dB_list))

    tau_bias_us = np.zeros(len(SNR_dB_list))
    A_bias_percent = np.zeros(len(SNR_dB_list))
    psi_bias_rad = np.zeros(len(SNR_dB_list))

    crb_tau_us = np.zeros(len(SNR_dB_list))
    crb_A_percent = np.zeros(len(SNR_dB_list))
    crb_psi_rad = np.zeros(len(SNR_dB_list))

    print(f"\nЗапуск Monte Carlo симуляции...")
    print(f"  Частоты: {f_m[0]} и {f_m[1]} Гц")
    print(f"  T_dur = {T_dur} с, fs = {fs} Гц")
    print(f"  N_samples = {T_dur * fs}")
    print(f"  SNR уровни: {SNR_dB_list}")
    print(f"  Итераций на SNR: {n_iterations}\n")

    for idx_snr, SNR_dB in enumerate(tqdm(SNR_dB_list, desc="Обработка SNR")):
        tau_errors_us = []
        A_errors_percent = []
        psi_errors_rad = []

        for iter_num in range(n_iterations):
            try:
                data = generate_test_signal_and_channel_response(
                    M, a_m, omega_m, phi_m,
                    N_true, np.array([A_true]), np.array([tau_true]), np.array([psi_true]),
                    fs, T1, T2, SNR_dB
                )

                Z_prime, x_m, y_m = compute_sufficient_statistics(
                    data['xi_observed'], data['t'], omega_m, a_m, N0, fs
                )

                tau_est_sec = estimate_delay(Z_prime, omega_m, a_m, phi_m, tau_true, T_dur, N0)

                tau_arr = np.array([tau_est_sec])
                C, S = compute_matrices_C_S(omega_m, phi_m, tau_arr)
                H_opt = compute_matrix_H(C, S)
                Q_opt = compute_matrix_Q(omega_m, a_m, tau_arr, N0)
                A_est, psi_est = estimate_amplitudes_and_phases(Z_prime, H_opt, Q_opt, T_dur)

                tau_err_us = (tau_est_sec - tau_true) * 1e6
                A_err_percent = (A_est[0] - A_true) / A_true * 100
                psi_err_rad = psi_est[0] - psi_true

                tau_errors_us.append(tau_err_us)
                A_errors_percent.append(A_err_percent)
                psi_errors_rad.append(psi_err_rad)

            except Exception as e:
                pass

        tau_errors_us = [e for e in tau_errors_us if not np.isnan(e)]
        A_errors_percent = [e for e in A_errors_percent if not np.isnan(e)]
        psi_errors_rad = [e for e in psi_errors_rad if not np.isnan(e)]

        if len(tau_errors_us) > 0:
            tau_rmse_us[idx_snr] = np.sqrt(np.mean(np.array(tau_errors_us)**2))
            tau_bias_us[idx_snr] = np.mean(tau_errors_us)
        if len(A_errors_percent) > 0:
            A_rmse_percent[idx_snr] = np.sqrt(np.mean(np.array(A_errors_percent)**2))
            A_bias_percent[idx_snr] = np.mean(A_errors_percent)
        if len(psi_errors_rad) > 0:
            psi_rmse_rad[idx_snr] = np.sqrt(np.mean(np.array(psi_errors_rad)**2))
            psi_bias_rad[idx_snr] = np.mean(psi_errors_rad)

        SNR_power = 10**(SNR_dB/10)
        crb_A_percent[idx_snr], crb_tau_us[idx_snr], crb_psi_rad[idx_snr] = \
            compute_crb_parameters(f_m, a_m, SNR_power, T_dur, fs, A_true)

    # ========== ВЫВОД РЕЗУЛЬТАТОВ ==========
    print("\n" + "=" * 80)
    print("РЕЗУЛЬТАТЫ")
    print("=" * 80)

    print(f"\n{'SNR (дБ)':>10} | {'τ RMSE (мкс)':>15} | {'A RMSE (%)':>15} | {'ψ RMSE (рад)':>15}")
    print("-" * 65)
    for idx, SNR_dB in enumerate(SNR_dB_list):
        print(f"{SNR_dB:10.0f} | {tau_rmse_us[idx]:15.4f} | {A_rmse_percent[idx]:15.4f} | {psi_rmse_rad[idx]:15.4f}")

    print("\n" + "=" * 80)
    print("ГРАНИЦА КРАМЕРА-РАО (CRB)")
    print("=" * 80)

    print(f"\n{'SNR (дБ)':>10} | {'CRB τ (мкс)':>18} | {'CRB A (%)':>15} | {'CRB ψ (рад)':>15}")
    print("-" * 65)
    for idx, SNR_dB in enumerate(SNR_dB_list):
        print(f"{SNR_dB:10.0f} | {crb_tau_us[idx]:18.6f} | {crb_A_percent[idx]:15.6f} | {crb_psi_rad[idx]:15.6f}")

    # ========== ГРАФИКИ В ЛИНЕЙНОМ МАСШТАБЕ (ТОЛЬКО ТРИ ГРАФИКА) ==========
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))  # 1 строка, 3 столбца, большой размер

    # 1. Задержка τ (ЛИНЕЙНЫЙ масштаб)
    ax = axes[0]
    ax.plot(SNR_dB_list, tau_rmse_us, 'bo-', linewidth=2, markersize=8, label='RMSE (МП оценка)')
    ax.plot(SNR_dB_list, crb_tau_us, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ', fontsize=12)
    ax.set_ylabel('RMSE τ, мкс', fontsize=12)
    ax.set_title('Ошибка оценки задержки τ', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    # 2. Амплитуда A (ЛИНЕЙНЫЙ масштаб)
    ax = axes[1]
    ax.plot(SNR_dB_list, A_rmse_percent, 'bo-', linewidth=2, markersize=8, label='RMSE (МП оценка)')
    ax.plot(SNR_dB_list, crb_A_percent, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ', fontsize=12)
    ax.set_ylabel('RMSE A, %', fontsize=12)
    ax.set_title('Ошибка оценки амплитуды A', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    # 3. Фаза ψ (ЛИНЕЙНЫЙ масштаб)
    ax = axes[2]
    ax.plot(SNR_dB_list, psi_rmse_rad, 'bo-', linewidth=2, markersize=8, label='RMSE (МП оценка)')
    ax.plot(SNR_dB_list, crb_psi_rad, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ', fontsize=12)
    ax.set_ylabel('RMSE ψ, рад', fontsize=12)
    ax.set_title('Ошибка оценки фазы ψ', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    plt.suptitle('Monte Carlo анализ оценки параметров канала (линейный масштаб)', fontsize=16)
    plt.tight_layout()
    plt.show()

    print("\n" + "=" * 80)
    print("МОДЕЛИРОВАНИЕ ЗАВЕРШЕНО")
    print("=" * 80)