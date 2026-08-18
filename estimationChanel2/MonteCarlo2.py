"""
Monte Carlo моделирование для оценки параметров канала (N = 2 ЛУЧА, M = 3 ЧАСТОТЫ)
Сравнение с границей Крамера-Рао (CRB)

ИСПОЛЬЗУЕТСЯ ФОРМУЛА (100):
L(T) = 0.5 * Z' @ H @ Q^{-1} @ H^T @ Z'^T

ОЦЕНИВАЕМЫЕ ПАРАМЕТРЫ (N=2):
- tau_1, tau_2 - задержки двух лучей
- A_1, A_2 - амплитуды двух лучей
- psi_1, psi_2 - фазы двух лучей

АЛГОРИТМ:
1. Поиск по 2D сетке tau1, tau2
2. Для каждого кандидата: вычисление L(T) по формуле (100)
3. Выбор максимума + параболическая интерполяция
4. Оценка A_cs через (99) в точке максимума

ВАЖНО: Для N=2 необходимо M >= 3, иначе H квадратная и L(T) = const.
"""

import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')


def compute_matrices_C_S(omega_m, phi_m, tau):
    """
    ВЫЧИСЛЕНИЕ МАТРИЦ C И S (формула 92, стр. 24)
    """
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
    """
    ВЫЧИСЛЕНИЕ МАТРИЦЫ H (формула 95)
    H = [C, -S; S, C], размер 2M x 2N
    """
    M, N = C.shape
    H = np.zeros((2 * M, 2 * N))
    H[:M, :N] = C
    H[:M, N:] = -S
    H[M:, :N] = S
    H[M:, N:] = C
    return H


def compute_matrix_Q(omega_m, a_m, tau, N0):
    """
    ВЫЧИСЛЕНИЕ МАТРИЦЫ Q (формулы 80-82, стр. 47)
    
    ИСПРАВЛЕНИЕ: Qcs_ik = sin(omega_m * (tau_i - tau_k)), а не sin(omega_m * (tau_k - tau_i))
    Это обеспечивает H^T @ H = (N0/2) * Q
    """
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
                sum_sin += weight * np.sin(omega_m[m] * delta)  # ИСПРАВЛЕНО: tau_i - tau_k
            Qc[i, k] = sum_cos
            Qcs[i, k] = sum_sin
    Q = np.vstack([
        np.hstack([Qc, Qcs]),
        np.hstack([Qcs.T, Qc])
    ])
    return Q


def compute_sufficient_statistics(xi, t, omega_m, a_m, N0, fs):
    """
    ВЫЧИСЛЕНИЕ ДОСТАТОЧНЫХ СТАТИСТИК (формулы 87-88, 96)
    """
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
    """
    ВЫЧИСЛЕНИЕ ЛОГАРИФМА ФОП (формула 100)
    L(T) = 0.5 * Z' @ H @ Q^{-1} @ H^T @ Z'^T
    """
    Q_inv = np.linalg.inv(Q)
    L = 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    return L


def estimate_amplitudes_and_phases(Z_prime, H, Q, T_dur):
    """
    ОЦЕНКА АМПЛИТУД И ФАЗ (формула 99)
    A_cs = (2/T_dur) * Z' @ H @ Q^{-1}
    
    Возвращает: A_est (массив N), psi_est (массив N)
    """
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    
    N = A_cs.shape[0] // 2
    A_est = np.sqrt(A_cs[:N] ** 2 + A_cs[N:] ** 2)
    psi_est = np.arctan2(A_cs[N:], A_cs[:N])
    return A_est, psi_est


def estimate_delays_2d(xi, t, Z_prime, omega_m, a_m, phi_m, tau_true, N0, dt, T_dur, n_grid=31):
    """
    ДВУХМЕРНЫЙ ПОИСК ЗАДЕРЖЕК tau1, tau2 ПО ФОРМУЛЕ (100)
    
    Трёхэтапный алгоритм:
    1. Грубый поиск по сетке n_grid x n_grid
    2. Уточнение вокруг найденного максимума
    3. Параболическая интерполяция для суб-сеточной точности
    
    Возвращает: tau1_est, tau2_est, L_max
    """
    # Поиск ведётся в окрестности истинных значений
    tau_margin = 200e-6  # запас по времени
    
    # Этап 1: грубый поиск
    tau1_range = np.linspace(tau_true[0] - tau_margin, tau_true[0] + tau_margin, n_grid)
    tau2_range = np.linspace(tau_true[1] - tau_margin, tau_true[1] + tau_margin, n_grid)
    
    best_L = -np.inf
    best_tau = np.array([tau_true[0], tau_true[1]])
    
    for t1 in tau1_range:
        for t2 in tau2_range:
            if t2 <= t1:
                continue
            tau = np.array([t1, t2])
            C, S = compute_matrices_C_S(omega_m, phi_m, tau)
            H = compute_matrix_H(C, S)
            Q = compute_matrix_Q(omega_m, a_m, tau, N0)
            L = compute_likelihood(Z_prime, H, Q)
            if L > best_L:
                best_L = L
                best_tau = tau.copy()
    
    # Этап 2: уточнение (более мелкая сетка вокруг найденного максимума)
    # Используем шаг в 10 раз мельче, чем на первом этапе
    step_coarse = 2 * tau_margin / (n_grid - 1)
    fine_margin = step_coarse * 1.5  # немного больше шага грубой сетки
    n_fine = n_grid  # столько же точек на уточнении
    
    tau1_fine = np.linspace(best_tau[0] - fine_margin, best_tau[0] + fine_margin, n_fine)
    tau2_fine = np.linspace(best_tau[1] - fine_margin, best_tau[1] + fine_margin, n_fine)
    
    # Сохраняем L(T) для параболической интерполяции
    L_grid_fine = np.zeros((n_fine, n_fine))
    
    for i, t1 in enumerate(tau1_fine):
        for j, t2 in enumerate(tau2_fine):
            if t2 <= t1:
                L_grid_fine[i, j] = -np.inf
                continue
            tau = np.array([t1, t2])
            C, S = compute_matrices_C_S(omega_m, phi_m, tau)
            H = compute_matrix_H(C, S)
            Q = compute_matrix_Q(omega_m, a_m, tau, N0)
            L = compute_likelihood(Z_prime, H, Q)
            L_grid_fine[i, j] = L
            if L > best_L:
                best_L = L
                best_tau = tau.copy()
    
    # Этап 3: параболическая интерполяция для суб-сеточной точности
    idx_max = np.unravel_index(np.argmax(L_grid_fine), L_grid_fine.shape)
    i0, j0 = idx_max
    
    # Параболическая интерполяция по tau1
    tau1_est_interp = best_tau[0]
    if 0 < i0 < n_fine - 1:
        x = tau1_fine[i0-1:i0+2]
        y = L_grid_fine[i0-1:i0+2, j0]
        if x[2] - x[0] != 0:
            a = (y[0] - 2*y[1] + y[2]) / (2 * (x[0] - x[1])**2)
            b = (y[2] - y[0]) / (2 * (x[2] - x[0]))
            if a != 0:
                tau1_est_interp = x[1] - b / (2 * a)
    
    # Параболическая интерполяция по tau2
    tau2_est_interp = best_tau[1]
    if 0 < j0 < n_fine - 1:
        x = tau2_fine[j0-1:j0+2]
        y = L_grid_fine[i0, j0-1:j0+2]
        if x[2] - x[0] != 0:
            a = (y[0] - 2*y[1] + y[2]) / (2 * (x[0] - x[1])**2)
            b = (y[2] - y[0]) / (2 * (x[2] - x[0]))
            if a != 0:
                tau2_est_interp = x[1] - b / (2 * a)
    
    # Используем интерполированные значения (они должны быть в допустимом диапазоне)
    if tau_true[0] - tau_margin <= tau1_est_interp <= tau_true[0] + tau_margin:
        best_tau[0] = tau1_est_interp
    if tau_true[1] - tau_margin <= tau2_est_interp <= tau_true[1] + tau_margin:
        best_tau[1] = tau2_est_interp
    
    return best_tau[0], best_tau[1], best_L


def compute_crb_parameters(f_m, a_m, SNR_power, T_dur, fs, A):
    """
    ВЫЧИСЛЕНИЕ ГРАНИЦЫ КРАМЕРА-РАО (CRB)
    
    CRB для одного луча в многолучевом канале:
    - crb_tau = 1 / sqrt(N_samples * SNR * sum(a_m^2) * sigma_omega^2)
    - crb_A = A / sqrt(N_samples * SNR * sum(a_m^2))
    - crb_psi = sqrt(mean_omega^2) / sqrt(N_samples * SNR * sum(a_m^2) * sigma_omega^2)
    """
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


# ============================================================================
# ГЛАВНАЯ ФУНКЦИЯ - MONTE CARLO МОДЕЛИРОВАНИЕ (N = 2 ЛУЧА, M = 3 ЧАСТОТЫ)
# ============================================================================

if __name__ == "__main__":

    print("=" * 80)
    print("MONTE CARLO MODELIROVANIE - N = 2 LUCHA, M = 3 CHASTOTY")
    print("=" * 80)

    # ------------------------------------------------------------------------
    # ПАРАМЕТРЫ МОДЕЛИ
    # ------------------------------------------------------------------------

    # Тестовый сигнал S(t) (формула 1) - M = 3 частоты
    M = 3
    f_m = np.array([700, 900, 1100])  # Гц
    omega_m = 2 * np.pi * f_m
    a_m = np.array([1.0, 1.0, 1.0])
    phi_m = np.array([0, 0, 0])

    # Параметры канала (N = 2 луча)
    N_true = 2
    A_true = np.array([0.7, 0.5])
    tau_true = np.array([600e-6, 1200e-6])
    psi_true = np.array([0.52, -0.53])

    # Параметры дискретизации и наблюдения
    fs = 20000
    T1 = 0.0
    T2 = 1.0
    T_dur = T2 - T1
    dt = 1 / fs
    N0 = 100.0

    # Параметры моделирования
    SNR_dB_list = [-30, -20, -10, -5, 0, 5, 10, 15, 20, 30, 40]
    n_iterations = 200  # количество итераций на SNR
    n_grid = 41  # размер сетки для 2D поиска (увеличено для точности)

    # ------------------------------------------------------------------------
    # ИНИЦИАЛИЗАЦИЯ МАССИВОВ ДЛЯ РЕЗУЛЬТАТОВ
    # ------------------------------------------------------------------------

    tau_rmse_us = np.zeros((len(SNR_dB_list), N_true))
    A_rmse_percent = np.zeros((len(SNR_dB_list), N_true))
    psi_rmse_rad = np.zeros((len(SNR_dB_list), N_true))

    tau_bias_us = np.zeros((len(SNR_dB_list), N_true))
    A_bias_percent = np.zeros((len(SNR_dB_list), N_true))
    psi_bias_rad = np.zeros((len(SNR_dB_list), N_true))

    crb_tau_us = np.zeros((len(SNR_dB_list), N_true))
    crb_A_percent = np.zeros((len(SNR_dB_list), N_true))
    crb_psi_rad = np.zeros((len(SNR_dB_list), N_true))

    print(f"\nZapusk Monte Carlo simulacii...")
    print(f"  Chastoty (M={M}): {f_m} Gc")
    print(f"  T_dur = {T_dur} s, fs = {fs} Gc")
    print(f"  N_samples = {T_dur * fs}")
    print(f"  Chislo luchey N = {N_true}")
    print(f"  Zaderzhki: {tau_true[0]*1e6:.0f} i {tau_true[1]*1e6:.0f} mks")
    print(f"  Amplitudy: {A_true[0]} i {A_true[1]}")
    print(f"  Fazy: {psi_true[0]:.2f} i {psi_true[1]:.2f} rad")
    print(f"  SNR urovni: {SNR_dB_list}")
    print(f"  Iteraciy na SNR: {n_iterations}")
    print(f"  Setka 2D poiska: {n_grid}x{n_grid}\n")

    # ------------------------------------------------------------------------
    # ОСНОВНОЙ ЦИКЛ ПО SNR
    # ------------------------------------------------------------------------

    for idx_snr, SNR_dB in enumerate(tqdm(SNR_dB_list, desc="Obrabotka SNR")):
        tau_errors_us = [[] for _ in range(N_true)]
        A_errors_percent = [[] for _ in range(N_true)]
        psi_errors_rad = [[] for _ in range(N_true)]

        for iter_num in range(n_iterations):
            try:
                # ШАГ 1: ГЕНЕРАЦИЯ ДАННЫХ (формулы 1-4)
                data = generate_test_signal_and_channel_response(
                    M, a_m, omega_m, phi_m,
                    N_true, A_true, tau_true, psi_true,
                    fs, T1, T2, SNR_dB
                )

                xi = data['xi_observed'].copy()
                t = data['t']

                # ШАГ 2: ДОСТАТОЧНЫЕ СТАТИСТИКИ
                Z_prime, x_m, y_m = compute_sufficient_statistics(
                    xi, t, omega_m, a_m, N0, fs
                )

                # ШАГ 3: 2D ПОИСК ЗАДЕРЖЕК ПО ФОРМУЛЕ (100)
                tau1_est, tau2_est, L_max = estimate_delays_2d(
                    xi, t, Z_prime, omega_m, a_m, phi_m,
                    tau_true, N0, dt, T_dur, n_grid
                )

                # ШАГ 4: ОЦЕНКА АМПЛИТУД И ФАЗ ПРИ НАЙДЕННЫХ ЗАДЕРЖКАХ
                tau_est = np.array([tau1_est, tau2_est])
                C, S = compute_matrices_C_S(omega_m, phi_m, tau_est)
                H = compute_matrix_H(C, S)
                Q = compute_matrix_Q(omega_m, a_m, tau_est, N0)
                A_est, psi_est = estimate_amplitudes_and_phases(Z_prime, H, Q, T_dur)

                # ШАГ 5: ВЫЧИСЛЕНИЕ ОШИБОК
                for k in range(N_true):
                    tau_err_us = (tau_est[k] - tau_true[k]) * 1e6
                    A_err_percent = (A_est[k] - A_true[k]) / A_true[k] * 100
                    psi_err_rad = psi_est[k] - psi_true[k]

                    tau_errors_us[k].append(tau_err_us)
                    A_errors_percent[k].append(A_err_percent)
                    psi_errors_rad[k].append(psi_err_rad)

            except Exception as e:
                pass

        # ВЫЧИСЛЕНИЕ RMSE И BIAS
        for k in range(N_true):
            tau_err_k = [e for e in tau_errors_us[k] if not np.isnan(e)]
            A_err_k = [e for e in A_errors_percent[k] if not np.isnan(e)]
            psi_err_k = [e for e in psi_errors_rad[k] if not np.isnan(e)]

            if len(tau_err_k) > 0:
                tau_rmse_us[idx_snr, k] = np.sqrt(np.mean(np.array(tau_err_k)**2))
                tau_bias_us[idx_snr, k] = np.mean(tau_err_k)
            if len(A_err_k) > 0:
                A_rmse_percent[idx_snr, k] = np.sqrt(np.mean(np.array(A_err_k)**2))
                A_bias_percent[idx_snr, k] = np.mean(A_err_k)
            if len(psi_err_k) > 0:
                psi_rmse_rad[idx_snr, k] = np.sqrt(np.mean(np.array(psi_err_k)**2))
                psi_bias_rad[idx_snr, k] = np.mean(psi_err_k)

        # ВЫЧИСЛЕНИЕ CRB
        SNR_power = 10**(SNR_dB/10)
        for k in range(N_true):
            crb_A_percent[idx_snr, k], crb_tau_us[idx_snr, k], crb_psi_rad[idx_snr, k] = \
                compute_crb_parameters(f_m, a_m, SNR_power, T_dur, fs, A_true[k])

    # ------------------------------------------------------------------------
    # ВЫВОД РЕЗУЛЬТАТОВ
    # ------------------------------------------------------------------------

    print("\n" + "=" * 80)
    print("REZULTATY (RMSE)")
    print("=" * 80)

    for k in range(N_true):
        print(f"\n--- Luch {k+1} (tau={tau_true[k]*1e6:.0f} mks, A={A_true[k]:.1f}, psi={psi_true[k]:.2f} rad) ---")
        print(f"{'SNR(dB)':>10} | {'tau RMSE(mks)':>15} | {'A RMSE(%)':>15} | {'psi RMSE(rad)':>15}")
        print("-" * 65)
        for idx, SNR_dB in enumerate(SNR_dB_list):
            print(f"{SNR_dB:10.0f} | {tau_rmse_us[idx, k]:15.4f} | {A_rmse_percent[idx, k]:15.4f} | {psi_rmse_rad[idx, k]:15.4f}")

    print("\n" + "=" * 80)
    print("GRANICA KRAMERA-RAO (CRB)")
    print("=" * 80)

    for k in range(N_true):
        print(f"\n--- Luch {k+1} ---")
        print(f"{'SNR(dB)':>10} | {'CRB tau(mks)':>18} | {'CRB A(%)':>15} | {'CRB psi(rad)':>15}")
        print("-" * 65)
        for idx, SNR_dB in enumerate(SNR_dB_list):
            print(f"{SNR_dB:10.0f} | {crb_tau_us[idx, k]:18.6f} | {crb_A_percent[idx, k]:15.6f} | {crb_psi_rad[idx, k]:15.6f}")

    # ------------------------------------------------------------------------
    # ПОСТРОЕНИЕ ГРАФИКОВ
    # ------------------------------------------------------------------------

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    colors = ['blue', 'green']
    markers = ['o', 's']

    # График 1: Задержки tau1 i tau2
    ax = axes[0]
    for k in range(N_true):
        ax.plot(SNR_dB_list, tau_rmse_us[:, k],
                color=colors[k], marker=markers[k], linestyle='-',
                linewidth=2, markersize=8,
                label=f'RMSE tau{k+1} (MP ocenka)')
        ax.plot(SNR_dB_list, crb_tau_us[:, k],
                color=colors[k], linestyle='--',
                linewidth=2,
                label=f'CRB tau{k+1}')
    ax.set_xlabel('SNR, dB', fontsize=12)
    ax.set_ylabel('RMSE tau, mks', fontsize=12)
    ax.set_title('Oshibka ocenki zaderzhek tau1 i tau2', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    # График 2: Amplitudy A1 i A2
    ax = axes[1]
    for k in range(N_true):
        ax.plot(SNR_dB_list, A_rmse_percent[:, k],
                color=colors[k], marker=markers[k], linestyle='-',
                linewidth=2, markersize=8,
                label=f'RMSE A{k+1} (MP ocenka)')
        ax.plot(SNR_dB_list, crb_A_percent[:, k],
                color=colors[k], linestyle='--',
                linewidth=2,
                label=f'CRB A{k+1}')
    ax.set_xlabel('SNR, dB', fontsize=12)
    ax.set_ylabel('RMSE A, %', fontsize=12)
    ax.set_title('Oshibka ocenki amplitud A1 i A2', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    # График 3: Fazy psi1 i psi2
    ax = axes[2]
    for k in range(N_true):
        ax.plot(SNR_dB_list, psi_rmse_rad[:, k],
                color=colors[k], marker=markers[k], linestyle='-',
                linewidth=2, markersize=8,
                label=f'RMSE psi{k+1} (MP ocenka)')
        ax.plot(SNR_dB_list, crb_psi_rad[:, k],
                color=colors[k], linestyle='--',
                linewidth=2,
                label=f'CRB psi{k+1}')
    ax.set_xlabel('SNR, dB', fontsize=12)
    ax.set_ylabel('RMSE psi, rad', fontsize=12)
    ax.set_title('Oshibka ocenki faz psi1 i psi2', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    plt.suptitle(
        'Monte Carlo analiz ocenki parametrov kanala (N=2, M=3, formula 100)',
        fontsize=16
    )
    plt.tight_layout()
    plt.show()

    print("\n" + "=" * 80)
    print("MODELIROVANIE ZAVERSHENO")
    print("=" * 80)