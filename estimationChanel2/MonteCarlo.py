"""
Monte Carlo моделирование для оценки параметров канала
Сравнение с границей Крамера-Рао (CRB) - ЛИНЕЙНЫЙ МАСШТАБ

ОСНОВНЫЕ ФОРМУЛЫ ИЗ СТАТЬИ:
- (1) Тестовый сигнал: S(t) = sum(a_m * cos(omega_m * t - phi_m))
- (3) Отклик канала: s_m(t) = sum(A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k))
- (4) Полезный сигнал: s(t) = sum(s_m(t))
- (31)-(32) Квадратурные корреляторы
- (87)-(88) Достаточные статистики
- (95) Матрица H
- (100) Логарифм ФОП
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

    ФОРМУЛА (92):
    C_mk = cos(omega_m * tau_k + phi_m)
    S_mk = sin(omega_m * tau_k + phi_m)

    Эти матрицы используются для формирования матрицы H (95)

    ПАРАМЕТРЫ:
    ----------
    omega_m : array (M,) - угловые частоты синусоид [рад/с]
    phi_m : array (M,) - начальные фазы синусоид [рад]
    tau : array (N,) - задержки лучей [с]

    ВОЗВРАЩАЕТ:
    -----------
    C : array (M, N) - матрица косинусов
    S : array (M, N) - матрица синусов
    """
    M = len(omega_m)
    N = len(tau)
    C = np.zeros((M, N))
    S = np.zeros((M, N))

    for m in range(M):
        for k in range(N):
            # ФОРМУЛА (92): фаза = omega_m * tau_k + phi_m
            phase = omega_m[m] * tau[k] + phi_m[m]
            C[m, k] = np.cos(phase)  # C_mk
            S[m, k] = np.sin(phase)  # S_mk

    return C, S


def compute_matrix_H(C, S):
    """
    ВЫЧИСЛЕНИЕ МАТРИЦЫ H (формула 95, стр. 24)

    ФОРМУЛА (95):
    H = [ C   -S ]
        [ S    C ]

    Размерность: (2M x 2N)

    Эта матрица используется в выражении для логарифма ФОП (98) и (100)

    ПАРАМЕТРЫ:
    ----------
    C : array (M, N) - матрица косинусов из (92)
    S : array (M, N) - матрица синусов из (92)

    ВОЗВРАЩАЕТ:
    -----------
    H : array (2M, 2N) - блочная матрица
    """
    M, N = C.shape
    H = np.zeros((2 * M, 2 * N))

    # Верхняя строка блоков: [C, -S]
    H[:M, :N] = C          # левый верхний блок
    H[:M, N:] = -S         # правый верхний блок

    # Нижняя строка блоков: [S, C]
    H[M:, :N] = S          # левый нижний блок
    H[M:, N:] = C          # правый нижний блок

    return H


def compute_matrix_Q(omega_m, a_m, tau, N0):
    """
    ВЫЧИСЛЕНИЕ МАТРИЦЫ Q (формулы 80-82, стр. 21)

    После ортогонализации гармоник, элементы матриц Qc, Qs, Qcs:

    ФОРМУЛА (80): Qc_ik = (T2-T1)/2 * sum_m(2*a_m^2/N0 * cos(omega_m*(tau_i - tau_k)))
    ФОРМУЛА (81): Qs_ik = (T2-T1)/2 * sum_m(2*a_m^2/N0 * cos(omega_m*(tau_i - tau_k)))
    ФОРМУЛА (82): Qcs_ik = (T2-T1)/2 * sum_m(2*a_m^2/N0 * sin(omega_m*(tau_k - tau_i)))

    Блочная матрица Q (47):
    Q = [ Qc    Qcs ]
        [ Qsc   Qs  ]

    где Qsc = Qcs^T

    ПАРАМЕТРЫ:
    ----------
    omega_m : array (M,) - угловые частоты
    a_m : array (M,) - амплитуды синусоид
    tau : array (N,) - задержки лучей
    N0 : float - спектральная плотность шума

    ВОЗВРАЩАЕТ:
    -----------
    Q : array (2N, 2N) - блочная матрица Q
    """
    M = len(omega_m)
    N = len(tau)

    # Инициализация блоков Qc и Qcs (N x N)
    Qc = np.zeros((N, N))
    Qcs = np.zeros((N, N))

    for i in range(N):
        for k in range(N):
            sum_cos = 0.0
            sum_sin = 0.0

            for m in range(M):
                # Множитель из формул (80)-(82): 2*a_m^2/N0
                weight = 2 * a_m[m] ** 2 / N0

                # Разность задержек для аргумента косинуса
                delta = tau[i] - tau[k]

                # ФОРМУЛА (80): cos(omega_m * (tau_i - tau_k))
                sum_cos += weight * np.cos(omega_m[m] * delta)

                # ФОРМУЛА (82): sin(omega_m * (tau_k - tau_i))
                sum_sin += weight * np.sin(omega_m[m] * (tau[k] - tau[i]))

            Qc[i, k] = sum_cos      # Qc_ik
            Qcs[i, k] = sum_sin     # Qcs_ik

    # ФОРМУЛА (47): блочная матрица Q
    # Q = [ Qc    Qcs ]
    #     [ Qsc   Qs  ]
    # где Qsc = Qcs^T, Qs = Qc (из формул 80-81)
    Q = np.vstack([
        np.hstack([Qc, Qcs]),
        np.hstack([Qcs.T, Qc])
    ])

    return Q


def compute_sufficient_statistics(xi, t, omega_m, a_m, N0, fs):
    """
    ВЫЧИСЛЕНИЕ ДОСТАТОЧНЫХ СТАТИСТИК (формулы 87-88, стр. 23)

    ФОРМУЛА (87): x_m = (2/N0) * integral(T1,T2) xi(t) * cos(omega_m * t) dt
    ФОРМУЛА (88): y_m = (2/N0) * integral(T1,T2) xi(t) * sin(omega_m * t) dt

    Эти статистики не зависят от оцениваемых параметров канала,
    что позволяет исключить многоканальность по tau_k.

    Вектор Z' (96): Z' = [a_1*x_1, ..., a_M*x_M, a_1*y_1, ..., a_M*y_M]

    ПАРАМЕТРЫ:
    ----------
    xi : array (K,) - наблюдаемый сигнал (2)
    t : array (K,) - временная сетка
    omega_m : array (M,) - угловые частоты
    a_m : array (M,) - амплитуды синусоид
    N0 : float - спектральная плотность шума
    fs : float - частота дискретизации

    ВОЗВРАЩАЕТ:
    -----------
    Z_prime : array (2M,) - вектор достаточных статистик (96)
    x_m : array (M,) - статистики x_m (87)
    y_m : array (M,) - статистики y_m (88)
    """
    M = len(omega_m)
    dt = 1 / fs  # шаг дискретизации для численного интегрирования

    x_m = np.zeros(M)
    y_m = np.zeros(M)

    for m in range(M):
        # Численное интегрирование по формуле прямоугольников

        # ФОРМУЛА (87): интеграл xi(t) * cos(omega_m * t) dt
        cos_sum = np.sum(xi * np.cos(omega_m[m] * t)) * dt

        # ФОРМУЛА (88): интеграл xi(t) * sin(omega_m * t) dt
        sin_sum = np.sum(xi * np.sin(omega_m[m] * t)) * dt

        # Множитель 2/N0 из формул (87)-(88)
        x_m[m] = (2.0 / N0) * cos_sum
        y_m[m] = (2.0 / N0) * sin_sum

    # ФОРМУЛА (96): Z' = || X' Y' ||, где X'_m = a_m * x_m, Y'_m = a_m * y_m
    Z_prime = np.concatenate([a_m * x_m, a_m * y_m])

    return Z_prime, x_m, y_m


def compute_likelihood(Z_prime, H, Q):
    """
    ВЫЧИСЛЕНИЕ ЛОГАРИФМА ФОП (формула 100, стр. 25)

    ФОРМУЛА (100): L(T) = (1/2) * Z' * H * Q^{-1} * H^T * Z'^T

    Это решающая статистика для оценки временных положений лучей.

    ПАРАМЕТРЫ:
    ----------
    Z_prime : array (2M,) - вектор достаточных статистик (96)
    H : array (2M, 2N) - матрица H (95)
    Q : array (2N, 2N) - матрица Q (47)

    ВОЗВРАЩАЕТ:
    -----------
    L : float - значение логарифма ФОП
    """
    try:
        # Вычисление Q^{-1}
        Q_inv = np.linalg.inv(Q)

        # ФОРМУЛА (100): L = 0.5 * Z' * H * Q^{-1} * H^T * Z'^T
        return 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    except:
        # Если матрица вырождена
        return -np.inf


def estimate_amplitudes_and_phases(Z_prime, H, Q, T_dur):
    """
    ОЦЕНКА АМПЛИТУД И ФАЗ ЛУЧЕЙ (формула 99 и 102, стр. 25)

    ФОРМУЛА (99): A_cs = Z' * H * Q^{-1}
    ФОРМУЛА (102): A_cs_hat = Z' * H_hat * Q_hat^{-1}

    После нахождения A_c и A_s, амплитуда и фаза восстанавливаются:
    A_k = sqrt(A_ck^2 + A_sk^2)
    psi_k = arctan2(A_sk, A_ck)

    ПАРАМЕТРЫ:
    ----------
    Z_prime : array (2M,) - вектор достаточных статистик (96)
    H : array (2M, 2N) - матрица H (95)
    Q : array (2N, 2N) - матрица Q (47)
    T_dur : float - длительность наблюдения T2 - T1

    ВОЗВРАЩАЕТ:
    -----------
    A_est : array (N,) - оценки амплитуд лучей
    psi_est : array (N,) - оценки фаз лучей [рад]
    """
    Q_inv = np.linalg.inv(Q)

    # ФОРМУЛА (99): A_cs = Z' * H * Q^{-1}
    A_cs = Z_prime @ H @ Q_inv

    # Нормировка на длительность (из формул 80-82, где T_dur = T2 - T1)
    A_cs = A_cs * (2.0 / T_dur)

    N = len(A_cs) // 2
    A_c = A_cs[:N]  # A_ck = A_k * cos(psi_k) (26)
    A_s = A_cs[N:]  # A_sk = A_k * sin(psi_k) (26)

    # Восстановление амплитуды и фазы
    A_est = np.sqrt(A_c ** 2 + A_s ** 2)
    psi_est = np.arctan2(A_s, A_c)

    return A_est, psi_est


def estimate_delay(Z_prime, omega_m, a_m, phi_m, tau_true, T_dur, N0):
    """
    ОЦЕНКА ЗАДЕРЖКИ ЛУЧА (формула 100-101, стр. 25)

    ФОРМУЛА (101): tau_hat = arg max L(T)

    Где L(T) определяется формулой (100):
    L(T) = (1/2) * Z' * H * Q^{-1} * H^T * Z'^T

    Поиск максимума производится на интервале tau_true ± 100 мкс.
    Для повышения точности используется параболическая интерполяция.

    ПАРАМЕТРЫ:
    ----------
    Z_prime : array (2M,) - вектор достаточных статистик (96)
    omega_m : array (M,) - угловые частоты
    a_m : array (M,) - амплитуды синусоид
    phi_m : array (M,) - начальные фазы
    tau_true : float - истинное значение задержки [с]
    T_dur : float - длительность наблюдения
    N0 : float - спектральная плотность шума

    ВОЗВРАЩАЕТ:
    -----------
    tau_est : float - оценка задержки [с]
    """
    # Интервал поиска: tau_true ± 100 мкс
    tau_range = np.linspace(tau_true - 100e-6, tau_true + 100e-6, 201)
    L_vals = []

    # Вычисление L(T) для каждого tau (формула 100)
    for tau_test in tau_range:
        tau = np.array([tau_test])

        # ФОРМУЛА (92): вычисление C и S
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)

        # ФОРМУЛА (95): вычисление H
        H = compute_matrix_H(C, S)

        # ФОРМУЛА (80-82, 47): вычисление Q
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)

        # ФОРМУЛА (100): L(T) = 0.5 * Z' * H * Q^{-1} * H^T * Z'^T
        L_vals.append(compute_likelihood(Z_prime, H, Q))

    L_vals = np.array(L_vals)

    # Поиск максимума L(T) (формула 101)
    idx_max = np.argmax(L_vals)
    tau_est = tau_range[idx_max]

    # Параболическая интерполяция для повышения точности
    if 0 < idx_max < len(L_vals) - 1:
        x = tau_range[idx_max-1:idx_max+2]
        y = L_vals[idx_max-1:idx_max+2]

        # Коэффициенты параболы y = a*x^2 + b*x + c
        a = (y[0] - 2*y[1] + y[2]) / (2 * (x[0] - x[1])**2)
        b = (y[2] - y[0]) / (2 * (x[2] - x[0]))

        # Вершина параболы: x_vertex = -b/(2*a)
        if a != 0:
            tau_est = x[1] - b / (2 * a)

    return tau_est


def compute_crb_parameters(f_m, a_m, SNR_power, T_dur, fs, A):
    """
    ВЫЧИСЛЕНИЕ ГРАНИЦЫ КРАМЕРА-РАО (CRB)

    CRB для параметров многолучевого канала:

    Для амплитуды A:
    CRB_A = A / sqrt(N_samples * SNR * sum(a_m^2))

    Для задержки tau:
    CRB_tau = 1 / sqrt(N_samples * SNR * sum(a_m^2) * sigma_omega^2)

    Для фазы psi:
    CRB_psi = sqrt(mean_omega^2) / sqrt(N_samples * SNR * sum(a_m^2) * sigma_omega^2)

    где:
    sum_a2 = sum(a_m^2)
    mean_omega = sum(omega_m * a_m^2) / sum_a2
    mean_omega2 = sum(omega_m^2 * a_m^2) / sum_a2
    sigma_omega2 = mean_omega2 - mean_omega^2

    ПАРАМЕТРЫ:
    ----------
    f_m : array (M,) - частоты синусоид [Гц]
    a_m : array (M,) - амплитуды синусоид
    SNR_power : float - отношение сигнал/шум по мощности (линейное)
    T_dur : float - длительность наблюдения [с]
    fs : float - частота дискретизации [Гц]
    A : float - амплитуда луча

    ВОЗВРАЩАЕТ:
    -----------
    crb_A_percent : float - CRB для амплитуды [%]
    crb_tau_us : float - CRB для задержки [мкс]
    crb_psi_rad : float - CRB для фазы [рад]
    """
    N_samples = T_dur * fs
    M = len(f_m)
    omega_m = 2 * np.pi * f_m

    # Вычисление моментов частот
    sum_a2 = np.sum(a_m**2)
    mean_omega = np.sum(omega_m * a_m**2) / sum_a2
    mean_omega2 = np.sum(omega_m**2 * a_m**2) / sum_a2
    sigma_omega2 = mean_omega2 - mean_omega**2

    # CRB для амплитуды
    crb_A_abs = A / np.sqrt(N_samples * SNR_power * sum_a2)
    crb_A_percent = crb_A_abs / A * 100

    # CRB для задержки
    crb_tau_sec = 1 / np.sqrt(N_samples * SNR_power * sum_a2 * sigma_omega2)
    crb_tau_us = crb_tau_sec * 1e6

    # CRB для фазы
    crb_psi_rad = np.sqrt(mean_omega2) / np.sqrt(N_samples * SNR_power * sum_a2 * sigma_omega2)

    return crb_A_percent, crb_tau_us, crb_psi_rad


# ============================================================================
# ГЛАВНАЯ ФУНКЦИЯ - MONTE CARLO МОДЕЛИРОВАНИЕ
# ============================================================================

if __name__ == "__main__":

    print("=" * 80)
    print("MONTE CARLO МОДЕЛИРОВАНИЕ - ЛИНЕЙНЫЙ МАСШТАБ")
    print("=" * 80)

    # ------------------------------------------------------------------------
    # ПАРАМЕТРЫ МОДЕЛИ (из раздела 1 статьи)
    # ------------------------------------------------------------------------

    # Тестовый сигнал S(t) (формула 1, стр. 2)
    M = 2                                    # количество синусоид
    f_m = np.array([900, 1100])              # частоты [Гц]
    omega_m = 2 * np.pi * f_m                # угловые частоты [рад/с]
    a_m = np.array([1.0, 1.0])               # амплитуды синусоид
    phi_m = np.array([0, 0])                 # начальные фазы [рад]

    # Параметры канала (лучи)
    N_true = 1                               # количество лучей
    A_true = 0.7                             # амплитуда луча
    tau_true = 600e-6                        # задержка [с] (600 мкс)
    psi_true = 0.52                          # фаза луча [рад]

    # Параметры дискретизации и наблюдения
    fs = 20000                               # частота дискретизации [Гц]
    T1 = 0.0                                 # начало наблюдения [с]
    T2 = 1.0                                 # конец наблюдения [с]
    T_dur = T2 - T1                          # длительность наблюдения
    N0 = 100.0                               # спектральная плотность шума

    # Параметры моделирования
    SNR_dB_list = [-30, -20, -10, -5, 0, 5, 10, 15, 20, 30, 40]
    n_iterations = 200                       # число итераций Monte Carlo

    # ------------------------------------------------------------------------
    # ИНИЦИАЛИЗАЦИЯ МАССИВОВ ДЛЯ РЕЗУЛЬТАТОВ
    # ------------------------------------------------------------------------

    # RMSE (Root Mean Square Error)
    tau_rmse_us = np.zeros(len(SNR_dB_list))
    A_rmse_percent = np.zeros(len(SNR_dB_list))
    psi_rmse_rad = np.zeros(len(SNR_dB_list))

    # Bias (систематическая ошибка)
    tau_bias_us = np.zeros(len(SNR_dB_list))
    A_bias_percent = np.zeros(len(SNR_dB_list))
    psi_bias_rad = np.zeros(len(SNR_dB_list))

    # CRB (граница Крамера-Рао)
    crb_tau_us = np.zeros(len(SNR_dB_list))
    crb_A_percent = np.zeros(len(SNR_dB_list))
    crb_psi_rad = np.zeros(len(SNR_dB_list))

    print(f"\nЗапуск Monte Carlo симуляции...")
    print(f"  Частоты: {f_m[0]} и {f_m[1]} Гц")
    print(f"  T_dur = {T_dur} с, fs = {fs} Гц")
    print(f"  N_samples = {T_dur * fs}")
    print(f"  SNR уровни: {SNR_dB_list}")
    print(f"  Итераций на SNR: {n_iterations}\n")

    # ------------------------------------------------------------------------
    # ОСНОВНОЙ ЦИКЛ ПО SNR
    # ------------------------------------------------------------------------

    for idx_snr, SNR_dB in enumerate(tqdm(SNR_dB_list, desc="Обработка SNR")):
        tau_errors_us = []
        A_errors_percent = []
        psi_errors_rad = []

        # --------------------------------------------------------------------
        # ЦИКЛ MONTE CARLO
        # --------------------------------------------------------------------
        for iter_num in range(n_iterations):
            try:
                # --------------------------------------------------------------
                # ШАГ 1: ГЕНЕРАЦИЯ ДАННЫХ (формулы 1-4, стр. 2)
                # --------------------------------------------------------------

                # ФОРМУЛА (1): S(t) = sum(a_m * cos(omega_m * t - phi_m))
                # ФОРМУЛА (3): s_m(t) = sum(A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k))
                # ФОРМУЛА (4): s(t) = sum(s_m(t))
                # ФОРМУЛА (2): xi(t) = s(t) + n(t)
                data = generate_test_signal_and_channel_response(
                    M, a_m, omega_m, phi_m,
                    N_true, np.array([A_true]), np.array([tau_true]), np.array([psi_true]),
                    fs, T1, T2, SNR_dB
                )

                # --------------------------------------------------------------
                # ШАГ 2: ВЫЧИСЛЕНИЕ ДОСТАТОЧНЫХ СТАТИСТИК (формулы 87-88, 96)
                # --------------------------------------------------------------

                # ФОРМУЛА (87): x_m = (2/N0) * integral(xi(t) * cos(omega_m * t)) dt
                # ФОРМУЛА (88): y_m = (2/N0) * integral(xi(t) * sin(omega_m * t)) dt
                # ФОРМУЛА (96): Z' = [a_1*x_1, ..., a_M*x_M, a_1*y_1, ..., a_M*y_M]
                Z_prime, x_m, y_m = compute_sufficient_statistics(
                    data['xi_observed'], data['t'], omega_m, a_m, N0, fs
                )

                # --------------------------------------------------------------
                # ШАГ 3: ОЦЕНКА ЗАДЕРЖКИ (формулы 100-101, стр. 25)
                # --------------------------------------------------------------

                # ФОРМУЛА (101): tau_hat = arg max L(T)
                # где L(T) из (100): L(T) = 0.5 * Z' * H * Q^{-1} * H^T * Z'^T
                tau_est_sec = estimate_delay(Z_prime, omega_m, a_m, phi_m, tau_true, T_dur, N0)

                # --------------------------------------------------------------
                # ШАГ 4: ОЦЕНКА АМПЛИТУДЫ И ФАЗЫ (формулы 99, 102)
                # --------------------------------------------------------------

                # ФОРМУЛА (92): C_mk = cos(omega_m * tau_k + phi_m), S_mk = sin(omega_m * tau_k + phi_m)
                tau_arr = np.array([tau_est_sec])
                C, S = compute_matrices_C_S(omega_m, phi_m, tau_arr)

                # ФОРМУЛА (95): H = [C, -S; S, C]
                H_opt = compute_matrix_H(C, S)

                # ФОРМУЛА (80-82, 47): Q = [Qc, Qcs; Qsc, Qs]
                Q_opt = compute_matrix_Q(omega_m, a_m, tau_arr, N0)

                # ФОРМУЛА (102): A_cs_hat = Z' * H_hat * Q_hat^{-1}
                # Затем: A_k = sqrt(A_ck^2 + A_sk^2), psi_k = arctan2(A_sk, A_ck)
                A_est, psi_est = estimate_amplitudes_and_phases(Z_prime, H_opt, Q_opt, T_dur)

                # --------------------------------------------------------------
                # ШАГ 5: ВЫЧИСЛЕНИЕ ОШИБОК
                # --------------------------------------------------------------

                tau_err_us = (tau_est_sec - tau_true) * 1e6
                A_err_percent = (A_est[0] - A_true) / A_true * 100
                psi_err_rad = psi_est[0] - psi_true

                tau_errors_us.append(tau_err_us)
                A_errors_percent.append(A_err_percent)
                psi_errors_rad.append(psi_err_rad)

            except Exception as e:
                # Пропускаем неудачные итерации
                pass

        # --------------------------------------------------------------------
        # ВЫЧИСЛЕНИЕ RMSE И BIAS
        # --------------------------------------------------------------------

        # Очистка от NaN значений
        tau_errors_us = [e for e in tau_errors_us if not np.isnan(e)]
        A_errors_percent = [e for e in A_errors_percent if not np.isnan(e)]
        psi_errors_rad = [e for e in psi_errors_rad if not np.isnan(e)]

        # RMSE (среднеквадратичная ошибка)
        if len(tau_errors_us) > 0:
            tau_rmse_us[idx_snr] = np.sqrt(np.mean(np.array(tau_errors_us)**2))
            tau_bias_us[idx_snr] = np.mean(tau_errors_us)
        if len(A_errors_percent) > 0:
            A_rmse_percent[idx_snr] = np.sqrt(np.mean(np.array(A_errors_percent)**2))
            A_bias_percent[idx_snr] = np.mean(A_errors_percent)
        if len(psi_errors_rad) > 0:
            psi_rmse_rad[idx_snr] = np.sqrt(np.mean(np.array(psi_errors_rad)**2))
            psi_bias_rad[idx_snr] = np.mean(psi_errors_rad)

        # --------------------------------------------------------------------
        # ВЫЧИСЛЕНИЕ CRB
        # --------------------------------------------------------------------

        SNR_power = 10**(SNR_dB/10)
        crb_A_percent[idx_snr], crb_tau_us[idx_snr], crb_psi_rad[idx_snr] = \
            compute_crb_parameters(f_m, a_m, SNR_power, T_dur, fs, A_true)

    # ------------------------------------------------------------------------
    # ВЫВОД РЕЗУЛЬТАТОВ
    # ------------------------------------------------------------------------

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

    # ------------------------------------------------------------------------
    # ПОСТРОЕНИЕ ГРАФИКОВ
    # ------------------------------------------------------------------------

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # График 1: Задержка τ
    ax = axes[0]
    ax.plot(SNR_dB_list, tau_rmse_us, 'bo-', linewidth=2, markersize=8, label='RMSE (МП оценка)')
    ax.plot(SNR_dB_list, crb_tau_us, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ', fontsize=12)
    ax.set_ylabel('RMSE τ, мкс', fontsize=12)
    ax.set_title('Ошибка оценки задержки τ', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    # График 2: Амплитуда A
    ax = axes[1]
    ax.plot(SNR_dB_list, A_rmse_percent, 'bo-', linewidth=2, markersize=8, label='RMSE (МП оценка)')
    ax.plot(SNR_dB_list, crb_A_percent, 'r--', linewidth=2, label='CRB')
    ax.set_xlabel('SNR, дБ', fontsize=12)
    ax.set_ylabel('RMSE A, %', fontsize=12)
    ax.set_title('Ошибка оценки амплитуды A', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    # График 3: Фаза ψ
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