"""
ТЕСТ: ОДИН ЛУЧ, ДВЕ ЧАСТОТЫ
Метод максимального правдоподобия (глава 4.2)
"""

import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response


def compute_matrices_C_S(omega_m, phi_m, tau):
    """Формула (92)"""
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
    """Формула (95)"""
    M, N = C.shape
    H = np.zeros((2 * M, 2 * N))
    H[:M, :N] = C
    H[:M, N:] = -S
    H[M:, :N] = S
    H[M:, N:] = C
    return H


def compute_matrix_Q(omega_m, a_m, tau, N0):
    """Формулы (80-82) - БЕЗ множителя T_dur/2"""
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


def compute_sufficient_statistics(xi, t, omega_m, a_m, N0, T_dur, fs):
    """Формулы (87-88) и (96)"""
    M = len(omega_m)
    dt = 1 / fs

    x_m = np.zeros(M)
    y_m = np.zeros(M)

    for m in range(M):
        cos_sum = np.sum(xi * np.cos(omega_m[m] * t)) * dt
        sin_sum = np.sum(xi * np.sin(omega_m[m] * t)) * dt
        x_m[m] = (2.0 / N0) * cos_sum
        y_m[m] = (2.0 / N0) * sin_sum

    # Z' = [a₁x₁, a₂x₂, ..., a_Mx_M, a₁y₁, a₂y₂, ..., a_My_M] (формула 96)
    Z_prime = np.concatenate([a_m * x_m, a_m * y_m])

    return Z_prime, x_m, y_m


def compute_likelihood(Z_prime, H, Q):
    """Формула (100): L(τ) = 0.5 * Z' * H * Q^{-1} * H^T * (Z')^T"""
    try:
        Q_inv = np.linalg.inv(Q)
        L = 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
        return float(L)
    except:
        return -np.inf


def estimate_amplitudes_and_phases(Z_prime, H, Q, T_dur):
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv

    # Правильная нормировка: умножаем на 2/T_dur
    A_cs = A_cs * (2.0 / T_dur)

    N = len(A_cs) // 2
    A_c = A_cs[:N]
    A_s = A_cs[N:]

    A_est = np.sqrt(A_c ** 2 + A_s ** 2)
    psi_est = np.arctan2(A_s, A_c)

    return A_est[0], psi_est[0]


if __name__ == "__main__":

    print("=" * 80)
    print("ТЕСТ: ОДИН ЛУЧ, ДВЕ ЧАСТОТЫ")
    print("Метод максимального правдоподобия (глава 4.2)")
    print("=" * 80)

    # Параметры
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
    SNR_dB = 10

    print(f"\nПАРАМЕТРЫ:")
    print(f"  f₁ = {f_m[0]} Гц, f₂ = {f_m[1]} Гц")
    print(f"  fs = {fs} Гц, T = {T_dur * 1000:.0f} мс")
    print(f"  SNR = {SNR_dB} дБ")
    print(f"  Истинные: τ = {tau_true * 1e6:.0f} мкс, ψ = {psi_true:.3f} рад, A = {A_true:.3f}")

    # Генерация данных
    print("\n1. ГЕНЕРАЦИЯ ДАННЫХ")
    data = generate_test_signal_and_channel_response(
        M, a_m, omega_m, phi_m,
        N_true, np.array([A_true]), np.array([tau_true]), np.array([psi_true]),
        fs, T1, T2, SNR_dB
    )

    # Достаточная статистика
    print("\n2. ДОСТАТОЧНАЯ СТАТИСТИКА (формулы 87-88, 96)")
    Z_prime, x_m, y_m = compute_sufficient_statistics(
        data['xi_observed'], data['t'], omega_m, a_m, N0, T_dur, fs
    )

    print(f"  x = {x_m}")
    print(f"  y = {y_m}")
    print(f"  Z' = {Z_prime}")

    # Поиск задержки
    print("\n3. ПОИСК ЗАДЕРЖКИ (формула 100)")
    tau_range = np.linspace(tau_true - 300e-6, tau_true + 300e-6, 201)
    L_vals = []

    for tau_test in tau_range:
        tau = np.array([tau_test])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L_vals.append(compute_likelihood(Z_prime, H, Q))

    L_vals = np.array(L_vals)
    tau_est = tau_range[np.argmax(L_vals)]
    L_max = np.max(L_vals)

    print(f"\n  РЕЗУЛЬТАТЫ ПОИСКА:")
    print(f"    Истинная τ = {tau_true * 1e6:.0f} мкс")
    print(f"    Оценка τ = {tau_est * 1e6:.1f} мкс")
    print(f"    Ошибка = {(tau_est - tau_true) * 1e6:.1f} мкс")
    print(f"    L_max = {L_max:.6f}")

    # Оценка амплитуды и фазы
    print("\n4. ОЦЕНКА АМПЛИТУДЫ И ФАЗЫ (формулы 99, 26)")
    C_opt, S_opt = compute_matrices_C_S(omega_m, phi_m, np.array([tau_est]))
    H_opt = compute_matrix_H(C_opt, S_opt)
    Q_opt = compute_matrix_Q(omega_m, a_m, np.array([tau_est]), N0)

    A_est, psi_est = estimate_amplitudes_and_phases(Z_prime, H_opt, Q_opt, T_dur)

    print(f"\n  РЕЗУЛЬТАТЫ:")
    print(
        f"    Амплитуда: истинная {A_true:.3f} → оценка {A_est:.3f} (ошибка {abs(A_est - A_true) / A_true * 100:.1f}%)")
    print(f"    Фаза ψ: истинная {psi_true:.3f} → оценка {psi_est:.3f} (ошибка {abs(psi_est - psi_true):.3f} рад)")

    # Визуализация
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # График L(τ)
    axes[0].plot(tau_range * 1e6, L_vals, 'b-', linewidth=1.5)
    axes[0].axvline(tau_true * 1e6, color='r', linestyle='--', label=f'Истинное τ = {tau_true * 1e6:.0f} мкс')
    axes[0].axvline(tau_est * 1e6, color='g', linestyle='--', label=f'Оценка τ = {tau_est * 1e6:.1f} мкс')
    axes[0].set_xlabel('τ, мкс')
    axes[0].set_ylabel('L(τ)')
    axes[0].set_title('Решающая статистика L(τ) (формула 100)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Комплексные амплитуды
    # Теоретические h
    h1_th = a_m[0] * A_true * np.exp(-1j * (omega_m[0] * tau_true + phi_m[0] + psi_true))
    h2_th = a_m[1] * A_true * np.exp(-1j * (omega_m[1] * tau_true + phi_m[1] + psi_true))

    # Оцененные h из достаточной статистики
    h1_est = x_m[0] - 1j * y_m[0]
    h2_est = x_m[1] - 1j * y_m[1]

    axes[1].plot(np.real(h1_th), np.imag(h1_th), 'r*', markersize=15, label='Теоретическая h₁')
    axes[1].plot(np.real(h2_th), np.imag(h2_th), 'b*', markersize=15, label='Теоретическая h₂')
    axes[1].plot(np.real(h1_est), np.imag(h1_est), 'ro', markersize=10, label='Оценка h₁')
    axes[1].plot(np.real(h2_est), np.imag(h2_est), 'bo', markersize=10, label='Оценка h₂')
    axes[1].axhline(0, color='k', linestyle='-', alpha=0.3)
    axes[1].axvline(0, color='k', linestyle='-', alpha=0.3)
    axes[1].set_xlabel('Re(h)')
    axes[1].set_ylabel('Im(h)')
    axes[1].set_title('Комплексные амплитуды h_m')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].axis('equal')

    plt.tight_layout()
    plt.show()

    print("\n" + "=" * 80)
    print("ТЕСТ ЗАВЕРШЕН")
    print("=" * 80)