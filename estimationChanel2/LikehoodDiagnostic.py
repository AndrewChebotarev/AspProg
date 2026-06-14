"""
Скрипт 4: ДИАГНОСТИКА РЕШАЮЩЕЙ СТАТИСТИКИ (формула 100)
Только L(τ) = 0.5 * Z' * H * Q^{-1} * H^T * (Z')^T
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
    H = np.zeros((2*M, 2*N))
    H[:M, :N] = C
    H[:M, N:] = -S
    H[M:, :N] = S
    H[M:, N:] = C
    return H


def compute_matrix_Q(omega_m, a_m, tau, T_dur, N0):
    """Формулы (80-82) - возвращает ТОЛЬКО матрицу Q"""
    M = len(omega_m)
    N = len(tau)
    coef = T_dur / 2

    Qc = np.zeros((N, N))
    Qcs = np.zeros((N, N))

    for i in range(N):
        for k in range(N):
            sum_cos = 0.0
            sum_sin = 0.0
            for m in range(M):
                weight = 2 * a_m[m]**2 / N0
                delta = tau[i] - tau[k]
                sum_cos += weight * np.cos(omega_m[m] * delta)
                sum_sin += weight * np.sin(omega_m[m] * (tau[k] - tau[i]))
            Qc[i, k] = coef * sum_cos
            Qcs[i, k] = coef * sum_sin

    Q = np.vstack([np.hstack([Qc, Qcs]), np.hstack([Qcs.T, Qc])])
    return Q


def compute_sufficient_statistics(xi, t, omega_m, a_m, N0, T_dur, fs):
    """Формулы (87-88) и (96)"""
    M = len(omega_m)
    dt = 1/fs

    x_m = np.zeros(M)
    y_m = np.zeros(M)

    for m in range(M):
        cos_sum = np.sum(xi * np.cos(omega_m[m] * t)) * dt
        sin_sum = np.sum(xi * np.sin(omega_m[m] * t)) * dt
        x_m[m] = (2.0 / N0) * cos_sum / T_dur
        y_m[m] = (2.0 / N0) * sin_sum / T_dur

    # ВАЖНО: минус перед y_m для согласования с H
    Z_prime = np.concatenate([x_m, -y_m])

    return Z_prime, x_m, y_m


def compute_likelihood(Z_prime, H, Q):
    """Формула (100): L(τ) = 0.5 * Z' * H * Q^{-1} * H^T * (Z')^T"""
    try:
        Q_inv = np.linalg.inv(Q)
        L = 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
        return float(L)
    except np.linalg.LinAlgError:
        return -np.inf


def grid_search_delays(Z_prime, omega_m, a_m, phi_m, tau1_range, tau2_range, T_dur, N0):
    """Поиск максимального L(τ) перебором по сетке для N=2"""
    n1, n2 = len(tau1_range), len(tau2_range)
    L_grid = np.zeros((n1, n2))

    print(f"  Поиск на сетке {n1} x {n2} = {n1*n2} точек...")

    for i, t1 in enumerate(tau1_range):
        for j, t2 in enumerate(tau2_range):
            tau_test = np.array([t1, t2])
            C, S = compute_matrices_C_S(omega_m, phi_m, tau_test)
            H = compute_matrix_H(C, S)
            Q = compute_matrix_Q(omega_m, a_m, tau_test, T_dur, N0)
            L_grid[i, j] = compute_likelihood(Z_prime, H, Q)

            # Прогресс
            if (i * n2 + j) % (n1*n2 // 10 + 1) == 0:
                print(f"    Прогресс: {100 * (i*n2 + j) / (n1*n2):.0f}%")

    return L_grid


# ============================================================================
# ОСНОВНОЙ СКРИПТ
# ============================================================================

if __name__ == "__main__":

    print("=" * 80)
    print("ДИАГНОСТИКА РЕШАЮЩЕЙ СТАТИСТИКИ L(τ) (ФОРМУЛА 100)")
    print("=" * 80)

    # ===== ПАРАМЕТРЫ (как в предыдущих скриптах) =====
    M = 3
    f_m = np.array([110, 210, 310])
    omega_m = 2 * np.pi * f_m
    a_m = np.array([1.0, 0.8, 0.6])
    phi_m = np.array([0, np.pi/4, np.pi/3])

    N_true = 2
    A_true = np.array([0.7, 0.5])
    tau_true = np.array([600e-6, 1200e-6])
    psi_true = np.array([0.52, -0.53])

    fs = 100 * np.max(f_m)  # 6000 Гц
    T1 = 0.0
    T2 = 0.05
    T_dur = T2 - T1
    N0 = 1.0
    SNR_dB = None  # Сначала без шума для проверки

    print(f"\nПАРАМЕТРЫ МОДЕЛИ:")
    print(f"  fs = {fs} Гц, T_dur = {T_dur*1000:.1f} мс")
    print(f"  M = {M}, N = {N_true}")
    print(f"  N0 = {N0}")
    print(f"  SNR = {SNR_dB if SNR_dB else 'без шума'}")
    print(f"  Истинные τ₁ = {tau_true[0]*1e6:.0f} мкс, τ₂ = {tau_true[1]*1e6:.0f} мкс")

    # ===== ГЕНЕРАЦИЯ ДАННЫХ =====
    print("\n" + "-" * 80)
    print("ШАГ 1: ГЕНЕРАЦИЯ ДАННЫХ")
    print("-" * 80)

    data = generate_test_signal_and_channel_response(
        M, a_m, omega_m, phi_m, N_true, A_true, tau_true, psi_true,
        fs, T1, T2, SNR_dB
    )

    # ===== ДОСТАТОЧНАЯ СТАТИСТИКА =====
    print("\n" + "-" * 80)
    print("ШАГ 2: ДОСТАТОЧНАЯ СТАТИСТИКА (формулы 87-88, 96)")
    print("-" * 80)

    Z_prime, x_m, y_m = compute_sufficient_statistics(
        data['xi_observed'], data['t'], omega_m, a_m, N0, T_dur, fs
    )

    print(f"  x_m = {x_m}")
    print(f"  y_m = {y_m}")
    print(f"  Z' = {Z_prime}")

    # ===== РЕШАЮЩАЯ СТАТИСТИКА =====
    print("\n" + "-" * 80)
    print("ШАГ 3: РЕШАЮЩАЯ СТАТИСТИКА L(τ) (формула 100)")
    print("-" * 80)

    # Сетка для поиска
    tau1_range = np.linspace(tau_true[0] - 100e-6, tau_true[0] + 100e-6, 41)
    tau2_range = np.linspace(tau_true[1] - 100e-6, tau_true[1] + 100e-6, 41)

    L_grid = grid_search_delays(
        Z_prime, omega_m, a_m, phi_m, tau1_range, tau2_range, T_dur, N0
    )

    # Находим максимум
    max_idx = np.unravel_index(np.argmax(L_grid), L_grid.shape)
    tau1_max = tau1_range[max_idx[0]]
    tau2_max = tau2_range[max_idx[1]]
    L_max = L_grid[max_idx]

    print(f"\n  РЕЗУЛЬТАТЫ ПОИСКА МАКСИМУМА L(τ):")
    print(f"    Истинные τ₁ = {tau_true[0]*1e6:.0f} мкс, τ₂ = {tau_true[1]*1e6:.0f} мкс")
    print(f"    Максимум L(τ) при τ₁ = {tau1_max*1e6:.1f} мкс, τ₂ = {tau2_max*1e6:.1f} мкс")
    print(f"    Ошибка по τ₁: {(tau1_max - tau_true[0])*1e6:.1f} мкс")
    print(f"    Ошибка по τ₂: {(tau2_max - tau_true[1])*1e6:.1f} мкс")
    print(f"    Значение L_max = {L_max:.6f}")

    # ===== ВИЗУАЛИЗАЦИЯ =====
    print("\n" + "-" * 80)
    print("ШАГ 4: ВИЗУАЛИЗАЦИЯ L(τ₁, τ₂)")
    print("-" * 80)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Контурный график
    im = axes[0].contourf(tau1_range*1e6, tau2_range*1e6, L_grid.T, levels=50, cmap='viridis')
    axes[0].set_xlabel('τ₁, мкс')
    axes[0].set_ylabel('τ₂, мкс')
    axes[0].set_title('Решающая статистика L(τ₁, τ₂)')
    plt.colorbar(im, ax=axes[0])
    axes[0].plot(tau_true[0]*1e6, tau_true[1]*1e6, 'r*', markersize=15, label='Истинное')
    axes[0].plot(tau1_max*1e6, tau2_max*1e6, 'wo', markersize=10, label='Максимум')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # 3D поверхность
    ax3d = fig.add_subplot(1, 2, 2, projection='3d')
    T1g, T2g = np.meshgrid(tau1_range*1e6, tau2_range*1e6)
    ax3d.plot_surface(T1g, T2g, L_grid, cmap='viridis', alpha=0.8)
    ax3d.set_xlabel('τ₁, мкс')
    ax3d.set_ylabel('τ₂, мкс')
    ax3d.set_zlabel('L(τ₁, τ₂)')
    ax3d.set_title('3D поверхность L(τ)')

    plt.tight_layout()
    plt.show()

    # ===== ПРОВЕРКА ПРИ ИСТИННЫХ τ =====
    print("\n" + "-" * 80)
    print("ШАГ 5: ПРОВЕРКА L(τ) ПРИ ИСТИННЫХ ЗАДЕРЖКАХ")
    print("-" * 80)

    C_true, S_true = compute_matrices_C_S(omega_m, phi_m, tau_true)
    H_true = compute_matrix_H(C_true, S_true)
    Q_true = compute_matrix_Q(omega_m, a_m, tau_true, T_dur, N0)
    L_true = compute_likelihood(Z_prime, H_true, Q_true)

    print(f"  L(τ_истинные) = {L_true:.6f}")
    print(f"  L(τ_максимум) = {L_max:.6f}")
    print(f"  Отношение L_max / L_true = {L_max/L_true:.4f}")

    if L_max >= L_true:
        print("  ✅ Максимум найден правильно (L_max ≥ L_true)")
    else:
        print("  ❌ Проблема: L_max < L_true, максимум не найден")

    print("\n" + "=" * 80)
    print("ДИАГНОСТИКА ЗАВЕРШЕНА")
    print("=" * 80)