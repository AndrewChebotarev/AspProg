"""
Диагностика матриц C, S, H, Q по статье 4.2
БЕЗ оценок амплитуд - только матрицы и их свойства
"""

import numpy as np
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
    """
    Формулы (80-82)
    Qc_ik = (T/2) * Σ (2a_m²/N0) * cos(ω_m (τ_i - τ_k))
    Qcs_ik = (T/2) * Σ (2a_m²/N0) * sin(ω_m (τ_k - τ_i))
    """
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
    return Q, Qc, Qcs


def print_matrix(name, matrix, precision=4):
    """Красивый вывод матрицы"""
    print(f"\n{name} (размер {matrix.shape}):")
    for i in range(matrix.shape[0]):
        row_str = "  ".join([f"{x:>{precision+4}.{precision}f}" for x in matrix[i, :]])
        print(f"  [{row_str}]")


# ============================================================================
# ОСНОВНОЙ СКРИПТ
# ============================================================================

if __name__ == "__main__":

    print("=" * 80)
    print("ДИАГНОСТИКА МАТРИЦ C, S, H, Q (ПО ФОРМУЛАМ СТАТЬИ)")
    print("=" * 80)

    # ===== ПАРАМЕТРЫ (как в первых двух скриптах) =====
    M = 3
    f_m = np.array([100, 200, 300])
    omega_m = 2 * np.pi * f_m
    a_m = np.array([1.0, 0.8, 0.6])
    phi_m = np.array([0, np.pi/4, np.pi/3])

    N_true = 2
    A_true = np.array([0.7, 0.5])
    tau_true = np.array([500e-6, 1200e-6])
    psi_true = np.array([0.52, -0.53])

    fs = 20 * np.max(f_m)  # 6000 Гц
    T1 = 0.0
    T2 = 0.01
    T_dur = T2 - T1
    N0 = 1.0

    print(f"\nПАРАМЕТРЫ МОДЕЛИ:")
    print(f"  fs = {fs} Гц")
    print(f"  T_dur = {T_dur*1000:.1f} мс")
    print(f"  M = {M}, N = {N_true}")
    print(f"  N0 = {N0}")
    print(f"  Истинные τ₁ = {tau_true[0]*1e6:.0f} мкс, τ₂ = {tau_true[1]*1e6:.0f} мкс")

    # ===== ВЫЧИСЛЕНИЕ МАТРИЦ =====
    print("\n" + "=" * 80)
    print("ВЫЧИСЛЕНИЕ МАТРИЦ ПО ФОРМУЛАМ СТАТЬИ")
    print("=" * 80)

    # Шаг 1: Матрицы C и S (формула 92)
    print("\n" + "-" * 80)
    print("ШАГ 1: МАТРИЦЫ C И S (формула 92)")
    print("-" * 80)

    C, S = compute_matrices_C_S(omega_m, phi_m, tau_true)

    print_matrix("C (cos(ω_m τ_k + φ_m))", C)
    print_matrix("S (sin(ω_m τ_k + φ_m))", S)

    # Проверка: C² + S² = 1
    check = C**2 + S**2
    print(f"\nПроверка C² + S² = 1:")
    for m in range(M):
        print(f"  m={m+1}: {check[m, :]}")

    # Шаг 2: Блочная матрица H (формула 95)
    print("\n" + "-" * 80)
    print("ШАГ 2: БЛОЧНАЯ МАТРИЦА H (формула 95)")
    print("-" * 80)

    H = compute_matrix_H(C, S)
    print_matrix("H = [[C, -S], [S, C]]", H)

    # Проверка структуры H
    M, N = C.shape  # Теперь N определено
    print(f"\nПроверка структуры H:")
    print(f"  Верхний левый блок (C): {H[:M, :N].shape}")
    print(f"  Верхний правый блок (-S): {H[:M, N:].shape}")
    print(f"  Нижний левый блок (S): {H[M:, :N].shape}")
    print(f"  Нижний правый блок (C): {H[M:, N:].shape}")

    # Шаг 3: Матрица Q (формулы 80-82)
    print("\n" + "-" * 80)
    print("ШАГ 3: МАТРИЦА Q (формулы 80-82)")
    print("-" * 80)

    Q, Qc, Qcs = compute_matrix_Q(omega_m, a_m, tau_true, T_dur, N0)

    print_matrix("Qc (блок косинусов)", Qc)
    print_matrix("Qcs (блок синусов)", Qcs)
    print_matrix("Q = [[Qc, Qcs], [Qcs^T, Qc]]", Q)

    # Свойства Q
    print(f"\nСВОЙСТВА МАТРИЦЫ Q:")
    print(f"  Симметричность Q: {np.allclose(Q, Q.T)}")
    print(f"  Определитель Q: {np.linalg.det(Q):.2e}")
    print(f"  Число обусловленности Q: {np.linalg.cond(Q):.2e}")
    print(f"  Ранг Q: {np.linalg.matrix_rank(Q)}")

    # Проверка: Qc должна быть симметричной
    print(f"\n  Qc симметрична: {np.allclose(Qc, Qc.T)}")
    print(f"  Qcs антисимметрична: {np.allclose(Qcs, -Qcs.T)}")

    # Дополнительная информация
    print("\n" + "-" * 80)
    print("ДОПОЛНИТЕЛЬНАЯ ИНФОРМАЦИЯ")
    print("-" * 80)

    print(f"\nИстинные значения для справки:")
    print(f"  A_true = {A_true}")
    print(f"  psi_true = {psi_true}")

    print("\nФормулы, использованные для вычислений:")
    print("  (92): C_mk = cos(ω_m τ_k + φ_m), S_mk = sin(ω_m τ_k + φ_m)")
    print("  (95): H = [[C, -S], [S, C]]")
    print("  (80-82): Qc_ik = (T/2)·Σ (2a_m²/N0)·cos(ω_m(τ_i-τ_k))")
    print("           Qcs_ik = (T/2)·Σ (2a_m²/N0)·sin(ω_m(τ_k-τ_i))")

    # Значения элементов Q для справки
    print("\nЗНАЧЕНИЯ ЭЛЕМЕНТОВ Q при истинных τ:")
    print(f"  Qc[0,0] = {Qc[0,0]:.6f}")
    print(f"  Qc[0,1] = {Qc[0,1]:.6f}")
    print(f"  Qc[1,0] = {Qc[1,0]:.6f}")
    print(f"  Qc[1,1] = {Qc[1,1]:.6f}")
    print(f"  Qcs[0,0] = {Qcs[0,0]:.6f}")
    print(f"  Qcs[0,1] = {Qcs[0,1]:.6f}")
    print(f"  Qcs[1,0] = {Qcs[1,0]:.6f}")
    print(f"  Qcs[1,1] = {Qcs[1,1]:.6f}")

    print("\n" + "=" * 80)
    print("ДИАГНОСТИКА ЗАВЕРШЕНА")
    print("=" * 80)