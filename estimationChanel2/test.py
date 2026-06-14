"""
ТЕСТ: ОДИН ЛУЧ, ДВЕ ЧАСТОТЫ
Оценка задержки τ и фазы ψ
"""

import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response


def compute_sufficient_statistics(xi, t, omega_m, a_m, N0, T_dur, fs):
    M = len(omega_m)
    dt = 1/fs

    x_m = np.zeros(M)
    y_m = np.zeros(M)

    for m in range(M):
        cos_sum = np.sum(xi * np.cos(omega_m[m] * t)) * dt
        sin_sum = np.sum(xi * np.sin(omega_m[m] * t)) * dt
        x_m[m] = (2.0 / N0) * cos_sum / T_dur
        y_m[m] = (2.0 / N0) * sin_sum / T_dur

    h_est = x_m - 1j * y_m

    return h_est, x_m, y_m


def estimate_delay_from_ratio(h1, h2, omega1, omega2):
    """
    Оценка задержки из отношения комплексных амплитуд
    h2/h1 = exp(-j(ω2-ω1)τ)
    """
    ratio = h2 / h1
    phase_diff = -np.angle(ratio)  # потому что exp(-jΔωτ)
    tau_est = phase_diff / (omega2 - omega1)
    return tau_est


if __name__ == "__main__":

    print("=" * 80)
    print("ТЕСТ: ОДИН ЛУЧ, ДВЕ ЧАСТОТЫ")
    print("Оценка задержки τ и фазы ψ")
    print("=" * 80)

    # Параметры
    M = 2
    f_m = np.array([900, 1100])  # две частоты
    omega_m = 2 * np.pi * f_m
    a_m = np.array([1.0, 1.0])   # равные амплитуды
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
    SNR_dB = 30

    print(f"\nПАРАМЕТРЫ:")
    print(f"  f₁ = {f_m[0]} Гц, f₂ = {f_m[1]} Гц")
    print(f"  fs = {fs} Гц, T = {T_dur*1000:.0f} мс")
    print(f"  SNR = {SNR_dB} дБ")
    print(f"  Истинные: τ = {tau_true*1e6:.0f} мкс, ψ = {psi_true:.3f} рад, A = {A_true:.3f}")

    # Теоретические h
    h1_theoretical = a_m[0] * A_true * np.exp(-1j * (omega_m[0] * tau_true + phi_m[0] + psi_true))
    h2_theoretical = a_m[1] * A_true * np.exp(-1j * (omega_m[1] * tau_true + phi_m[1] + psi_true))

    print(f"\nТеоретические:")
    print(f"  h₁ = {h1_theoretical:.6f}")
    print(f"  h₂ = {h2_theoretical:.6f}")

    # Генерация данных
    print("\n1. ГЕНЕРАЦИЯ ДАННЫХ")
    data = generate_test_signal_and_channel_response(
        M, a_m, omega_m, phi_m,
        N_true, np.array([A_true]), np.array([tau_true]), np.array([psi_true]),
        fs, T1, T2, SNR_dB
    )

    # Оценка h
    print("\n2. ОЦЕНКА КОМПЛЕКСНЫХ АМПЛИТУД")
    h_est, x_m, y_m = compute_sufficient_statistics(
        data['xi_observed'], data['t'], omega_m, a_m, N0, T_dur, fs
    )

    print(f"  h₁ = {h_est[0]:.6f}")
    print(f"  h₂ = {h_est[1]:.6f}")

    # Оценка задержки из отношения
    tau_est = estimate_delay_from_ratio(h_est[0], h_est[1], omega_m[0], omega_m[1])

    print(f"\n3. ОЦЕНКА ЗАДЕРЖКИ:")
    print(f"  Истинная τ = {tau_true*1e6:.0f} мкс")
    print(f"  Оценка τ = {tau_est*1e6:.1f} мкс")
    print(f"  Ошибка = {(tau_est - tau_true)*1e6:.1f} мкс")

    # Оценка амплитуды (усреднение)
    A_est = (np.abs(h_est[0]) + np.abs(h_est[1])) / 2

    # Оценка фазы ψ из первого канала
    # Коррекция фазы
    psi_est = -np.angle(h_est[0]) - omega_m[0] * tau_est
    # Приводим к диапазону [-π, π]
    psi_est = np.arctan2(np.sin(psi_est), np.cos(psi_est))

    print(f"\n4. ОЦЕНКА АМПЛИТУДЫ И ФАЗЫ:")
    print(f"  Амплитуда: истинная {A_true:.3f} → оценка {A_est:.3f} (ошибка {abs(A_est-A_true)/A_true*100:.1f}%)")
    print(f"  Фаза ψ: истинная {psi_true:.3f} → оценка {psi_est:.3f} (ошибка {abs(psi_est-psi_true):.3f} рад)")

    # Визуализация
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Комплексные амплитуды
    axes[0].plot(np.real(h1_theoretical), np.imag(h1_theoretical), 'r*', markersize=15, label='Теоретическая h₁')
    axes[0].plot(np.real(h2_theoretical), np.imag(h2_theoretical), 'b*', markersize=15, label='Теоретическая h₂')
    axes[0].plot(np.real(h_est[0]), np.imag(h_est[0]), 'ro', markersize=10, label='Оценка h₁')
    axes[0].plot(np.real(h_est[1]), np.imag(h_est[1]), 'bo', markersize=10, label='Оценка h₂')
    axes[0].axhline(0, color='k', linestyle='-', alpha=0.3)
    axes[0].axvline(0, color='k', linestyle='-', alpha=0.3)
    axes[0].set_xlabel('Re(h)')
    axes[0].set_ylabel('Im(h)')
    axes[0].set_title('Комплексные амплитуды')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].axis('equal')

    # Фазовая разность
    phase_diff_true = -np.angle(h2_theoretical / h1_theoretical)
    phase_diff_est = -np.angle(h_est[1] / h_est[0])

    axes[1].bar([0, 1], [phase_diff_true, phase_diff_est], tick_label=['Теория', 'Оценка'])
    axes[1].set_ylabel('Δφ = (ω₂-ω₁)τ, рад')
    axes[1].set_title('Разность фаз между частотами')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    print("\n" + "=" * 80)
    print("ТЕСТ ЗАВЕРШЕН")
    print("=" * 80)