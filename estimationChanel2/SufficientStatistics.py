"""
Этап 1: Вычисление достаточной статистики (формулы 87, 88)
С ПРАВИЛЬНЫМ масштабированием и знаком
"""

import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response


def compute_sufficient_statistics_fixed(xi_observed, t, omega_m, N0, fs):
    """
    Вычисление достаточной статистики по формулам (87)-(88)

    x_m = (2/N0) * (1/T) * ∫ ξ(t) cos(ω_m t) dt
    y_m = (2/N0) * (1/T) * ∫ ξ(t) sin(ω_m t) dt

    В комплексной форме: h_m = x_m - j*y_m (ВАЖНО: минус!)
    """

    M = len(omega_m)
    x_m = np.zeros(M)
    y_m = np.zeros(M)

    dt = 1/fs
    T_dur = t[-1] - t[0]

    print("\n" + "=" * 80)
    print("ВЫЧИСЛЕНИЕ ДОСТАТОЧНОЙ СТАТИСТИКИ")
    print("=" * 80)
    print(f"  N0 = {N0}")
    print(f"  dt = {dt*1e6:.2f} мкс")
    print(f"  T_dur = {T_dur*1000:.2f} мс")
    print(f"  Множитель 2/(N0*T_dur) = {2/(N0*T_dur):.4f}")
    print(f"  ВНИМАНИЕ: h_m = x_m - j*y_m")

    for m in range(M):
        cos_sum = np.sum(xi_observed * np.cos(omega_m[m] * t)) * dt
        sin_sum = np.sum(xi_observed * np.sin(omega_m[m] * t)) * dt

        x_m[m] = (2.0 / N0) * cos_sum / T_dur
        y_m[m] = (2.0 / N0) * sin_sum / T_dur

        f_hz = omega_m[m] / (2 * np.pi)
        print(f"\n  Гармоника {m+1} (f={f_hz:.0f} Гц):")
        print(f"    x_{m+1} = {x_m[m]:.6f}")
        print(f"    y_{m+1} = {y_m[m]:.6f}")

    return x_m, y_m


def compute_theoretical_h(a_m, omega_m, phi_m, A_k, tau_k, psi_k):
    """Теоретическое значение h_m по формуле (13)"""
    M = len(a_m)
    N = len(A_k)
    h_theoretical = np.zeros(M, dtype=complex)

    print("\n" + "=" * 80)
    print("ТЕОРЕТИЧЕСКИЕ ЗНАЧЕНИЯ h_m (формула 13)")
    print("=" * 80)

    for m in range(M):
        for k in range(N):
            phase = omega_m[m] * tau_k[k] + phi_m[m] + psi_k[k]
            h_theoretical[m] += a_m[m] * A_k[k] * np.exp(-1j * phase)

        f_hz = omega_m[m]/(2*np.pi)
        print(f"\n  Гармоника {m+1} (f={f_hz:.0f} Гц):")
        print(f"    h_m = {h_theoretical[m]:.6f}")
        print(f"    |h_m| = {np.abs(h_theoretical[m]):.6f}")
        print(f"    arg(h_m) = {np.angle(h_theoretical[m]):.6f} рад")

    return h_theoretical


# ========== ОСНОВНОЙ СКРИПТ ==========
if __name__ == "__main__":

    # Параметры тестового сигнала
    M = 3
    f_m = np.array([100, 200, 300])
    omega_m = 2 * np.pi * f_m
    a_m = np.array([1.0, 0.8, 0.6])
    phi_m = np.array([0, np.pi/4, np.pi/3])

    # Параметры канала
    N = 2
    A_k = np.array([0.7, 0.5])
    tau_k = np.array([500e-6, 1200e-6])
    psi_k = np.array([0.52, -0.53])

    # Параметры дискретизации
    max_freq = np.max(f_m)
    fs = 20 * max_freq
    T1 = 0.0
    T2 = 0.01

    N0 = 1.0  # Нормировка

    # ТЕСТ 1: БЕЗ ШУМА
    print("\n" + "=" * 80)
    print("ТЕСТ 1: БЕЗ ШУМА (проверка корректности)")
    print("=" * 80)

    data_no_noise = generate_test_signal_and_channel_response(
        M=M, a_m=a_m, omega_m=omega_m, phi_m=phi_m,
        N=N, A_k=A_k, tau_k=tau_k, psi_k=psi_k,
        fs=fs, T1=T1, T2=T2,
        SNR_dB=None
    )

    x_m, y_m = compute_sufficient_statistics_fixed(
        xi_observed=data_no_noise['xi_observed'],
        t=data_no_noise['t'],
        omega_m=omega_m,
        N0=N0,
        fs=fs
    )

    h_theoretical = compute_theoretical_h(a_m, omega_m, phi_m, A_k, tau_k, psi_k)

    # ВАЖНО: h_estimated = x_m - j*y_m (минус!)
    h_estimated = x_m - 1j * y_m

    print("\n" + "=" * 80)
    print("СРАВНЕНИЕ (без шума, N0=1)")
    print("=" * 80)

    errors = []
    for m in range(M):
        error = np.abs(h_theoretical[m] - h_estimated[m])
        rel_error = error / np.abs(h_theoretical[m]) * 100
        errors.append(error)
        print(f"\n  Гармоника {m+1}:")
        print(f"    Теория: {h_theoretical[m]:.6f}")
        print(f"    Оценка: {h_estimated[m]:.6f}")
        print(f"    Ошибка: {error:.8f} ({rel_error:.4f}%)")

    print(f"\n  Средняя ошибка: {np.mean(errors):.8f}")

    # ТЕСТ 2: С ШУМОМ
    print("\n" + "=" * 80)
    print("ТЕСТ 2: С ШУМОМ (SNR = 15 дБ)")
    print("=" * 80)

    SNR_dB = 15

    data_with_noise = generate_test_signal_and_channel_response(
        M=M, a_m=a_m, omega_m=omega_m, phi_m=phi_m,
        N=N, A_k=A_k, tau_k=tau_k, psi_k=psi_k,
        fs=fs, T1=T1, T2=T2,
        SNR_dB=SNR_dB
    )

    x_m_noise, y_m_noise = compute_sufficient_statistics_fixed(
        xi_observed=data_with_noise['xi_observed'],
        t=data_with_noise['t'],
        omega_m=omega_m,
        N0=N0,
        fs=fs
    )

    h_estimated_noise = x_m_noise - 1j * y_m_noise

    print("\n" + "=" * 80)
    print(f"СРАВНЕНИЕ (с шумом, SNR={SNR_dB} дБ, N0=1)")
    print("=" * 80)

    errors_noise = []
    for m in range(M):
        error = np.abs(h_theoretical[m] - h_estimated_noise[m])
        rel_error = error / np.abs(h_theoretical[m]) * 100
        errors_noise.append(error)
        print(f"\n  Гармоника {m+1}:")
        print(f"    Теория (без шума): {h_theoretical[m]:.6f}")
        print(f"    Оценка (с шумом):  {h_estimated_noise[m]:.6f}")
        print(f"    Ошибка: {error:.6f} ({rel_error:.2f}%)")

    print(f"\n  Средняя ошибка: {np.mean(errors_noise):.6f}")

    # Визуализация
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    indices = np.arange(1, M+1)
    width = 0.35

    # 1. x_m и y_m (с шумом)
    axes[0, 0].bar(indices - width/2, x_m_noise, width, label='x_m', alpha=0.7, color='blue')
    axes[0, 0].bar(indices + width/2, y_m_noise, width, label='y_m', alpha=0.7, color='red')
    axes[0, 0].set_xlabel('Номер гармоники m')
    axes[0, 0].set_ylabel('Значение')
    axes[0, 0].set_title(f'Достаточная статистика (SNR={SNR_dB} дБ)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_xticks(indices)

    # 2. Сравнение |h_m|
    axes[0, 1].bar(indices - width/2, np.abs(h_estimated_noise), width, label='Оценка (с шумом)', alpha=0.7, color='green')
    axes[0, 1].bar(indices + width/2, np.abs(h_theoretical), width, label='Теория (без шума)', alpha=0.7, color='orange')
    axes[0, 1].set_xlabel('Номер гармоники m')
    axes[0, 1].set_ylabel('|h_m|')
    axes[0, 1].set_title('Амплитуда комплексной огибающей')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_xticks(indices)

    # 3. Сравнение фазы
    axes[1, 0].bar(indices - width/2, np.angle(h_estimated_noise), width, label='Оценка (с шумом)', alpha=0.7, color='purple')
    axes[1, 0].bar(indices + width/2, np.angle(h_theoretical), width, label='Теория (без шума)', alpha=0.7, color='brown')
    axes[1, 0].set_xlabel('Номер гармоники m')
    axes[1, 0].set_ylabel('Фаза, рад')
    axes[1, 0].set_title('Фаза комплексной огибающей')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_xticks(indices)

    # 4. Вектор Z'
    Z_prime_x = a_m * x_m_noise
    Z_prime_y = a_m * y_m_noise
    Z_prime = np.concatenate([Z_prime_x, Z_prime_y])
    axes[1, 1].stem(np.arange(1, 2*M+1), Z_prime, basefmt=" ", linefmt='b-', markerfmt='bo')
    axes[1, 1].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    axes[1, 1].set_xlabel('Индекс компонента')
    axes[1, 1].set_ylabel('Значение')
    axes[1, 1].set_title("Вектор Z' = [a₁x₁, ..., a_Mx_M, a₁y₁, ..., a_My_M]")
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    print("\n" + "=" * 80)
    print("ГОТОВО! Данные для этапа 2 сохранены")
    print("=" * 80)