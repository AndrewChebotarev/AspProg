import numpy as np
import matplotlib.pyplot as plt


def generate_test_signal_and_channel_response(M, a_m, omega_m, phi_m, N, A_k, tau_k, psi_k, fs, T1, T2, SNR_dB=None):
    """
    Генерация тестового сигнала и отклика многолучевого канала (формулы 1-4 из статьи)

    МАТЕМАТИЧЕСКАЯ МОДЕЛЬ:

    1. Тестовый сигнал (1): S(t) = sum_{m=1}^{M} a_m * cos(omega_m * t - phi_m)

    2. Отклик канала (4):
       s(t) = sum_{m=1}^{M} sum_{k=1}^{N} A_k * a_m * cos(omega_m*(t - tau_k) - phi_m - psi_k)

    3. Наблюдаемый сигнал (2): xi(t) = s(t) + n(t)

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
        - 'S_test': тестовый сигнал S(t)
        - 's_useful': полезный сигнал s(t) (без шума)
        - 'xi_observed': наблюдаемый сигнал xi(t) (с шумом)
        - 'params': параметры модели
        - 'SNR_dB': заданный SNR
    """

    # Временная сетка
    t = np.arange(T1, T2, 1 / fs)
    n_samples = len(t)

    print("=" * 80)
    print("ЭТАП 1: ГЕНЕРАЦИЯ ДАННЫХ")
    print("=" * 80)
    print(f"\nПАРАМЕТРЫ:")
    print(f"  - Интервал: [{T1}, {T2}] с")
    print(f"  - Частота дискретизации: {fs} Гц")
    print(f"  - Отсчетов: {n_samples}")
    print(f"  - Синусоид: M = {M}")
    print(f"  - Лучей: N = {N}")
    print(f"  - SNR: {SNR_dB if SNR_dB else 'без шума'} дБ")

    # Проверка Котельникова
    max_freq = np.max(omega_m) / (2 * np.pi)
    if fs < 2 * max_freq:
        print(f"  ВНИМАНИЕ: fs={fs} < {2 * max_freq:.0f} Гц (теорема Котельникова)")

    print(f"\nТЕСТОВЫЙ СИГНАЛ (формула 1):")
    for m in range(M):
        f_hz = omega_m[m] / (2 * np.pi)
        print(f"  m={m + 1}: a={a_m[m]:.3f}, f={f_hz:.0f} Гц, φ={phi_m[m]:.3f} рад")

    print(f"\nКАНАЛ (лучи):")
    for k in range(N):
        print(f"  k={k + 1}: A={A_k[k]:.3f}, τ={tau_k[k] * 1e6:.0f} мкс, ψ={psi_k[k]:.3f} рад")

    # ===== 1. Тестовый сигнал S(t) =====
    S_test = np.zeros(n_samples)
    for m in range(M):
        S_test += a_m[m] * np.cos(omega_m[m] * t - phi_m[m])

    # ===== 2. Полезный сигнал s(t) =====
    s_useful = np.zeros(n_samples)
    for m in range(M):
        for k in range(N):
            s_useful += A_k[k] * a_m[m] * np.cos(omega_m[m] * (t - tau_k[k]) - phi_m[m] - psi_k[k])

    # ===== 3. Добавление шума =====
    if SNR_dB is not None:
        Ps = np.mean(s_useful ** 2)
        Pn = Ps / (10 ** (SNR_dB / 10))
        noise = np.sqrt(Pn) * np.random.randn(n_samples)
        xi_observed = s_useful + noise
        actual_SNR = 10 * np.log10(Ps / np.mean(noise ** 2))
        print(f"\nЭНЕРГИЯ: Ps={Ps:.4f}, Pn={Pn:.4f}, SNR_факт={actual_SNR:.1f} дБ")
    else:
        xi_observed = s_useful.copy()
        print(f"\nЭНЕРГИЯ: режим без шума")

    """
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))

    n_disp = min(n_samples, 500)
    t_disp = t[:n_disp] * 1e3

    axes[0].plot(t_disp, S_test[:n_disp], 'b-', lw=1.5)
    axes[0].set_ylabel('Амплитуда')
    axes[0].set_title(f'1. Тестовый сигнал S(t) (M={M})')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t_disp, s_useful[:n_disp], 'g-', lw=1.5)
    axes[1].set_ylabel('Амплитуда')
    axes[1].set_title(f'2. Полезный сигнал s(t) (N={N} лучей)')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t_disp, xi_observed[:n_disp], 'r-', lw=1, alpha=0.8, label='С шумом')
    axes[2].plot(t_disp, s_useful[:n_disp], 'g--', lw=0.8, alpha=0.5, label='Без шума')
    axes[2].set_xlabel('Время, мс')
    axes[2].set_ylabel('Амплитуда')
    snr_str = f'SNR={SNR_dB}дБ' if SNR_dB else 'без шума'
    axes[2].set_title(f'3. Наблюдаемый сигнал ξ(t) ({snr_str})')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    
    fig2, ax = plt.subplots(1, 1, figsize=(12, 5))
    freq = np.fft.fftfreq(n_samples, 1 / fs)
    spec = np.abs(np.fft.fft(xi_observed))
    ax.plot(freq[:n_samples // 2], spec[:n_samples // 2], 'b-', lw=1)
    ax.set_xlabel('Частота, Гц')
    ax.set_ylabel('|Спектр|')
    ax.set_title('Спектр наблюдаемого сигнала')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, max_freq * 1.5])

    for m in range(M):
        f_hz = omega_m[m] / (2 * np.pi)
        ax.axvline(x=f_hz, color='r', ls='--', alpha=0.5, label=f'f_{m + 1}' if m == 0 else '')
    if M > 0:
        ax.legend()

    plt.tight_layout()
    plt.show()"""

    result = {
        't': t,
        'S_test': S_test,
        's_useful': s_useful,
        'xi_observed': xi_observed,
        'params': {
            'M': M, 'N': N,
            'a_m': a_m, 'omega_m': omega_m, 'phi_m': phi_m,
            'A_k': A_k, 'tau_k': tau_k, 'psi_k': psi_k,
            'fs': fs, 'T1': T1, 'T2': T2
        },
        'SNR_dB': SNR_dB
    }

    print("\n" + "=" * 80)
    print("ГОТОВО! Данные сохранены в словаре result")
    print("=" * 80)
    print(f"\nКлючи result: {list(result.keys())}")
    print(f"  - t: {len(result['t'])} отсчетов")
    print(f"  - S_test: {result['S_test'].shape}")
    print(f"  - s_useful: {result['s_useful'].shape}")
    #print(f"  - xi_observed: {result['xi_observed'].shape}")=====

    return result


# ========== ИСПОЛЬЗОВАНИЕ ==========
if __name__ == "__main__":
    # ----- ПАРАМЕТРЫ (меняйте здесь) -----

    # Тестовый сигнал
    M = 3
    f_m = np.array([100, 200, 300])  # Гц
    omega_m = 2 * np.pi * f_m
    a_m = np.array([1.0, 0.8, 0.6])
    phi_m = np.array([0, np.pi / 4, np.pi / 3])

    # Канал (лучи)
    N = 2
    A_k = np.array([0.7, 0.5])
    tau_k = np.array([500e-6, 1200e-6])  # 500 и 1200 микросекунд
    psi_k = np.array([0.52, -0.53])

    # Дискретизация
    fs = 20 * np.max(f_m)  # 6000 Гц
    T1 = 0.0
    T2 = 0.01  # 10 мс

    # Шум (None - без шума, или число: 30, 20, 10, 0, -10 дБ)
    SNR_dB = 2

    # --------------------------------------

    print("\n" + "=" * 80)
    print("ГЕНЕРАТОР ТЕСТОВЫХ ДАННЫХ ДЛЯ ОЦЕНКИ ПАРАМЕТРОВ КАНАЛА")
    print("=" * 80)
    print(f"Частота дискретизации: {fs} Гц")
    print(f"Отсчетов: {fs * (T2 - T1):.0f}")
    print(f"SNR: {SNR_dB if SNR_dB else 'без шума'} дБ")
    print("=" * 80 + "\n")

    # Генерируем данные
    data = generate_test_signal_and_channel_response(
        M=M, a_m=a_m, omega_m=omega_m, phi_m=phi_m,
        N=N, A_k=A_k, tau_k=tau_k, psi_k=psi_k,
        fs=fs, T1=T1, T2=T2,
        SNR_dB=SNR_dB
    )