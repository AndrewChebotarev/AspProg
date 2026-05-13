# ============================================================================
# МОНТЕ-КАРЛО АНАЛИЗ МП-ОЦЕНИВАНИЯ ПАРАМЕТРОВ КАНАЛА
# Запуск основного алгоритма при разных ОСШ с накоплением статистики
# ============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.linalg import sqrtm
import pickle
from datetime import datetime
import os
import warnings

warnings.filterwarnings('ignore')

# Импортируем функции из вашего основного скрипта
# (предполагается, что он называется channel_estimation1.py)
from channel_estimation import (
    generate_signal, compute_sufficient_statistics,
    build_Q_matrix, build_H_matrix, calc_L, search_N,
    t, omegas, a_m, phi_m, tau_true, A_true, psi_true,
    Fs, T1, T2, M, N_true
)


# ============================================================================
# КЛАСС ДЛЯ СБОРА СТАТИСТИКИ
# ============================================================================
class MonteCarloStats:
    """Сбор и анализ статистики оценок"""

    def __init__(self, snr_values, n_realizations=1000):
        self.snr_values = snr_values
        self.n_realizations = n_realizations
        self.results = {snr: {
            'N_est': [],  # Оценки числа лучей
            'tau_est': [],  # Оценки задержек
            'A_est': [],  # Оценки амплитуд
            'psi_est': [],  # Оценки фаз
            'tau_errors': [],  # Ошибки задержек
            'A_errors': [],  # Ошибки амплитуд
            'psi_errors': [],  # Ошибки фаз
            'L_max': []  # Макс. значение ФОП
        } for snr in snr_values}

    def add_result(self, snr, N_est, tau_est, A_est, psi_est, L_max):
        """Добавление результатов одной реализации"""
        # Сортируем лучи по задержке для сопоставления с истинными
        sort_idx = np.argsort(tau_est)
        tau_est = tau_est[sort_idx]
        A_est = A_est[sort_idx]
        psi_est = psi_est[sort_idx]

        self.results[snr]['N_est'].append(N_est)
        self.results[snr]['tau_est'].append(tau_est)
        self.results[snr]['A_est'].append(A_est)
        self.results[snr]['psi_est'].append(psi_est)
        self.results[snr]['L_max'].append(L_max)

        # Вычисляем ошибки (только если N_est совпадает с истинным)
        if N_est == N_true:
            tau_errors = tau_est - tau_true[:N_est]
            A_errors = A_est - A_true[:N_est]
            psi_errors = psi_est - psi_true[:N_est]
        else:
            # Если число лучей оценено неверно, ошибка большая
            tau_errors = np.full(N_true, np.nan)
            A_errors = np.full(N_true, np.nan)
            psi_errors = np.full(N_true, np.nan)

        self.results[snr]['tau_errors'].append(tau_errors)
        self.results[snr]['A_errors'].append(A_errors)
        self.results[snr]['psi_errors'].append(psi_errors)

    def compute_statistics(self):
        """Вычисление статистик по всем реализациям"""
        stats = {}

        for snr in self.snr_values:
            n_real = len(self.results[snr]['N_est'])

            # Вероятность правильного определения числа лучей
            prob_correct = np.mean(np.array(self.results[snr]['N_est']) == N_true)

            # Собираем только реализации с правильным N
            valid_idx = [i for i, N in enumerate(self.results[snr]['N_est']) if N == N_true]
            n_valid = len(valid_idx)

            if n_valid > 0:
                tau_errors_valid = np.array([self.results[snr]['tau_errors'][i] for i in valid_idx])
                A_errors_valid = np.array([self.results[snr]['A_errors'][i] for i in valid_idx])
                psi_errors_valid = np.array([self.results[snr]['psi_errors'][i] for i in valid_idx])

                # MSE для задержек, амплитуд и фаз
                mse_tau = np.nanmean(tau_errors_valid ** 2, axis=0)
                mse_A = np.nanmean(A_errors_valid ** 2, axis=0)
                mse_psi = np.nanmean(psi_errors_valid ** 2, axis=0)

                # Смещение (bias)
                bias_tau = np.nanmean(tau_errors_valid, axis=0)
                bias_A = np.nanmean(A_errors_valid, axis=0)
                bias_psi = np.nanmean(psi_errors_valid, axis=0)

                # Стандартное отклонение
                std_tau = np.nanstd(tau_errors_valid, axis=0)
                std_A = np.nanstd(A_errors_valid, axis=0)
                std_psi = np.nanstd(psi_errors_valid, axis=0)
            else:
                mse_tau = np.full(N_true, np.nan)
                mse_A = np.full(N_true, np.nan)
                mse_psi = np.full(N_true, np.nan)
                bias_tau = np.full(N_true, np.nan)
                bias_A = np.full(N_true, np.nan)
                bias_psi = np.full(N_true, np.nan)
                std_tau = np.full(N_true, np.nan)
                std_A = np.full(N_true, np.nan)
                std_psi = np.full(N_true, np.nan)

            stats[snr] = {
                'prob_correct': prob_correct,
                'n_valid': n_valid,
                'mse_tau': mse_tau,
                'mse_A': mse_A,
                'mse_psi': mse_psi,
                'bias_tau': bias_tau,
                'bias_A': bias_A,
                'bias_psi': bias_psi,
                'std_tau': std_tau,
                'std_A': std_A,
                'std_psi': std_psi
            }

        return stats


# ============================================================================
# ФУНКЦИЯ ДЛЯ ВЫЧИСЛЕНИЯ ГРАНИЦЫ КРАМЕРА-РАО
# ============================================================================
def compute_cramer_rao_bound(snr_db, tau_true, A_true, psi_true, omegas, a_m, T, Fs, N0):
    """
    Вычисление границы Крамера-Рао для оценок задержек

    Для многолучевого канала с гармоническим тестовым сигналом
    """
    snr_linear = 10 ** (snr_db / 10)
    M = len(omegas)
    N = len(tau_true)

    # Информационная матрица Фишера для задержек (упрощённая оценка)
    # При ортогональных сигналах элементы матрицы:
    # F_{kl} ≈ (2 * SNR * M) * sinc(2*pi*f_max*(tau_k - tau_l))

    f_max = max(omegas) / (2 * np.pi)  # Максимальная частота
    FIM = np.zeros((N, N))

    for k in range(N):
        for l in range(N):
            delta_tau = tau_true[k] - tau_true[l]
            if k == l:
                # Дисперсия оценки задержки для одного луча
                # CRB_tau = 1 / (8 * pi^2 * SNR * M * f_rms^2 * T)
                FIM[k, k] = 8 * np.pi ** 2 * snr_linear * M * (f_max ** 2) * T
            else:
                # Взаимная информация между лучами
                FIM[k, l] = FIM[k, k] * np.sinc(2 * f_max * delta_tau)

    # Граница Крамера-Рао - обратная матрица Фишера
    try:
        CRB = np.linalg.inv(FIM)
        crb_tau = np.sqrt(np.diag(CRB))  # СКО
    except:
        crb_tau = np.ones(N) * np.nan

    return crb_tau


# ============================================================================
# ОСНОВНАЯ ФУНКЦИЯ МОНТЕ-КАРЛО
# ============================================================================
def run_monte_carlo(snr_linear_values, n_realizations=100, N_max=3, save_results=True, signal_power=None):
    """
    Запуск Монте-Карло моделирования с линейными ОСШ

    Параметры:
    - snr_linear_values: список линейных ОСШ [0, 0.5, 1, ...]
    - signal_power: мощность сигнала (для расчёта шума)
    """
    print("=" * 80)
    print(f"🎲 МОНТЕ-КАРЛО МОДЕЛИРОВАНИЕ")
    print(f"   ЧИСЛО ЛУЧЕЙ: {N_true}")
    print(f"   ОСШ (линейные): {snr_linear_values}")
    print(f"   РЕАЛИЗАЦИЙ: {n_realizations}")
    print("=" * 80)

    stats_collector = MonteCarloStats(snr_linear_values, n_realizations)

    # Основной цикл по ОСШ
    for snr_idx, snr_linear in enumerate(snr_linear_values):
        print(
            f"\n🔍 ОСШ = {snr_linear} (≈ {10 * np.log10(snr_linear) if snr_linear > 0 else -np.inf:.1f} дБ) [{snr_idx + 1}/{len(snr_linear_values)}]")

        # Цикл по реализациям
        for run in range(n_realizations):
            if (run + 1) % 100 == 0:
                print(f"  Выполнено {run + 1}/{n_realizations} реализаций...")

            # 🔧 Генерация данных с текущим линейным SNR
            s_true, xi = generate_signal(
                t, omegas, a_m, phi_m,
                tau_true, A_true, psi_true,
                snr_linear=snr_linear,
                signal_power=signal_power
            )

            # Вычисление достаточной статистики
            # 🔧 Для compute_sufficient_statistics нужно вычислить N0 из SNR
            if snr_linear > 0:
                N0_current = signal_power * 2 / (snr_linear * Fs)
            else:
                N0_current = 1e-2  # Запасной вариант

            Z_prime = compute_sufficient_statistics(xi, t, omegas, a_m, N0_current)

            # Поиск оценок для разных N
            L_N = np.zeros(N_max)
            tau_est_all = []

            for N_guess in range(1, N_max + 1):
                tau_est, L_val = search_N(N_guess, Z_prime, omegas, a_m, phi_m, N0_current, T2 - T1)
                L_N[N_guess - 1] = L_val
                tau_est_all.append(tau_est)

            # Выбор числа лучей по критерию
            criterion = (L_N ** 2) / np.arange(1, N_max + 1)
            N_est = np.argmax(criterion) + 1

            # Восстановление амплитуд и фаз
            best_tau = np.array(tau_est_all[N_est - 1])

            if len(best_tau) > 0:
                Q = build_Q_matrix(best_tau, omegas, a_m, N0_current, T2 - T1)
                H = build_H_matrix(best_tau, omegas, phi_m)

                try:
                    v = np.linalg.solve(Q + 1e-10 * np.eye(2 * len(best_tau)), H.T @ Z_prime)
                    Ac, As = v[:len(best_tau)], v[len(best_tau):]
                    A_est = np.sqrt(Ac ** 2 + As ** 2)
                    psi_est = np.arctan2(As, Ac)
                except:
                    A_est = np.zeros(len(best_tau))
                    psi_est = np.zeros(len(best_tau))
            else:
                A_est = np.array([])
                psi_est = np.array([])

            # Сохранение результатов
            stats_collector.add_result(snr_linear, N_est, best_tau, A_est, psi_est, L_N[N_est - 1])

    # Вычисление статистик
    print("\n📈 Вычисление статистик...")
    statistics = stats_collector.compute_statistics()

    # Вывод результатов (с линейными SNR)
    print("\n" + "=" * 80)
    print("📊 РЕЗУЛЬТАТЫ МОНТЕ-КАРЛО МОДЕЛИРОВАНИЯ")
    print("=" * 80)
    print(f"{'SNR (лин)':<12} {'SNR, дБ':<10} {'P(N̂=N)':<12} {'MSE_τ₁, нс²':<15} {'MSE_τ₂, нс²':<15}")
    print("-" * 80)

    for snr in snr_linear_values:
        snr_db = 10 * np.log10(snr) if snr > 0 else -np.inf
        prob = statistics[snr]['prob_correct']
        mse_tau = statistics[snr]['mse_tau']

        if not np.isnan(mse_tau[0]):
            print(f"{snr:<12.2f} {snr_db:<10.1f} {prob:<12.3f} {mse_tau[0] * 1e12:<15.2f} {mse_tau[1] * 1e12:<15.2f}")
        else:
            print(f"{snr:<12.2f} {snr_db:<10.1f} {prob:<12.3f} {'N/A':<15} {'N/A':<15}")

    # Сохранение результатов
    if save_results:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"monte_carlo_results_linear_SNR_{timestamp}.pkl"
        with open(filename, 'wb') as f:
            pickle.dump({
                'stats_collector': stats_collector,
                'statistics': statistics,
                'snr_linear_values': snr_linear_values,
                'n_realizations': n_realizations
            }, f)
        print(f"\n💾 Результаты сохранены в {filename}")

    return stats_collector, statistics


# ============================================================================
# ФУНКЦИИ ДЛЯ ВИЗУАЛИЗАЦИИ
# ============================================================================
def plot_results(stats_collector, statistics, snr_values, n_realizations):
    """Построение всех графиков"""

    fig = plt.figure(figsize=(16, 12))

    # 1. Вероятность правильного определения числа лучей
    ax1 = plt.subplot(3, 3, 1)
    prob_correct = [statistics[snr]['prob_correct'] for snr in snr_values]
    ax1.semilogy(snr_values, prob_correct, 'bo-', linewidth=2, markersize=8)
    ax1.set_xlabel('ОСШ, дБ')
    ax1.set_ylabel('P(N̂ = N)')
    ax1.set_title('Вероятность правильного определения\nчисла лучей')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([0, 1.05])

    # 2. MSE задержек (первый луч)
    ax2 = plt.subplot(3, 3, 2)
    mse_tau1 = [statistics[snr]['mse_tau'][0] for snr in snr_values]
    ax2.loglog(snr_values, mse_tau1, 'ro-', linewidth=2, markersize=8, label='MSE τ₁ (Монте-Карло)')

    # Добавляем границу Крамера-Рао
    crb_tau1 = []
    for snr in snr_values:
        crb = compute_cramer_rao_bound(snr, tau_true, A_true, psi_true, omegas, a_m, T2 - T1, Fs, None)
        crb_tau1.append(crb[0] ** 2)
    ax2.loglog(snr_values, crb_tau1, 'b--', linewidth=2, label='Граница Крамера-Рао')

    ax2.set_xlabel('ОСШ, дБ')
    ax2.set_ylabel('MSE, с²')
    ax2.set_title('Среднеквадратичная ошибка\nоценки задержки τ₁')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. MSE задержек (второй луч)
    ax3 = plt.subplot(3, 3, 3)
    mse_tau2 = [statistics[snr]['mse_tau'][1] for snr in snr_values]
    ax3.loglog(snr_values, mse_tau2, 'ro-', linewidth=2, markersize=8, label='MSE τ₂ (Монте-Карло)')

    crb_tau2 = []
    for snr in snr_values:
        crb = compute_cramer_rao_bound(snr, tau_true, A_true, psi_true, omegas, a_m, T2 - T1, Fs, None)
        crb_tau2.append(crb[1] ** 2)
    ax3.loglog(snr_values, crb_tau2, 'b--', linewidth=2, label='Граница Крамера-Рао')

    ax3.set_xlabel('ОСШ, дБ')
    ax3.set_ylabel('MSE, с²')
    ax3.set_title('Среднеквадратичная ошибка\nоценки задержки τ₂')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. MSE амплитуд
    ax4 = plt.subplot(3, 3, 4)
    mse_A1 = [statistics[snr]['mse_A'][0] for snr in snr_values]
    mse_A2 = [statistics[snr]['mse_A'][1] for snr in snr_values]
    ax4.loglog(snr_values, mse_A1, 'ro-', linewidth=2, markersize=8, label='MSE A₁')
    ax4.loglog(snr_values, mse_A2, 'bs-', linewidth=2, markersize=8, label='MSE A₂')
    ax4.set_xlabel('ОСШ, дБ')
    ax4.set_ylabel('MSE')
    ax4.set_title('Среднеквадратичная ошибка\nоценки амплитуд')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Гистограмма оценок задержек при низком ОСШ
    ax5 = plt.subplot(3, 3, 5)
    low_snr = snr_values[0]
    tau_ests = [est[0] * 1e6 for est in stats_collector.results[low_snr]['tau_est']
                if len(est) > 0 and stats_collector.results[low_snr]['N_est'][
                    stats_collector.results[low_snr]['tau_est'].index(est)] == N_true]
    if len(tau_ests) > 0:
        ax5.hist(tau_ests, bins=30, alpha=0.7, color='blue', edgecolor='black')
        ax5.axvline(tau_true[0] * 1e6, color='red', linestyle='--', linewidth=2, label='Истинное значение')
        ax5.set_xlabel('τ₁, мкс')
        ax5.set_ylabel('Частота')
        ax5.set_title(f'Распределение оценок τ₁\nОСШ = {low_snr} дБ')
        ax5.legend()
        ax5.grid(True, alpha=0.3)

    # 6. Диаграмма рассеяния (оценки vs истинные значения)
    ax6 = plt.subplot(3, 3, 6)
    for snr in snr_values[::2]:  # Берём каждый второй для наглядности
        tau_est_all = []
        tau_true_all = []
        for i, est in enumerate(stats_collector.results[snr]['tau_est']):
            if len(est) > 0 and stats_collector.results[snr]['N_est'][i] == N_true:
                tau_est_all.extend(est[:2] * 1e6)
                tau_true_all.extend(tau_true[:2] * 1e6)

        if len(tau_est_all) > 0:
            ax6.scatter(tau_true_all, tau_est_all, alpha=0.3, s=10, label=f'{snr} дБ')

    ax6.plot([0, 500], [0, 500], 'k--', linewidth=1, label='Идеальная оценка')
    ax6.set_xlabel('Истинная задержка, мкс')
    ax6.set_ylabel('Оценённая задержка, мкс')
    ax6.set_title('Диаграмма рассеяния оценок задержек')
    ax6.legend(fontsize=8)
    ax6.grid(True, alpha=0.3)

    # 7. Смещение (bias) оценок
    ax7 = plt.subplot(3, 3, 7)
    bias_tau1 = [statistics[snr]['bias_tau'][0] * 1e9 for snr in snr_values]
    bias_tau2 = [statistics[snr]['bias_tau'][1] * 1e9 for snr in snr_values]
    ax7.semilogy(snr_values, np.abs(bias_tau1), 'ro-', linewidth=2, markersize=8, label='|Bias τ₁|, нс')
    ax7.semilogy(snr_values, np.abs(bias_tau2), 'bs-', linewidth=2, markersize=8, label='|Bias τ₂|, нс')
    ax7.set_xlabel('ОСШ, дБ')
    ax7.set_ylabel('|Смещение|, нс')
    ax7.set_title('Абсолютное смещение оценок задержек')
    ax7.legend()
    ax7.grid(True, alpha=0.3)

    # 8. Эффективность оценок (отношение CRB к MSE)
    ax8 = plt.subplot(3, 3, 8)
    efficiency1 = []
    efficiency2 = []
    for i, snr in enumerate(snr_values):
        if mse_tau1[i] > 0:
            efficiency1.append(crb_tau1[i] / mse_tau1[i])
            efficiency2.append(crb_tau2[i] / mse_tau2[i])
        else:
            efficiency1.append(np.nan)
            efficiency2.append(np.nan)

    ax8.semilogy(snr_values, efficiency1, 'ro-', linewidth=2, markersize=8, label='Эффективность τ₁')
    ax8.semilogy(snr_values, efficiency2, 'bs-', linewidth=2, markersize=8, label='Эффективность τ₂')
    ax8.axhline(y=1, color='k', linestyle='--', linewidth=1, label='Граница Крамера-Рао')
    ax8.set_xlabel('ОСШ, дБ')
    ax8.set_ylabel('CRB / MSE')
    ax8.set_title('Эффективность оценок\n(CRB/MSE → 1 при высоком ОСШ)')
    ax8.legend()
    ax8.grid(True, alpha=0.3)

    # 9. Количество валидных реализаций
    ax9 = plt.subplot(3, 3, 9)
    n_valid = [statistics[snr]['n_valid'] for snr in snr_values]
    ax9.bar(range(len(snr_values)), n_valid, alpha=0.7, color='green', edgecolor='black')
    ax9.set_xticks(range(len(snr_values)))
    ax9.set_xticklabels([f'{snr}' for snr in snr_values])
    ax9.set_xlabel('ОСШ, дБ')
    ax9.set_ylabel('Количество')
    ax9.set_title('Число реализаций с\nправильным определением N')
    ax9.grid(True, alpha=0.3, axis='y')

    plt.suptitle(f'Монте-Карло анализ МП-оценивания (N={N_true} лучей, {n_realizations} реализаций)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('monte_carlo_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()


# ============================================================================
# ЗАПУСК АНАЛИЗА (исправлено для линейных SNR)
# ============================================================================
if __name__ == "__main__":
    # 🔧 Линейные значения ОСШ (как вы хотели)
    snr_linear_values = [0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 10]

    # 🔧 1000 реализаций для статистической достоверности
    n_realizations = 1000

    # Предварительный расчёт мощности сигнала (один раз)
    s_pilot, _ = generate_signal(t, omegas, a_m, phi_m, tau_true, A_true, psi_true,
                                 snr_linear=1e10, signal_power=None)  # Очень высокий SNR ≈ без шума
    signal_power = np.mean(s_pilot ** 2)

    print("\n" + "=" * 80)
    print("🚀 ЗАПУСК МОНТЕ-КАРЛО АНАЛИЗА")
    print("=" * 80)
    print(f"📊 Параметры:")
    print(f"   - Истинное число лучей: {N_true}")
    print(f"   - Задержки: {tau_true * 1e3} мс")
    print(f"   - Амплитуды: {A_true}")
    print(f"   - ОСШ (линейные): {snr_linear_values}")
    print(f"   - Мощность сигнала: {signal_power:.4f}")
    print(f"   - Реализаций на ОСШ: {n_realizations}")
    print(f"   - Всего вычислений: {len(snr_linear_values) * n_realizations}")
    print("=" * 80 + "\n")

    # 🔧 Запуск Монте-Карло с передачей signal_power
    stats_collector, statistics = run_monte_carlo(
        snr_linear_values,
        n_realizations,
        N_max=3,
        signal_power=signal_power  # 🔧 Новый параметр
    )

    # Построение графиков (используем линейные SNR для подписей)
    plot_results(stats_collector, statistics, snr_linear_values, n_realizations, use_linear_snr=True)

    print("\n✅ Анализ завершён!")