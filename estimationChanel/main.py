import numpy as np
import matplotlib.pyplot as plt
from channel_estimation import (
    run_single_estimation, generate_signal, t, Fs, T1, T2,
    omegas, a_m, phi_m, tau_true, A_true, psi_true, N_true
)
import warnings
import time
from datetime import datetime
import os

warnings.filterwarnings('ignore')

# ============================================================================
# НАСТРОЙКИ ЭКСПЕРИМЕНТА
# ============================================================================

target_SNR_linear = np.array([0.1, 0.5, 1, 2, 3, 5, 6, 7, 8, 9, 10])
N_MONTE_CARLO = 1000  # Увеличьте для лучшей статистики

RESULTS_DIR = "estimation_results"
os.makedirs(RESULTS_DIR, exist_ok=True)


# ============================================================================
# КЛАСС ДЛЯ ОТСЛЕЖИВАНИЯ ПРОГРЕССА
# ============================================================================
class ProgressTracker:
    def __init__(self, total_snr, total_runs_per_snr):
        self.total_snr = total_snr
        self.total_runs_per_snr = total_runs_per_snr
        self.current_snr_idx = 0
        self.current_run = 0
        self.start_time = None
        self.snr_start_time = None

    def start_experiment(self):
        self.start_time = time.time()
        print("\n" + "=" * 80)
        print(f"🚀 НАЧАЛО ЭКСПЕРИМЕНТА: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 80)

    def start_snr(self, snr_lin, snr_db, n0):
        self.current_snr_idx += 1
        self.snr_start_time = time.time()
        print("\n" + "=" * 80)
        print(f"📊 [SNR {self.current_snr_idx}/{self.total_snr}]")
        print(f"   Линейный SNR: {snr_lin:.2f}")
        print(f"   SNR (дБ):     {snr_db:.1f} дБ")
        print(f"   N0:           {n0:.2e}")
        print(f"   Время начала: {datetime.now().strftime('%H:%M:%S')}")
        print("=" * 80)

    def update_run(self, run, total_runs, tau_est=None, n_est=None):
        self.current_run = run + 1
        percent = (self.current_run / total_runs) * 100
        bar_length = 30
        filled = int(bar_length * self.current_run // total_runs)
        bar = '█' * filled + '░' * (bar_length - filled)
        elapsed = time.time() - self.snr_start_time
        if self.current_run > 1:
            time_per_run = elapsed / self.current_run
            remaining = time_per_run * (total_runs - self.current_run)
            remaining_str = f"{remaining:.0f}с"
        else:
            remaining_str = "оценка..."
        status = f"\r   Прогон [{bar}] {percent:.1f}% ({self.current_run}/{total_runs}) | Осталось: {remaining_str}"
        if tau_est is not None and n_est is not None:
            tau_str = ", ".join([f"{t * 1e3:.2f}" for t in tau_est[:2]]) if len(tau_est) > 0 else "---"
            status += f" | N̂={n_est} | τ̂=[{tau_str}] мс"
        print(status, end="", flush=True)

    def finish_snr(self, results_dict):
        elapsed = time.time() - self.snr_start_time
        print(f"\n   ✅ Завершено за {elapsed:.1f}с")
        print(f"   📈 P(N̂={N_true}) = {results_dict.get('p_correct', 0):.3f}")
        print(f"   📉 MSE(τ) = {results_dict.get('mse_tau', 'N/A')} мкс²")
        print(f"   📉 MSE(A) = {results_dict.get('mse_A', 'N/A')}")
        print(f"   📉 MSE(ψ) = {results_dict.get('mse_psi', 'N/A')} рад²")

    def finish_experiment(self):
        total_elapsed = time.time() - self.start_time
        print("\n" + "=" * 80)
        print(f"✨ ЭКСПЕРИМЕНТ ЗАВЕРШЁН: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"⏱️  Общее время: {total_elapsed:.1f}с ({total_elapsed / 60:.1f} мин)")
        print("=" * 80)


def save_checkpoint(results, checkpoint_file):
    try:
        checkpoint = {
            'SNR_linear': results['SNR_linear'],
            'SNR_dB': results['SNR_dB'],
            'N0_actual': results['N0_actual'],
            'tau_mse': results['tau_mse'],
            'A_mse': results['A_mse'],
            'psi_mse': results['psi_mse'],
            'timestamp': datetime.now().isoformat()
        }
        np.savez(checkpoint_file, **checkpoint)
        print(f"\n   💾 Чекпоинт сохранён: {checkpoint_file}")
    except Exception as e:
        print(f"\n   ⚠️ Не удалось сохранить чекпоинт: {e}")


def calculate_crb(snr_linear, tau_true, A_true, omegas, a_m, Fs, T1, T2):
    freqs = omegas / (2 * np.pi)
    mean_f = np.mean(freqs)
    f_rms_sq = np.mean((freqs - mean_f) ** 2)
    T_len = T2 - T1
    crb_tau_single = 1 / (2 * snr_linear * (2 * np.pi) ** 2 * f_rms_sq * T_len)
    crb_tau = crb_tau_single * np.ones(len(tau_true))
    crb_A = 1 / (2 * snr_linear) * np.ones(len(A_true))
    crb_psi = 1 / (2 * snr_linear) * np.ones(len(psi_true))
    return crb_tau, crb_A, crb_psi


def get_N0_for_SNR(target_snr_linear, signal_power_ref, Fs):
    return 2 * signal_power_ref / (target_snr_linear * Fs)


# ============================================================================
# ФУНКЦИИ ВИЗУАЛИЗАЦИИ (сокращены для краткости)
# ============================================================================
def plot_all_results(results, crb_tau_values, crb_A_values, crb_psi_values, N_true):
    """Построение всех графиков"""
    snr_db = np.array(results['SNR_dB'])

    # Создаём фигуру с 6 подграфиками
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Анализ точности МП-оценок параметров многолучевого канала', fontsize=14, fontweight='bold')

    # 1. Вероятность правильного определения числа лучей
    ax = axes[0, 0]
    p_correct = [np.mean(np.array(dist) == N_true) for dist in results['N_est_distribution']]
    ax.plot(snr_db, p_correct, 'b-o', linewidth=2, markersize=8)
    ax.fill_between(snr_db, 0, p_correct, alpha=0.3)
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('P(N̂ = N_true)')
    ax.set_title('Вероятность правильного определения\nчисла лучей')
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1.05])

    # 2. MSE задержек
    ax = axes[0, 1]
    valid = [not np.isnan(x) and x > 0 for x in results['tau_mse']]
    if any(valid):
        valid_snr = snr_db[valid]
        valid_mse = np.array(results['tau_mse'])[valid]
        valid_crb = np.array(crb_tau_values)[valid] / 1e6
        ax.loglog(valid_snr, valid_mse, 'r-s', linewidth=2, markersize=8, label='MSE(τ)')
        ax.loglog(valid_snr, valid_crb, 'k--', linewidth=2, label='CRB(τ)')
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('MSE(τ), мкс²')
    ax.set_title('Ошибка оценки задержек лучей')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. MSE амплитуд
    ax = axes[0, 2]
    valid = [not np.isnan(x) and x > 0 for x in results['A_mse']]
    if any(valid):
        valid_snr = snr_db[valid]
        valid_mse = np.array(results['A_mse'])[valid]
        valid_crb = np.array(crb_A_values)[valid]
        ax.loglog(valid_snr, valid_mse, 'g-d', linewidth=2, markersize=8, label='MSE(A)')
        ax.loglog(valid_snr, valid_crb, 'k--', linewidth=2, label='CRB(A)')
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('MSE(A)')
    ax.set_title('Ошибка оценки амплитуд лучей')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. MSE фаз
    ax = axes[1, 0]
    valid = [not np.isnan(x) and x > 0 for x in results['psi_mse']]
    if any(valid):
        valid_snr = snr_db[valid]
        valid_mse = np.array(results['psi_mse'])[valid]
        valid_crb = np.array(crb_psi_values)[valid]
        ax.loglog(valid_snr, valid_mse, 'm-o', linewidth=2, markersize=8, label='MSE(ψ)')
        ax.loglog(valid_snr, valid_crb, 'k--', linewidth=2, label='CRB(ψ)')
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('MSE(ψ), рад²')
    ax.set_title('Ошибка оценки фаз лучей')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 5. Относительная эффективность
    ax = axes[1, 1]
    valid = [not np.isnan(x) and x > 0 for x in results['tau_mse']]
    if any(valid):
        valid_snr = snr_db[valid]
        valid_mse = np.array(results['tau_mse'])[valid]
        valid_crb = np.array(crb_tau_values)[valid] / 1e6
        ratio = valid_mse / valid_crb
        ax.semilogy(valid_snr, ratio, 'purple', linewidth=2, marker='o')
        ax.axhline(1, color='r', linestyle='--', label='CRB')
        ax.axhline(10, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('SNR, дБ')
    ax.set_ylabel('MSE / CRB')
    ax.set_title('Относительная эффективность оценки\n(1 = оптимальная оценка)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 6. Гистограмма распределения N̂
    ax = axes[1, 2]
    mid_idx = len(snr_db) // 2
    N_est_mid = results['N_est_distribution'][mid_idx]
    unique_n, counts = np.unique(N_est_mid, return_counts=True)
    ax.bar(unique_n, counts / len(N_est_mid) * 100, color='steelblue', edgecolor='black')
    ax.set_xlabel('Оценка числа лучей N̂')
    ax.set_ylabel('Частота, %')
    ax.set_title(f'Распределение N̂ при SNR = {snr_db[mid_idx]:.1f} дБ')
    ax.set_xticks(unique_n)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'estimation_analysis.png'), dpi=300, bbox_inches='tight')
    plt.show()
    print(f"   ✅ Графики сохранены")


# ============================================================================
# ОСНОВНАЯ ФУНКЦИЯ
# ============================================================================
def main():
    tracker = ProgressTracker(len(target_SNR_linear), N_MONTE_CARLO)
    tracker.start_experiment()

    # Калибровка
    print("📡 Калибровка мощности сигнала...")
    s_ref, _ = generate_signal(t, omegas, a_m, phi_m, tau_true, A_true, psi_true, 1e-4, Fs)
    signal_power_ref = np.mean(s_ref ** 2)
    print(f"   Мощность сигнала: {signal_power_ref:.4f}")

    # Инициализация результатов
    results = {
        'SNR_linear': [], 'SNR_dB': [], 'N0_actual': [],
        'N_est_distribution': [],
        'tau_errors': [], 'A_errors': [], 'psi_errors': [],
        'tau_mse': [], 'A_mse': [], 'psi_mse': [],
    }

    # Параметры для отбрасывания выбросов
    max_tau_error_us = 200  # Максимальная разумная ошибка задержки (мкс)
    tau_spacing_us = np.diff(tau_true)[0] * 1e6  # 400 мкс

    for snr_idx, target_snr_lin in enumerate(target_SNR_linear):
        target_snr_db = 10 * np.log10(target_snr_lin)
        N0_target = get_N0_for_SNR(target_snr_lin, signal_power_ref, Fs)

        tracker.start_snr(target_snr_lin, target_snr_db, N0_target)

        tau_errors_all = []
        A_errors_all = []
        psi_errors_all = []
        N_estimates = []
        run_times = []

        for run in range(N_MONTE_CARLO):
            run_start = time.time()

            try:
                result = run_single_estimation(N0_target)
                run_times.append(time.time() - run_start)
                N_estimates.append(result['N_est'])

                if result['N_est'] == N_true:
                    # Сортируем оценки по возрастанию
                    tau_est_sorted = np.sort(result['tau_est'])

                    for k in range(N_true):
                        # Вычисляем ошибку
                        tau_error = (tau_est_sorted[k] - tau_true[k]) * 1e6
                        A_error = result['A_est'][k] - result['A_true'][k]
                        psi_error = result['psi_est'][k] - result['psi_true'][k]

                        # ====================================================
                        # ИСПРАВЛЕННАЯ ОБРАБОТКА ВЫБРОСОВ
                        # ====================================================

                        # Логическая проверка: если ошибка слишком большая - это выброс
                        if abs(tau_error) > max_tau_error_us:
                            # Пытаемся определить, что пошло не так
                            if k == 0 and tau_est_sorted[0] * 1e6 > tau_spacing_us / 2:
                                # Первый луч далеко от 0 - возможно, лучи перепутаны
                                if len(tau_est_sorted) > 1:
                                    # Меняем местами оценки
                                    tau_est_swapped = np.array([tau_est_sorted[1], tau_est_sorted[0]])
                                    tau_error = (tau_est_swapped[k] - tau_true[k]) * 1e6
                                    if abs(tau_error) <= max_tau_error_us:
                                        # Успешно переставили
                                        continue

                            # Если коррекция не помогла - пропускаем эту ошибку
                            # Не добавляем выброс в статистику
                            continue

                        # Нормальная ошибка - добавляем
                        tau_errors_all.append(tau_error)
                        A_errors_all.append(A_error)
                        psi_errors_all.append(psi_error)

            except Exception as e:
                print(f"\n   ⚠️ Ошибка в прогоне {run}: {e}")
                N_estimates.append(0)
                continue

            tracker.update_run(run, N_MONTE_CARLO,
                               tau_est=result.get('tau_est', None),
                               n_est=result.get('N_est', None))

        print()  # Новая строка

        # Преобразуем в массивы
        tau_errors_all = np.array(tau_errors_all) if tau_errors_all else np.array([])
        A_errors_all = np.array(A_errors_all) if A_errors_all else np.array([])
        psi_errors_all = np.array(psi_errors_all) if psi_errors_all else np.array([])

        # Сохраняем результаты
        results['SNR_linear'].append(target_snr_lin)
        results['SNR_dB'].append(target_snr_db)
        results['N0_actual'].append(N0_target)
        results['N_est_distribution'].append(N_estimates)

        # Расчёт MSE (только если есть данные)
        if len(tau_errors_all) > 0:
            results['tau_errors'].append(tau_errors_all)
            results['A_errors'].append(A_errors_all)
            results['psi_errors'].append(psi_errors_all)
            results['tau_mse'].append(np.mean(tau_errors_all ** 2))
            results['A_mse'].append(np.mean(A_errors_all ** 2))
            results['psi_mse'].append(np.mean(psi_errors_all ** 2))
        else:
            results['tau_errors'].append([])
            results['A_errors'].append([])
            results['psi_errors'].append([])
            results['tau_mse'].append(np.nan)
            results['A_mse'].append(np.nan)
            results['psi_mse'].append(np.nan)

        # Статистика
        p_correct = np.mean(np.array(N_estimates) == N_true)
        mse_tau = results['tau_mse'][-1] if not np.isnan(results['tau_mse'][-1]) else np.nan
        mse_A = results['A_mse'][-1] if not np.isnan(results['A_mse'][-1]) else np.nan
        mse_psi = results['psi_mse'][-1] if not np.isnan(results['psi_mse'][-1]) else np.nan
        avg_time = np.mean(run_times) if run_times else 0

        tracker.finish_snr({
            'p_correct': p_correct,
            'mse_tau': mse_tau,
            'mse_A': mse_A,
            'mse_psi': mse_psi
        })

        print(f"   ⏱️  Среднее время: {avg_time:.2f}с")
        print(f"   📊 Валидных выборок: {len(tau_errors_all)}/{N_MONTE_CARLO * N_true}")
        print(f"   🧹 Отброшено выбросов: {N_MONTE_CARLO * N_true - len(tau_errors_all)}")

        save_checkpoint(results, os.path.join(RESULTS_DIR, f'checkpoint_snr_{snr_idx + 1}.npz'))

    tracker.finish_experiment()

    # Расчёт CRB
    print("\n" + "=" * 80)
    print("📐 РАСЧЁТ ГРАНИЦЫ КРАМЕРА-РАО")
    print("=" * 80)

    crb_tau_values = []
    crb_A_values = []
    crb_psi_values = []

    for snr_lin in results['SNR_linear']:
        crb_tau, crb_A, crb_psi = calculate_crb(snr_lin, tau_true, A_true, omegas, a_m, Fs, T1, T2)
        crb_tau_values.append(np.mean(crb_tau) * 1e12)
        crb_A_values.append(np.mean(crb_A))
        crb_psi_values.append(np.mean(crb_psi))

    # Визуализация
    print("\n🎨 Построение графиков...")
    plot_all_results(results, crb_tau_values, crb_A_values, crb_psi_values, N_true)

    # Вывод таблицы
    print("\n" + "=" * 80)
    print("📊 СВОДНАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ")
    print("=" * 80)
    print(
        f"{'SNR лин.':<10} {'SNR, дБ':<10} {'N0':<12} {'P(N̂=2)':<10} {'MSE(τ), мкс²':<18} {'MSE(A)':<12} {'MSE(ψ), рад²':<14}")
    print("-" * 90)

    for i in range(len(results['SNR_linear'])):
        snr_lin = results['SNR_linear'][i]
        snr_db = results['SNR_dB'][i]
        p_corr = np.mean(np.array(results['N_est_distribution'][i]) == N_true)
        mse_tau = results['tau_mse'][i] if not np.isnan(results['tau_mse'][i]) else np.nan
        mse_A = results['A_mse'][i] if not np.isnan(results['A_mse'][i]) else np.nan
        mse_psi = results['psi_mse'][i] if not np.isnan(results['psi_mse'][i]) else np.nan

        # Заменяем NaN и 0 на разумные значения для отображения
        mse_tau_display = mse_tau if not np.isnan(mse_tau) and mse_tau > 0 else np.nan

        print(
            f"{snr_lin:<10.2f} {snr_db:<10.1f} {results['N0_actual'][i]:<12.2e} {p_corr:<10.3f} {mse_tau_display:<18.3e} {mse_A:<12.3e} {mse_psi:<14.3e}")

    print("=" * 90)

    # Сохранение
    np.savez(os.path.join(RESULTS_DIR, 'estimation_results.npz'),
             SNR_linear=results['SNR_linear'],
             SNR_dB=results['SNR_dB'],
             N0_values=results['N0_actual'],
             N_est_distribution=np.array(results['N_est_distribution'], dtype=object),
             tau_mse=results['tau_mse'],
             A_mse=results['A_mse'],
             psi_mse=results['psi_mse'],
             crb_tau=crb_tau_values,
             crb_A=crb_A_values,
             crb_psi=crb_psi_values)

    print(f"\n💾 Результаты сохранены в '{RESULTS_DIR}/'")
    print("✨ Анализ завершён!")


if __name__ == "__main__":
    main()