import numpy as np

# Загружаем файл
data = np.load('estimation_results.npz', allow_pickle=True)

print("=" * 80)
print("ПОЛНЫЙ АНАЛИЗ ФАЙЛА estimation_results.npz")
print("=" * 80)

# 1. Количество прогонов из N_est_distribution
n_estimates = data['N_est_distribution']
print(f"\n1. Количество прогонов (из N_est_distribution):")
print(f"   Форма массива: {n_estimates.shape}")
print(f"   Всего прогонов: {n_estimates.shape[0] * n_estimates.shape[1]}")
print(f"   - SNR точек: {n_estimates.shape[0]}")
print(f"   - Прогонов на SNR: {n_estimates.shape[1]}")

# 2. Проверяем другие поля на количество данных
print(f"\n2. Размеры других полей:")
for key in ['tau_mse', 'A_mse', 'psi_mse', 'crb_tau', 'crb_A', 'crb_psi']:
    if key in data:
        print(f"   {key}: {data[key].shape} - {data[key].shape[0]} значений")

# 3. Проверяем, есть ли ошибки (tau_errors и т.д.)
print(f"\n3. Наличие полных данных об ошибках:")
if 'tau_errors' in data:
    tau_errors = data['tau_errors']
    print(f"   tau_errors есть, форма: {tau_errors.shape}")
    total_estimates = 0
    for i, errors in enumerate(tau_errors):
        if len(errors) > 0:
            total_estimates += len(errors)
    print(f"   Всего корректных оценок: {total_estimates}")
else:
    print("   ⚠️ Поля tau_errors, A_errors, psi_errors отсутствуют!")

# 4. Проверяем метаданные
print(f"\n4. Метаданные:")
if 'timestamp' in data:
    print(f"   Время создания: {data['timestamp']}")
else:
    print("   ⚠️ Время создания не сохранено")

# 5. Подробный вывод N_est_distribution
print(f"\n5. Распределение оценок числа лучей (N̂) для каждого SNR:")
snr_db = data['SNR_dB']
n_est_dist = data['N_est_distribution']

for i in range(len(snr_db)):
    estimates = n_est_dist[i]
    unique, counts = np.unique(estimates, return_counts=True)
    print(f"\n   SNR = {snr_db[i]:6.2f} дБ:")
    print(f"      Всего прогонов: {len(estimates)}")
    print(f"      Распределение N̂:")
    for n, cnt in zip(unique, counts):
        print(f"         N̂={n}: {cnt} раз ({cnt/len(estimates)*100:.1f}%)")
    print(f"      P(N̂=2): {np.mean(estimates == 2):.3f}")

# 6. Вывод MSE
print(f"\n6. Среднеквадратичные ошибки:")
print(f"   {'SNR, дБ':<10} {'MSE(τ), мкс²':<18} {'MSE(A)':<12} {'MSE(ψ), рад²':<14}")
print(f"   {'-'*10} {'-'*18} {'-'*12} {'-'*14}")
for i in range(len(snr_db)):
    print(f"   {snr_db[i]:<10.2f} {data['tau_mse'][i]:<18.3e} {data['A_mse'][i]:<12.3e} {data['psi_mse'][i]:<14.3e}")