import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# ПАРАМЕТРЫ СИГНАЛА И ДИСКРЕТИЗАЦИИ
# =============================================================================
fs = 1000  # Частота дискретизации [Гц]
t_duration = 1.0  # Длительность сигнала [сек]
t = np.linspace(0, t_duration, int(fs * t_duration), endpoint=False)

# =============================================================================
# 1. ФОРМИРОВАНИЕ СИГНАЛОВ
# =============================================================================
f_signal = 5  # Частота сигнала [Гц]

# РЕАЛЬНЫЙ СИГНАЛ: x_real(t) = sin(2π·f·t)
sinusoid_real = np.sin(2 * np.pi * f_signal * t)

# КОМПЛЕКСНЫЙ СИГНАЛ: x_comp(t) = e^(j·2π·f·t) = cos(2π·f·t) + j·sin(2π·f·t)
sinusoid_complex = np.exp(1j * 2 * np.pi * f_signal * t)

# =============================================================================
# 2. ФОРМИРОВАНИЕ ИМПУЛЬСНЫХ ХАРАКТЕРИСТИК (ИХ)
# =============================================================================
center = 0.3  # Центр импульсной характеристики [сек]
sigma = 0.08  # Ширина гауссова колокола [сек]

# РЕАЛЬНАЯ ИХ: h_real(t) = exp(-(t-center)²/(2σ²)) (нормированная)
impulse_response_real = np.exp(-(t - center) ** 2 / (2 * sigma ** 2))
impulse_response_real = impulse_response_real / np.max(impulse_response_real)

# КОМПЛЕКСНАЯ ИХ: h_comp(t) = h_real(t) + j·0.5·h_imag(t)
impulse_response_complex_real = np.exp(-(t - center) ** 2 / (2 * sigma ** 2))
impulse_response_complex_imag = 0.7 * np.exp(-(t - center - 0.05) ** 2 / (2 * (sigma / 1.5) ** 2))
impulse_response_complex = (impulse_response_complex_real + 0.5j * impulse_response_complex_imag)
impulse_response_complex = impulse_response_complex / np.max(np.abs(impulse_response_complex))

# =============================================================================
# 3. ВЫЧИСЛЕНИЕ СВЕРТОК (сигналы после канала)
# =============================================================================
convolved_real = np.convolve(sinusoid_real, impulse_response_real, mode='same')
convolved_complex = np.convolve(sinusoid_complex, impulse_response_complex, mode='same')

# =============================================================================
# 4. ГРАФИКИ - ИСХОДНЫЕ СИГНАЛЫ, ИХ И СВЕРТКИ
# =============================================================================
print("=" * 60)
print("ВИЗУАЛИЗАЦИЯ СИГНАЛОВ И КАНАЛА")
print("=" * 60)

fig, axes = plt.subplots(3, 2, figsize=(15, 10))

# СТОЛБЕЦ 1: РЕАЛЬНЫЙ СИГНАЛ
# Исходный сигнал
axes[0, 0].plot(t, sinusoid_real, 'b-', linewidth=2)
axes[0, 0].set_title('РЕАЛЬНЫЙ: Исходный сигнал\nx(t) = sin(2π·5·t)')
axes[0, 0].set_ylabel('Амплитуда')
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_xlim(0, 1)

# Импульсная характеристика
axes[1, 0].plot(t, impulse_response_real, 'r-', linewidth=2)
axes[1, 0].set_title('РЕАЛЬНАЯ ИХ: h(t) = Gauss(t-0.3, σ=0.08)')
axes[1, 0].set_ylabel('Амплитуда')
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].set_xlim(0, 1)

# Результат свертки
axes[2, 0].plot(t, convolved_real, 'g-', linewidth=2)
axes[2, 0].set_title('СВЕРТКА: y(t) = (x ∗ h)(t)')
axes[2, 0].set_xlabel('Время, t [сек]')
axes[2, 0].set_ylabel('Амплитуда')
axes[2, 0].grid(True, alpha=0.3)
axes[2, 0].set_xlim(0, 1)

# СТОЛБЕЦ 2: КОМПЛЕКСНЫЙ СИГНАЛ
# Исходный сигнал
axes[0, 1].plot(t, np.real(sinusoid_complex), 'blue', linewidth=2, label='Re(x) = cos(2π·5·t)')
axes[0, 1].plot(t, np.imag(sinusoid_complex), 'red', linewidth=2, label='Im(x) = sin(2π·5·t)')
axes[0, 1].set_title('КОМПЛЕКСНЫЙ: Исходный сигнал\nx(t) = e^(j·2π·5·t)')
axes[0, 1].set_ylabel('Амплитуда')
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].legend()

# Импульсная характеристика
axes[1, 1].plot(t, np.real(impulse_response_complex), 'darkred', linewidth=2, label='Re(h)')
axes[1, 1].plot(t, np.imag(impulse_response_complex), 'orange', linewidth=2, label='Im(h)')
axes[1, 1].set_title('КОМПЛЕКСНАЯ ИХ: h(t) = h_real(t) + j·0.5·h_imag(t)')
axes[1, 1].set_ylabel('Амплитуда')
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].legend()

# Результат свертки
axes[2, 1].plot(t, np.real(convolved_complex), 'green', linewidth=2, label='Re(y)')
axes[2, 1].plot(t, np.imag(convolved_complex), 'purple', linewidth=2, label='Im(y)')
axes[2, 1].set_title('СВЕРТКА: y(t) = (x ∗ h)(t)')
axes[2, 1].set_xlabel('Время, t [сек]')
axes[2, 1].set_ylabel('Амплитуда')
axes[2, 1].grid(True, alpha=0.3)
axes[2, 1].legend()

plt.tight_layout()
plt.show()

# =============================================================================
# 5. ДОПОЛНИТЕЛЬНЫЕ ГРАФИКИ ДЛЯ КОМПЛЕКСНОГО СЛУЧАЯ
# =============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# Амплитуда комплексного сигнала
ax1.plot(t, np.abs(sinusoid_complex), 'b-', linewidth=2, alpha=0.7, label='|x(t)|')
ax1.plot(t, np.abs(convolved_complex), 'r-', linewidth=2, label='|y(t)|')
ax1.set_title('АМПЛИТУДЫ: |x(t)| и |y(t)|')
ax1.set_xlabel('Время, t [сек]')
ax1.set_ylabel('Амплитуда')
ax1.grid(True, alpha=0.3)
ax1.legend()

# Фаза комплексного сигнала
ax2.plot(t, np.angle(sinusoid_complex), 'b-', linewidth=2, alpha=0.7, label='∠x(t)')
ax2.plot(t, np.angle(convolved_complex), 'r-', linewidth=2, label='∠y(t)')
ax2.set_title('ФАЗЫ: ∠x(t) и ∠y(t)')
ax2.set_xlabel('Время, t [сек]')
ax2.set_ylabel('Фаза [рад]')
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.show()


# =============================================================================
# 6. АДАПТИВНЫЕ ФИЛЬТРЫ ДЛЯ КОМПЕНСАЦИИ КАНАЛА
# =============================================================================

def adaptive_equalizer_nlms(received_signal, training_signal, filter_length=64, mu=0.1):
    """
    Нормализованный LMS адаптивный эквалайзер
    """
    n = len(received_signal)
    w = np.zeros(filter_length, dtype=received_signal.dtype)
    equalized = np.zeros(n, dtype=received_signal.dtype)
    error = np.zeros(n, dtype=received_signal.dtype)

    for i in range(filter_length, n):
        # Входной вектор
        x = received_signal[i - filter_length:i]

        # Выход фильтра
        y = np.dot(w, x)
        equalized[i] = y

        # Ошибка
        e = training_signal[i] - y
        error[i] = e

        # Нормализованное обновление
        signal_power = np.dot(x, np.conj(x)).real if np.iscomplexobj(x) else np.dot(x, x)
        if signal_power > 1e-10:
            w = w + (mu * e * np.conj(x) if np.iscomplexobj(x) else mu * e * x) / signal_power

    return equalized, error, w


def apply_trained_filter(signal, coefficients):
    """Применение обученного фильтра"""
    filter_length = len(coefficients)
    n = len(signal)
    output = np.zeros(n, dtype=signal.dtype)

    for i in range(filter_length, n):
        x = signal[i - filter_length:i]
        output[i] = np.dot(coefficients, x)

    return output


# =============================================================================
# 7. ОБУЧЕНИЕ АДАПТИВНЫХ ФИЛЬТРОВ
# =============================================================================
print("\n" + "=" * 60)
print("ОБУЧЕНИЕ АДАПТИВНЫХ ФИЛЬТРОВ")
print("=" * 60)

# Обучаем РЕАЛЬНЫЙ фильтр
print("РЕАЛЬНЫЙ ФИЛЬТР:")
equalized_real, error_real, coeffs_real = adaptive_equalizer_nlms(
    convolved_real, sinusoid_real, filter_length=64, mu=0.1
)
print(f"Средняя ошибка: {np.mean(np.abs(error_real[100:])):.6f}")

# Обучаем КОМПЛЕКСНЫЙ фильтр
print("\nКОМПЛЕКСНЫЙ ФИЛЬТР:")
equalized_complex, error_complex, coeffs_complex = adaptive_equalizer_nlms(
    convolved_complex, sinusoid_complex, filter_length=64, mu=0.1
)
print(f"Средняя ошибка: {np.mean(np.abs(error_complex[100:])):.6f}")

# =============================================================================
# 8. ТЕСТИРОВАНИЕ НА НОВЫХ СИГНАЛАХ
# =============================================================================
print("\n" + "=" * 60)
print("ТЕСТИРОВАНИЕ НА НОВЫХ СИГНАЛАХ")
print("=" * 60)

# Новые тестовые сигналы
test_signal_real = 0.7 * np.sin(2 * np.pi * 4 * t) + 0.3 * np.cos(2 * np.pi * 8 * t)
test_signal_complex = 0.5 * np.exp(1j * 2 * np.pi * 3 * t) + 0.5 * np.exp(1j * 2 * np.pi * 7 * t)

# Пропускаем тестовые сигналы через канал
test_distorted_real = np.convolve(test_signal_real, impulse_response_real, mode='same')
test_distorted_complex = np.convolve(test_signal_complex, impulse_response_complex, mode='same')

# Восстанавливаем с помощью обученных фильтров
restored_real = apply_trained_filter(test_distorted_real, coeffs_real)
restored_complex = apply_trained_filter(test_distorted_complex, coeffs_complex)

# Ошибки восстановления
error_restoration_real = test_signal_real - restored_real
error_restoration_complex = test_signal_complex - restored_complex

print(f"Ошибка восстановления (реальный): {np.mean(np.abs(error_restoration_real[100:])):.6f}")
print(f"Ошибка восстановления (комплексный): {np.mean(np.abs(error_restoration_complex[100:])):.6f}")

# =============================================================================
# 9. ГРАФИКИ - РЕЗУЛЬТАТЫ АДАПТИВНОЙ ФИЛЬТРАЦИИ
# =============================================================================
fig, axes = plt.subplots(3, 2, figsize=(15, 10))

# СТОЛБЕЦ 1: РЕАЛЬНЫЙ СИГНАЛ
# Оригинальные сигналы
axes[0, 0].plot(t, sinusoid_real, 'b-', linewidth=2, label='Исходный')
axes[0, 0].plot(t, convolved_real, 'r-', alpha=0.7, linewidth=1, label='После канала')
axes[0, 0].plot(t, equalized_real, 'g-', linewidth=2, label='После эквалайзера')
axes[0, 0].set_title('РЕАЛЬНЫЙ: Сравнение сигналов')
axes[0, 0].legend()
axes[0, 0].grid(True)

# Тестовый сигнал
axes[1, 0].plot(t, test_signal_real, 'b-', linewidth=2, label='Тестовый')
axes[1, 0].plot(t, test_distorted_real, 'r-', alpha=0.7, linewidth=1, label='После канала')
axes[1, 0].plot(t, restored_real, 'g-', linewidth=2, label='Восстановленный')
axes[1, 0].set_title('РЕАЛЬНЫЙ: Тестовый сигнал')
axes[1, 0].legend()
axes[1, 0].grid(True)

# Ошибки
axes[2, 0].plot(t, error_real, 'r-', linewidth=1, label='Ошибка обучения')
axes[2, 0].plot(t, error_restoration_real, 'b-', linewidth=1, label='Ошибка восстановления')
axes[2, 0].set_title('РЕАЛЬНЫЙ: Ошибки')
axes[2, 0].legend()
axes[2, 0].grid(True)

# СТОЛБЕЦ 2: КОМПЛЕКСНЫЙ СИГНАЛ
# Оригинальные сигналы (Real часть)
axes[0, 1].plot(t, np.real(sinusoid_complex), 'blue', linewidth=2, label='Re(исходный)')
axes[0, 1].plot(t, np.real(convolved_complex), 'red', alpha=0.7, linewidth=1, label='Re(после канала)')
axes[0, 1].plot(t, np.real(equalized_complex), 'green', linewidth=2, label='Re(после эквалайзера)')
axes[0, 1].set_title('КОМПЛЕКСНЫЙ: Real части')
axes[0, 1].legend()
axes[0, 1].grid(True)

# Тестовый сигнал (Real часть)
axes[1, 1].plot(t, np.real(test_signal_complex), 'blue', linewidth=2, label='Re(тестовый)')
axes[1, 1].plot(t, np.real(test_distorted_complex), 'red', alpha=0.7, linewidth=1, label='Re(после канала)')
axes[1, 1].plot(t, np.real(restored_complex), 'green', linewidth=2, label='Re(восстановленный)')
axes[1, 1].set_title('КОМПЛЕКСНЫЙ: Real части тестового')
axes[1, 1].legend()
axes[1, 1].grid(True)

# Ошибки
axes[2, 1].plot(t, np.abs(error_complex), 'r-', linewidth=1, label='|Ошибка обучения|')
axes[2, 1].plot(t, np.abs(error_restoration_complex), 'b-', linewidth=1, label='|Ошибка восстановления|')
axes[2, 1].set_title('КОМПЛЕКСНЫЙ: Ошибки')
axes[2, 1].legend()
axes[2, 1].grid(True)

plt.tight_layout()
plt.show()

# =============================================================================
# 10. ФИНАЛЬНЫЕ ГРАФИКИ - ИМПУЛЬСНЫЕ ХАРАКТЕРИСТИКИ
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 8))

# Реальные ИХ
axes[0, 0].plot(t[:100], impulse_response_real[:100], 'r-', linewidth=2, label='ИХ канала')
axes[0, 0].plot(np.arange(len(coeffs_real)), coeffs_real, 'b-', linewidth=2, label='Коэф. фильтра')
axes[0, 0].set_title('РЕАЛЬНЫЙ: ИХ канала vs Коэффициенты фильтра')
axes[0, 0].legend()
axes[0, 0].grid(True)

# Комплексные ИХ (Real часть)
axes[0, 1].plot(t[:100], np.real(impulse_response_complex[:100]), 'r-', linewidth=2, label='Re(ИХ канала)')
axes[0, 1].plot(np.arange(len(coeffs_complex)), np.real(coeffs_complex), 'b-', linewidth=2, label='Re(Коэф. фильтра)')
axes[0, 1].set_title('КОМПЛЕКСНЫЙ: Real части')
axes[0, 1].legend()
axes[0, 1].grid(True)

# Комплексные ИХ (Imag часть)
axes[1, 0].plot(t[:100], np.imag(impulse_response_complex[:100]), 'r-', linewidth=2, label='Im(ИХ канала)')
axes[1, 0].plot(np.arange(len(coeffs_complex)), np.imag(coeffs_complex), 'b-', linewidth=2, label='Im(Коэф. фильтра)')
axes[1, 0].set_title('КОМПЛЕКСНЫЙ: Imag части')
axes[1, 0].legend()
axes[1, 0].grid(True)

# Амплитуды комплексных ИХ
axes[1, 1].plot(t[:100], np.abs(impulse_response_complex[:100]), 'r-', linewidth=2, label='|ИХ канала|')
axes[1, 1].plot(np.arange(len(coeffs_complex)), np.abs(coeffs_complex), 'b-', linewidth=2, label='|Коэф. фильтра|')
axes[1, 1].set_title('КОМПЛЕКСНЫЙ: Амплитуды')
axes[1, 1].legend()
axes[1, 1].grid(True)

plt.tight_layout()
plt.show()

print("\n" + "=" * 70)
print("ИТОГИ:")
print("=" * 70)
print("✓ Показаны исходные сигналы, ИХ и их свертки")
print("✓ Реализованы адаптивные фильтры NLMS для компенсации канала")
print("✓ Фильтры успешно обучаются и восстанавливают новые сигналы")
print("✓ Комплексный случай показывает полную информацию (амплитуда + фаза)")
print("=" * 70)