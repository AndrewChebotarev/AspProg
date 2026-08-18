"""
Диагностика: проверяем сигнал и достаточные статистики
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

# Генерируем сигнал с двумя лучами
data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']
s_useful = data['s_useful']

print(f"Сигнал с 2 лучами, SNR=40 дБ:")
print(f"  xi: min={xi.min():.6f}, max={xi.max():.6f}, mean={xi.mean():.6f}, std={xi.std():.6f}")
print(f"  s_useful: min={s_useful.min():.6f}, max={s_useful.max():.6f}, std={s_useful.std():.6f}")

# Вычисляем мощность сигнала и шума
noise = xi - s_useful
signal_power = np.mean(s_useful**2)
noise_power = np.mean(noise**2)
print(f"\nМощность сигнала: {signal_power:.10f}")
print(f"Мощность шума: {noise_power:.10f}")
print(f"SNR фактический: {10*np.log10(signal_power/noise_power):.2f} дБ")

# Достаточные статистики
Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)
print(f"\nДостаточные статистики (N0={N0}):")
print(f"  x_m = {x_m}")
print(f"  y_m = {y_m}")
print(f"  Z_prime = {Z_prime}")

# Проверяем: что если уменьшить N0?
for N0_test in [100.0, 10.0, 1.0, 0.1]:
    Zp,_,_ = compute_sufficient_statistics(xi,t,omega_m,a_m,N0_test,fs)
    print(f"\nN0={N0_test}: Z_prime = {Zp}")

# Для сравнения: сигнал с 1 лучом
data1=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,1,np.array([A_true[0]]),np.array([tau_true[0]]),np.array([psi_true[0]]),fs,T1,T2,40)
xi1=data1['xi_observed'].copy()
s1 = data1['s_useful']
print(f"\nСигнал с 1 лучом, SNR=40 дБ:")
print(f"  xi: min={xi1.min():.6f}, max={xi1.max():.6f}, std={xi1.std():.6f}")
print(f"  s_useful: min={s1.min():.6f}, max={s1.max():.6f}, std={s1.std():.6f}")

noise1 = xi1 - s1
sp1 = np.mean(s1**2)
np1 = np.mean(noise1**2)
print(f"  Мощность сигнала: {sp1:.10f}")
print(f"  Мощность шума: {np1:.10f}")
print(f"  SNR: {10*np.log10(sp1/np1):.2f} дБ")

for N0_test in [100.0, 10.0, 1.0, 0.1]:
    Zp,_,_ = compute_sufficient_statistics(xi1,t,omega_m,a_m,N0_test,fs)
    print(f"  N0={N0_test}: Z_prime = {Zp}")

# Ключевой тест: L(T) для N=1 из состава N=2 при разных N0
print("\n\n=== L(T) для N=1 из состава N=2 при разных N0 ===")
for N0_test in [100.0, 10.0, 1.0, 0.1]:
    Zp,_,_ = compute_sufficient_statistics(xi,t,omega_m,a_m,N0_test,fs)
    tau_range = np.linspace(0, 2000e-6, 2001)
    L_vals = []
    for tau_test in tau_range:
        tau = np.array([tau_test])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0_test)
        L_vals.append(compute_likelihood(Zp, H, Q))
    L_vals = np.array(L_vals)
    idx_max = np.argmax(L_vals)
    print(f"N0={N0_test}: max L(T) при tau={tau_range[idx_max]*1e6:.1f} мкс, L_max={L_vals[idx_max]:.6f}")
    print(f"  L(600 мкс)={L_vals[np.argmin(np.abs(tau_range-600e-6))]:.6f}")
    print(f"  L(1200 мкс)={L_vals[np.argmin(np.abs(tau_range-1200e-6))]:.6f}")