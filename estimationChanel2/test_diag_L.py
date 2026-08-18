"""
Диагностика: сравниваем L(T) для N=1 (работает) и N=1 из состава N=2 (не работает)
"""
import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

# Генерируем сигнал с двумя лучами
data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

# Достаточные статистики
Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)
print(f"Z_prime = {Z_prime}")

# Сканируем L(T) для одного луча в широком диапазоне
tau_range = np.linspace(0, 2000e-6, 2001)
L_vals = []
for tau_test in tau_range:
    tau = np.array([tau_test])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    L_vals.append(compute_likelihood(Z_prime, H, Q))
L_vals = np.array(L_vals)

idx_max = np.argmax(L_vals)
print(f"\nМаксимум L(T) для N=1 (из состава N=2): tau={tau_range[idx_max]*1e6:.1f} мкс")
print(f"  L_max = {L_vals[idx_max]:.4f}")
print(f"  Истинная задержка первого луча: {tau_true[0]*1e6:.0f} мкс")

# Для сравнения: генерируем сигнал ТОЛЬКО с одним лучом
data1 = generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,1,np.array([A_true[0]]),np.array([tau_true[0]]),np.array([psi_true[0]]),fs,T1,T2,40)
xi1 = data1['xi_observed'].copy()

Z_prime1,_,_ = compute_sufficient_statistics(xi1,t,omega_m,a_m,N0,fs)
print(f"\nZ_prime (1 луч) = {Z_prime1}")

L_vals1 = []
for tau_test in tau_range:
    tau = np.array([tau_test])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    L_vals1.append(compute_likelihood(Z_prime1, H, Q))
L_vals1 = np.array(L_vals1)

idx_max1 = np.argmax(L_vals1)
print(f"\nМаксимум L(T) для N=1 (только 1 луч): tau={tau_range[idx_max1]*1e6:.1f} мкс")
print(f"  L_max = {L_vals1[idx_max1]:.4f}")

# График
plt.figure(figsize=(12, 6))
plt.plot(tau_range*1e6, L_vals, 'b-', label='N=1 из состава N=2 (сигнал с 2 лучами)')
plt.plot(tau_range*1e6, L_vals1, 'r-', label='N=1 (только 1 луч)')
plt.axvline(tau_true[0]*1e6, color='g', linestyle='--', label=f'Истинная τ₁={tau_true[0]*1e6:.0f} мкс')
plt.axvline(tau_true[1]*1e6, color='orange', linestyle='--', label=f'Истинная τ₂={tau_true[1]*1e6:.0f} мкс')
plt.xlabel('τ, мкс')
plt.ylabel('L(T)')
plt.title('Сравнение L(T) для N=1: сигнал с 1 лучом vs сигнал с 2 лучами')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Дополнительно: посмотрим на спектр L(T) вблизи
print("\n--- Детальный анализ L(T) для N=1 из состава N=2 ---")
tau_fine = np.linspace(400e-6, 800e-6, 401)
L_fine = []
for tau_test in tau_fine:
    tau = np.array([tau_test])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    L_fine.append(compute_likelihood(Z_prime, H, Q))
L_fine = np.array(L_fine)
idx_fine = np.argmax(L_fine)
print(f"Максимум на [400, 800] мкс: tau={tau_fine[idx_fine]*1e6:.1f} мкс, L={L_fine[idx_fine]:.4f}")
print(f"L при tau=600 мкс: {L_fine[np.argmin(np.abs(tau_fine-600e-6))]:.4f}")