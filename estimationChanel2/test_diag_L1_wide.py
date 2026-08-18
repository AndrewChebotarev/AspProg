"""
L(T) dlya N=1 na signale s 2 luchami v shirokom diapazone
"""
import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; dt=1/fs; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,_,_ = compute_sufficient_statistics(xi, t, omega_m, a_m, N0, fs)

# L(T) dlya N=1 v shirokom diapazone
tau_range = np.linspace(0, T_dur, 2001)
L_vals = []
for tau_test in tau_range:
    tau = np.array([tau_test])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    
    A_est = np.sqrt(A_cs[0]**2 + A_cs[1]**2)
    psi_est = np.arctan2(A_cs[1], A_cs[0])
    
    s = np.zeros_like(t)
    for m in range(M):
        s += A_est * a_m[m] * np.cos(omega_m[m] * (t - tau_test) - phi_m[m] - psi_est)
    
    int_xi_s = np.sum(xi * s) * dt
    int_s2 = np.sum(s**2) * dt
    L = (2.0/N0) * int_xi_s - (1.0/N0) * int_s2
    L_vals.append(L)

L_vals = np.array(L_vals)
idx_max = np.argmax(L_vals)
print(f"Max L(T) pri tau={tau_range[idx_max]*1e6:.1f} mks")
print(f"L_max = {L_vals[idx_max]:.10f}")
print(f"Istinnye tau: tau1={tau_true[0]*1e6:.0f}, tau2={tau_true[1]*1e6:.0f} mks")

# Grafik
plt.figure(figsize=(12, 6))
plt.plot(tau_range*1e6, L_vals, 'b-', linewidth=1)
plt.axvline(tau_true[0]*1e6, color='g', linestyle='--', label=f'tau1={tau_true[0]*1e6:.0f} mks')
plt.axvline(tau_true[1]*1e6, color='orange', linestyle='--', label=f'tau2={tau_true[1]*1e6:.0f} mks')
plt.axvline(tau_range[idx_max]*1e6, color='r', linestyle=':', label=f'max={tau_range[idx_max]*1e6:.1f} mks')
plt.xlabel('tau, mks')
plt.ylabel('L(T)')
plt.title('L(T) dlya N=1 na signale s 2 luchami')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# To zhe dlya signala s 1 luchom
data1=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,1,np.array([A_true[0]]),np.array([tau_true[0]]),np.array([psi_true[0]]),fs,T1,T2,40)
xi1=data1['xi_observed'].copy()
Z_prime1,_,_ = compute_sufficient_statistics(xi1, t, omega_m, a_m, N0, fs)

L_vals1 = []
for tau_test in tau_range:
    tau = np.array([tau_test])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime1 @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    
    A_est = np.sqrt(A_cs[0]**2 + A_cs[1]**2)
    psi_est = np.arctan2(A_cs[1], A_cs[0])
    
    s = np.zeros_like(t)
    for m in range(M):
        s += A_est * a_m[m] * np.cos(omega_m[m] * (t - tau_test) - phi_m[m] - psi_est)
    
    int_xi_s = np.sum(xi1 * s) * dt
    int_s2 = np.sum(s**2) * dt
    L = (2.0/N0) * int_xi_s - (1.0/N0) * int_s2
    L_vals1.append(L)

L_vals1 = np.array(L_vals1)
idx_max1 = np.argmax(L_vals1)
print(f"\nSignal s 1 luchom: max L(T) pri tau={tau_range[idx_max1]*1e6:.1f} mks")

plt.figure(figsize=(12, 6))
plt.plot(tau_range*1e6, L_vals1, 'r-', linewidth=1)
plt.axvline(tau_true[0]*1e6, color='g', linestyle='--', label=f'tau1={tau_true[0]*1e6:.0f} mks')
plt.axvline(tau_range[idx_max1]*1e6, color='r', linestyle=':', label=f'max={tau_range[idx_max1]*1e6:.1f} mks')
plt.xlabel('tau, mks')
plt.ylabel('L(T)')
plt.title('L(T) dlya N=1 na signale s 1 luchom')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()