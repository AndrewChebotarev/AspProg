"""
Proverka: rabotaet li L(T) dlya N=2 s M=3 chastotami?
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

# M=3, N=2
M=3; a_m=np.array([1.0,1.0,1.0]); omega_m=2*np.pi*np.array([700, 900, 1100]); phi_m=np.array([0,0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)
print(f"Z_prime (M=3) = {Z_prime}")
print(f"Razmer Z_prime = {len(Z_prime)} (2M={2*M})")

# 2D scan
N_grid = 51
tau1_range = np.linspace(0, 2000e-6, N_grid)
tau2_range = np.linspace(0, 2000e-6, N_grid)
L_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            L_grid[i,j] = -np.inf
            continue
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L_grid[i,j] = compute_likelihood(Z_prime, H, Q)

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
print(f"\nMax at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.6f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"L(true) = {L_grid[i_true, j_true]:.6f}")

# Sravnenie dlya neskolkih tochek
print("\n=== L(T) dlya raznyh (tau1, tau2) ===")
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L = compute_likelihood(Z_prime, H, Q)
        print(f"  tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.6f}")