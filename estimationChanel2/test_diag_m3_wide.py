"""
Proverka: M=3, shirokiy raznos zaderzhek
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

# M=3, N=2, shirokiy raznos
M=3; a_m=np.array([1.0,1.0,1.0]); omega_m=2*np.pi*np.array([700, 900, 1100]); phi_m=np.array([0,0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([100e-6, 1900e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)
print(f"Z_prime (M=3) = {Z_prime}")

# 2D scan
print("\n=== L(T) dlya M=3, shirokiy raznos ===")
for t1 in [50e-6, 100e-6, 150e-6]:
    for t2 in [1850e-6, 1900e-6, 1950e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L = compute_likelihood(Z_prime, H, Q)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.10f}")

# Poisk maksimuma
print("\n=== Poisk maksimuma ===")
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
print(f"Max at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.10f}")

i_true = np.argmin(np.abs(tau1_range - 100e-6))
j_true = np.argmin(np.abs(tau2_range - 1900e-6))
print(f"L(true) = {L_grid[i_true, j_true]:.10f}")

# Sravnenie s M=2
print("\n=== Sravnenie s M=2, shirokiy raznos ===")
M2=2; a_m2=np.array([1.0,1.0]); omega_m2=2*np.pi*np.array([900,1100]); phi_m2=np.array([0,0])
data2=generate_test_signal_and_channel_response(M2,a_m2,omega_m2,phi_m2,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi2=data2['xi_observed'].copy()
Z2,_,_=compute_sufficient_statistics(xi2,t,omega_m2,a_m2,N0,fs)

for t1 in [50e-6, 100e-6, 150e-6]:
    for t2 in [1850e-6, 1900e-6, 1950e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m2, phi_m2, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m2, a_m2, tau, N0)
        L = compute_likelihood(Z2, H, Q)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.10f}")