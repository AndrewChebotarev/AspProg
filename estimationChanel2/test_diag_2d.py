"""
Diagnostika: L(T1,T2) dlya N=2 - sovmestnaya ocenka dvuh zaderzhek
"""
import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)
print(f"Z_prime = {Z_prime}")

# 2D scan of L(T1,T2)
N_grid = 101
tau1_range = np.linspace(0, 2000e-6, N_grid)
tau2_range = np.linspace(0, 2000e-6, N_grid)

L_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:  # symmetry: tau2 > tau1
            L_grid[i, j] = -np.inf
            continue
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L_grid[i, j] = compute_likelihood(Z_prime, H, Q)

# Find maximum
idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
t1_est = tau1_range[idx_max[0]]
t2_est = tau2_range[idx_max[1]]
print(f"\n2D max L(T1,T2) at: tau1={t1_est*1e6:.1f} us, tau2={t2_est*1e6:.1f} us")
print(f"True: tau1={tau_true[0]*1e6:.0f} us, tau2={tau_true[1]*1e6:.0f} us")
print(f"L_max = {L_grid[idx_max]:.6f}")

# L at true values
i_true = np.argmin(np.abs(tau1_range - tau_true[0]))
j_true = np.argmin(np.abs(tau2_range - tau_true[1]))
print(f"L(true) = {L_grid[i_true, j_true]:.6f}")

# Plot
plt.figure(figsize=(10, 8))
plt.imshow(L_grid, extent=[0, 2000, 0, 2000], origin='lower', aspect='auto', cmap='viridis')
plt.colorbar(label='L(T1,T2)')
plt.plot(tau_true[0]*1e6, tau_true[1]*1e6, 'r*', markersize=15, label='True')
plt.plot(t1_est*1e6, t2_est*1e6, 'wx', markersize=10, label='Estimated')
plt.xlabel('tau1, us')
plt.ylabel('tau2, us')
plt.title('L(T1,T2) for N=2')
plt.legend()
plt.tight_layout()
plt.show()

# Also check: what if we fix tau1 and scan tau2?
print("\n=== Fix tau1=600 us, scan tau2 ===")
tau1_fixed = 600e-6
tau2_scan = np.linspace(0, 2000e-6, 2001)
L_tau2 = []
for t2 in tau2_scan:
    if t2 <= tau1_fixed:
        L_tau2.append(-np.inf)
        continue
    tau = np.array([tau1_fixed, t2])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    L_tau2.append(compute_likelihood(Z_prime, H, Q))
L_tau2 = np.array(L_tau2)
idx = np.argmax(L_tau2)
print(f"Max at tau2={tau2_scan[idx]*1e6:.1f} us (true={tau_true[1]*1e6:.0f} us)")
print(f"L(tau2=1200)={L_tau2[np.argmin(np.abs(tau2_scan-1200e-6))]:.6f}")

# Fix tau2=1200 us, scan tau1
print("\n=== Fix tau2=1200 us, scan tau1 ===")
tau2_fixed = 1200e-6
tau1_scan = np.linspace(0, 2000e-6, 2001)
L_tau1 = []
for t1 in tau1_scan:
    if t1 >= tau2_fixed:
        L_tau1.append(-np.inf)
        continue
    tau = np.array([t1, tau2_fixed])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    L_tau1.append(compute_likelihood(Z_prime, H, Q))
L_tau1 = np.array(L_tau1)
idx = np.argmax(L_tau1)
print(f"Max at tau1={tau1_scan[idx]*1e6:.1f} us (true={tau_true[0]*1e6:.0f} us)")
print(f"L(tau1=600)={L_tau1[np.argmin(np.abs(tau1_scan-600e-6))]:.6f}")