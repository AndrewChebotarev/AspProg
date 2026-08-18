"""
Proverka: formula (98) + (99) dlya N=2
L(tau1,tau2) = Z' @ H @ A_cs_hat^T - 0.5 * A_cs_hat @ Q @ A_cs_hat^T
gde A_cs_hat = Z' @ H @ Q^{-1} (formula 99)
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

def compute_likelihood_98(Z_prime, H, Q):
    """Formula (98): L = Z' @ H @ A_cs^T - 0.5 * A_cs @ Q @ A_cs^T
       gde A_cs = Z' @ H @ Q^{-1} (formula 99)"""
    try:
        Q_inv = np.linalg.inv(Q)
        A_cs = Z_prime @ H @ Q_inv  # formula (99)
        M_val = Z_prime @ H @ A_cs  # = A_cs @ Q @ A_cs^T (esli Z' @ H = A_cs @ Q)
        Q_val = A_cs @ Q @ A_cs
        return M_val - 0.5 * Q_val
    except:
        return -np.inf

print("=== Formula (98)+(99): L = Z'@H@A_cs^T - 0.5*A_cs@Q@A_cs^T ===")
print()

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L = compute_likelihood_98(Z_prime, H, Q)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.10f}")

# Poisk maksimuma
print("\n=== Poisk maksimuma L(T1,T2) po formule (98)+(99) ===")
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
        L_grid[i,j] = compute_likelihood_98(Z_prime, H, Q)

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
print(f"Max at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.10f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"L(true) = {L_grid[i_true, j_true]:.10f}")

# Sravnenie s formulo (100)
print("\n=== Sravnenie s formulo (100) ===")
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L100 = compute_likelihood(Z_prime, H, Q)
        L98 = compute_likelihood_98(Z_prime, H, Q)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L(100)={L100:.10f}, L(98)={L98:.10f}")