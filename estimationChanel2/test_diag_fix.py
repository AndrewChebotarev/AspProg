"""
Ispravlenie: Qcs_ik = sin(omega_m * (tau_i - tau_k)) a ne sin(omega_m * (tau_k - tau_i))
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
t=data['t']
xi=data['xi_observed'].copy()

# Ispravlennaya Q
def compute_Q_fixed(omega_m, a_m, tau, N0):
    M = len(omega_m)
    N = len(tau)
    Qc = np.zeros((N, N))
    Qcs = np.zeros((N, N))
    for i in range(N):
        for k in range(N):
            sum_cos = 0.0
            sum_sin = 0.0
            for m in range(M):
                weight = 2 * a_m[m] ** 2 / N0
                delta = tau[i] - tau[k]
                sum_cos += weight * np.cos(omega_m[m] * delta)
                sum_sin += weight * np.sin(omega_m[m] * delta)  # tau[i] - tau[k], a ne tau[k] - tau[i]
            Qc[i, k] = sum_cos
            Qcs[i, k] = sum_sin
    Q = np.vstack([np.hstack([Qc, Qcs]), np.hstack([Qcs.T, Qc])])
    return Q

# Proverim: H^T @ H == (N0/2) * Q_fixed ?
from MonteCarlo import compute_matrices_C_S, compute_matrix_H

tau = np.array([600e-6, 1200e-6])
C, S = compute_matrices_C_S(omega_m, phi_m, tau)
H = compute_matrix_H(C, S)
Q_fixed = compute_Q_fixed(omega_m, a_m, tau, N0)

HtH = H.T @ H
print("H^T @ H =")
print(HtH)
print("\nQ_fixed * N0/2 =")
print(Q_fixed * N0 / 2)
print(f"\nRaznica = {np.max(np.abs(HtH - Q_fixed * N0 / 2)):.15f}")

# Sovpadaet!
print("\n=== S ISPRAVLENNOY Q ===")

from MonteCarlo2 import compute_sufficient_statistics, compute_likelihood, compute_matrices_C_S, compute_matrix_H

Z_prime,x_m,y_m = compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

def compute_likelihood_fixed(Z_prime, H, Q):
    try:
        Q_inv = np.linalg.inv(Q)
        return 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    except:
        return -np.inf

print("\nL(T) s ispravlennoy Q:")
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_Q_fixed(omega_m, a_m, tau, N0)
        L = compute_likelihood_fixed(Z_prime, H, Q)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.10f}")

# Poisk maksimuma
print("\n=== Poisk maksimuma s ispravlennoy Q ===")
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
        Q = compute_Q_fixed(omega_m, a_m, tau, N0)
        L_grid[i,j] = compute_likelihood_fixed(Z_prime, H, Q)

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
print(f"Max at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.10f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"L(true) = {L_grid[i_true, j_true]:.10f}")