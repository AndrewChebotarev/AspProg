"""
Proverka: H @ Q^{-1} @ H^T dlya ispravlennoy Q
"""
import numpy as np
from MonteCarlo import compute_matrices_C_S, compute_matrix_H

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N0=100.0

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
                sum_sin += weight * np.sin(omega_m[m] * delta)
            Qc[i, k] = sum_cos
            Qcs[i, k] = sum_sin
    Q = np.vstack([np.hstack([Qc, Qcs]), np.hstack([Qcs.T, Qc])])
    return Q

print("=== Proverka: H @ Q^{-1} @ H^T ===")
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_Q_fixed(omega_m, a_m, tau, N0)
        Q_inv = np.linalg.inv(Q)
        HQH = H @ Q_inv @ H.T
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): diag(HQH) = {np.diag(HQH)}")

# Sravnenie s originalnoy Q
print("\n=== S originalnoy Q ===")
from MonteCarlo2 import compute_matrix_Q as compute_Q_orig

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_Q_orig(omega_m, a_m, tau, N0)
        Q_inv = np.linalg.inv(Q)
        HQH = H @ Q_inv @ H.T
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): diag(HQH) = {np.diag(HQH)}")

# Proverim: chto daet formula (100) s ispravlennoy Q?
# L = 0.5 * Z' @ H @ Q^{-1} @ H^T @ Z'^T
# Esli H @ Q^{-1} @ H^T = const, to L = const dlya vseh tau
print("\n=== H @ Q^{-1} @ H^T (fixed Q) ===")
tau = np.array([600e-6, 1200e-6])
C, S = compute_matrices_C_S(omega_m, phi_m, tau)
H = compute_matrix_H(C, S)
Q = compute_Q_fixed(omega_m, a_m, tau, N0)
Q_inv = np.linalg.inv(Q)
HQH = H @ Q_inv @ H.T
print(HQH)

print("\n=== H @ Q^{-1} @ H^T (original Q) ===")
from MonteCarlo2 import compute_matrix_Q as compute_Q_orig
Q_orig = compute_Q_orig(omega_m, a_m, tau, N0)
Q_orig_inv = np.linalg.inv(Q_orig)
HQH_orig = H @ Q_orig_inv @ H.T
print(HQH_orig)