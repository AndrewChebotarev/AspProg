"""
Proverka: H^T @ H vs Q dlya N=1 i N=2
"""
import numpy as np
from MonteCarlo import compute_matrices_C_S, compute_matrix_H, compute_matrix_Q

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N0=100.0

print("=== N=1: H^T @ H vs Q ===")
tau1 = np.array([600e-6])
C1,S1 = compute_matrices_C_S(omega_m, phi_m, tau1)
H1 = compute_matrix_H(C1, S1)
Q1 = compute_matrix_Q(omega_m, a_m, tau1, N0)

HtH1 = H1.T @ H1
print(f"H^T @ H =\n{HtH1}")
print(f"Q =\n{Q1}")
print(f"Raznica = {np.max(np.abs(HtH1 - Q1)):.15f}")

# Dlya N=1, H^T @ H = Q * (N0/2) * (2/sum(a_m^2))? 
# Q = diag(sum(2*a_m^2/N0), sum(2*a_m^2/N0))
# H^T @ H = ?
print(f"\nH1.T @ H1 diagonal: {np.diag(HtH1)}")
print(f"Q1 diagonal: {np.diag(Q1)}")
print(f"sum(a_m^2) = {np.sum(a_m**2)}")
print(f"2*sum(a_m^2)/N0 = {2*np.sum(a_m**2)/N0}")

print("\n=== N=2: H^T @ H vs Q ===")
tau2 = np.array([600e-6, 1200e-6])
C2,S2 = compute_matrices_C_S(omega_m, phi_m, tau2)
H2 = compute_matrix_H(C2, S2)
Q2 = compute_matrix_Q(omega_m, a_m, tau2, N0)

HtH2 = H2.T @ H2
print(f"H^T @ H =\n{HtH2}")
print(f"Q =\n{Q2}")
print(f"Raznica = {np.max(np.abs(HtH2 - Q2)):.15f}")

# Sravnenie elementov
print(f"\nH^T @ H diagonal: {np.diag(HtH2)}")
print(f"Q diagonal: {np.diag(Q2)}")

# Proverim: H^T @ H dolzhno byt ravno Q * (N0/2) * (T_dur/2)?
# Dlya N=1: H^T @ H = diag(2, 2), Q = diag(0.04, 0.04)
# 2 / 0.04 = 50 = N0/2 * ? 
print(f"\nH^T @ H / Q = {HtH2 / Q2}")

# Esche odna proverka: formula dlya Q
# Qc_ik = sum_m(2*a_m^2/N0 * cos(omega_m*(tau_i - tau_k)))
# (H^T @ H)_ik = ?
# H = [C, -S; S, C]
# H^T @ H = [C^T@C + S^T@S, C^T@(-S) + S^T@C; (-S)^T@C + C^T@S, (-S)^T@(-S) + C^T@C]
# (C^T@C + S^T@S)_ik = sum_m(cos(phase_mi)*cos(phase_mk) + sin(phase_mi)*sin(phase_mk))
# = sum_m(cos(phase_mi - phase_mk)) = sum_m(cos(omega_m*(tau_i - tau_k)))
# No v Q est mnozhitel 2*a_m^2/N0, a v H^T@H net!
print(f"\nBez a_m: H^T @ H (bez a_m) =")
# Esli a_m = [1,1], to H soderzhit tolko cos i sin
print(HtH2)
print(f"\nS a_m^2:")
# (H * diag(a_m))^T @ (H * diag(a_m))
H2_weighted = H2.copy()
H2_weighted[:M] *= a_m[:, np.newaxis]  # a_m dlya C
H2_weighted[M:] *= a_m[:, np.newaxis]  # a_m dlya S
HtH_weighted = H2_weighted.T @ H2_weighted
print(HtH_weighted)
print(f"\nQ * N0/2 =")
print(Q2 * N0 / 2)
print(f"Raznica = {np.max(np.abs(HtH_weighted - Q2 * N0 / 2)):.15f}")