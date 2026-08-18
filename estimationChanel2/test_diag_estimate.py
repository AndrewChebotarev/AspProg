"""
Proverka: rabotaet li ocenka A_cs cherez (99) pri istinnyh tau?
I est' li maksimum L po A_cs pri istinnyh A_cs?
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
t=data['t']
xi=data['xi_observed']
dt=1/fs

# Schitaem Z_prime
Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

# Istinnye A_cs
A_cs_true = np.array([A_true[0]*np.cos(psi_true[0]), A_true[1]*np.cos(psi_true[1]), 
                       A_true[0]*np.sin(psi_true[0]), A_true[1]*np.sin(psi_true[1])])
print(f"Istinnye A_cs = {A_cs_true}")

# Ocenka A_cs cherez (99) pri istinnyh tau
tau = np.array([600e-6, 1200e-6])
C, S = compute_matrices_C_S(omega_m, phi_m, tau)
H = compute_matrix_H(C, S)
Q = compute_matrix_Q(omega_m, a_m, tau, N0)

Q_inv = np.linalg.inv(Q)
A_cs_est = Z_prime @ H @ Q_inv
A_cs_est = A_cs_est * (2.0 / T_dur)  # normirovka
print(f"Ocenka A_cs (99) = {A_cs_est}")

# Vosstanavlivaem A i psi
N = 2
A_c = A_cs_est[:N]
A_s = A_cs_est[N:]
A_est = np.sqrt(A_c**2 + A_s**2)
psi_est = np.arctan2(A_s, A_c)
print(f"A_est = {A_est}, A_true = {A_true}")
print(f"psi_est = {psi_est}, psi_true = {psi_true}")

# Teper: ispolzuem pryamuyu FOP dlya ocenki
def compute_signal(t, A, tau, psi):
    s = np.zeros_like(t)
    for k in range(len(A)):
        for m in range(M):
            s += A[k] * a_m[m] * np.cos(omega_m[m] * (t - tau[k]) - phi_m[m] - psi[k])
    return s

def likelihood_direct(xi, t, A, tau, psi, N0, dt):
    s = compute_signal(t, A, tau, psi)
    int_xi_s = np.sum(xi * s) * dt
    int_s2 = np.sum(s**2) * dt
    return (2.0/N0) * int_xi_s - (1.0/N0) * int_s2

# Sravnenie: L pryamaya vs L cherez (100) dlya raznyh tau
print("\n=== Sravnenie: pryamaya FOP vs formula (100) ===")
print(f"{'tau1':>6} {'tau2':>6} | {'L_direct':>12} {'L_100':>12} {'L_98':>12}")
print("-" * 55)

for t1 in [400e-6, 500e-6, 600e-6, 700e-6, 800e-6]:
    for t2 in [1000e-6, 1100e-6, 1200e-6, 1300e-6, 1400e-6]:
        if t2 <= t1:
            continue
        # Pryamaya FOP s istinnymi A, psi
        L_dir = likelihood_direct(xi, t, A_true, np.array([t1, t2]), psi_true, N0, dt)
        
        # Formula (100)
        tau_arr = np.array([t1, t2])
        C2, S2 = compute_matrices_C_S(omega_m, phi_m, tau_arr)
        H2 = compute_matrix_H(C2, S2)
        Q2 = compute_matrix_Q(omega_m, a_m, tau_arr, N0)
        L_100 = compute_likelihood(Z_prime, H2, Q2)
        
        # Formula (98) s ocenkoj A_cs
        Q2_inv = np.linalg.inv(Q2)
        A_cs_hat = Z_prime @ H2 @ Q2_inv
        M_val = Z_prime @ H2 @ A_cs_hat
        Q_val = A_cs_hat @ Q2 @ A_cs_hat
        L_98 = M_val - 0.5 * Q_val
        
        print(f"{t1*1e6:6.0f} {t2*1e6:6.0f} | {L_dir:12.8f} {L_100:12.8f} {L_98:12.8f}")