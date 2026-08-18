"""
Proverka skorosti: pryamaya FOP dlya N=2
"""
import numpy as np
import time
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import compute_sufficient_statistics, compute_matrices_C_S, compute_matrix_H

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
t=data['t']
xi=data['xi_observed'].copy()
dt=1/fs

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

def compute_signal_from_A_cs(t, A_cs, omega_m, a_m, phi_m, tau):
    N = len(tau)
    A_c = A_cs[:N]
    A_s = A_cs[N:]
    A = np.sqrt(A_c**2 + A_s**2)
    psi = np.arctan2(A_s, A_c)
    
    s = np.zeros_like(t)
    for k in range(N):
        for m in range(M):
            s += A[k] * a_m[m] * np.cos(omega_m[m] * (t - tau[k]) - phi_m[m] - psi[k])
    return s

def estimate_one_point(tau, xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur):
    """Ocenka dlya odnogo kandidata (tau1, tau2)"""
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_Q_fixed(omega_m, a_m, tau, N0)
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    
    s = compute_signal_from_A_cs(t, A_cs, omega_m, a_m, phi_m, tau)
    int_xi_s = np.sum(xi * s) * dt
    int_s2 = np.sum(s**2) * dt
    L = (2.0/N0) * int_xi_s - (1.0/N0) * int_s2
    
    N = len(tau)
    A_est = np.sqrt(A_cs[:N]**2 + A_cs[N:]**2)
    psi_est = np.arctan2(A_cs[N:], A_cs[:N])
    
    return L, A_est, psi_est

# Proverka skorosti
Z_prime,_,_ = compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

print("Zamer skorosti...")
n_tests = 50
start = time.time()
for _ in range(n_tests):
    tau = np.array([600e-6, 1200e-6])
    L, A, psi = estimate_one_point(tau, xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur)
end = time.time()
print(f"Srednee vremya na odnu tochku: {(end-start)/n_tests*1000:.2f} ms")
print(f"Ocenok v minutu: {60/((end-start)/n_tests):.0f}")

# Dlya setki 51x51=2601 tochek:
n_grid = 51
total_points = n_grid * (n_grid - 1) // 2
print(f"\nDlya setki {n_grid}x{n_grid} = {total_points} tochek:")
print(f"  Vremya: {total_points * (end-start)/n_tests:.1f} s")

# Proverim rezultat dlya istinnyh tau
L, A, psi = estimate_one_point(np.array([600e-6, 1200e-6]), xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur)
print(f"\nIstinnye tau: L={L:.10f}, A={A}, psi={psi}")

# Dlya sosednih
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        if t2 <= t1:
            continue
        L, A, psi = estimate_one_point(np.array([t1, t2]), xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.10f}, A={A}, psi={psi}")