"""
Proverka: ocenka A_cs cherez psevdoobrashenie H vmesto Q^{-1}
Formula: Z' = (T_dur/2) * H * A_cs^T  => A_cs = (2/T_dur) * pinv(H) * Z'
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

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

A_cs_true = np.array([A_true[0]*np.cos(psi_true[0]), A_true[1]*np.cos(psi_true[1]), 
                       A_true[0]*np.sin(psi_true[0]), A_true[1]*np.sin(psi_true[1])])

def estimate_A_cs_pinv(Z_prime, omega_m, a_m, phi_m, tau, T_dur):
    """Ocenka A_cs cherez pinv(H)"""
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    H_pinv = np.linalg.pinv(H)
    A_cs = (2.0 / T_dur) * Z_prime @ H_pinv.T
    return A_cs

def estimate_A_cs_Q(Z_prime, omega_m, a_m, phi_m, tau, N0, T_dur):
    """Ocenka A_cs cherez Q^{-1} (formula 99)"""
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    return A_cs

print("=== Sravnenie ocenok A_cs ===")
print(f"Istina: A_cs = {A_cs_true}")
print()

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        A_pinv = estimate_A_cs_pinv(Z_prime, omega_m, a_m, phi_m, tau, T_dur)
        A_Q = estimate_A_cs_Q(Z_prime, omega_m, a_m, phi_m, tau, N0, T_dur)
        
        N = 2
        A_pinv_est = np.sqrt(A_pinv[:N]**2 + A_pinv[N:]**2)
        A_Q_est = np.sqrt(A_Q[:N]**2 + A_Q[N:]**2)
        
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}):")
        print(f"  pinv: A_cs={A_pinv}, A={A_pinv_est}")
        print(f"  Q:    A_cs={A_Q}, A={A_Q_est}")

# Teper: L(T) cherez pryamuyu FOP s ocenkoj A_cs cherez pinv
print("\n=== L(T) cherez pryamuyu FOP s A_cs ot pinv(H) ===")

def compute_signal_from_A_cs(t, A_cs, omega_m, a_m, phi_m, tau):
    """Vosstanovlenie signala iz A_cs"""
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

def likelihood_with_A_cs(xi, t, Z_prime, omega_m, a_m, phi_m, tau, N0, dt, T_dur):
    """Pryamaya FOP s ocenkoj A_cs cherez pinv"""
    A_cs = estimate_A_cs_pinv(Z_prime, omega_m, a_m, phi_m, tau, T_dur)
    s = compute_signal_from_A_cs(t, A_cs, omega_m, a_m, phi_m, tau)
    int_xi_s = np.sum(xi * s) * dt
    int_s2 = np.sum(s**2) * dt
    return (2.0/N0) * int_xi_s - (1.0/N0) * int_s2, A_cs

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau = np.array([t1, t2])
        L, A_cs = likelihood_with_A_cs(xi, t, Z_prime, omega_m, a_m, phi_m, tau, N0, dt, T_dur)
        N = 2
        A_est = np.sqrt(A_cs[:N]**2 + A_cs[N:]**2)
        psi_est = np.arctan2(A_cs[N:], A_cs[:N])
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.10f}, A={A_est}, psi={psi_est}")

# Poisk maksimuma
print("\n=== Poisk maksimuma L(T1,T2) cherez pryamuyu FOP + pinv ===")
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
        L, _ = likelihood_with_A_cs(xi, t, Z_prime, omega_m, a_m, phi_m, tau, N0, dt, T_dur)
        L_grid[i,j] = L

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
print(f"Max at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.10f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"L(true) = {L_grid[i_true, j_true]:.10f}")