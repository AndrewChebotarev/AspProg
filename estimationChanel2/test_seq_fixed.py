"""
Posledovatelnaya ocenka s ispravlennoy Q
"""
import numpy as np
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

def estimate_1ray_fixed(xi, t, Z_prime, omega_m, a_m, phi_m, tau_search, N0, dt, T_dur):
    """Ocenka 1 lucha s ispravlennoy Q"""
    best_L = -np.inf
    best_tau = None
    best_A_cs = None
    
    for tau_test in tau_search:
        tau = np.array([tau_test])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_Q_fixed(omega_m, a_m, tau, N0)
        Q_inv = np.linalg.inv(Q)
        A_cs = Z_prime @ H @ Q_inv
        A_cs = A_cs * (2.0 / T_dur)
        
        # L cherez (98)
        M_val = Z_prime @ H @ A_cs
        Q_val = A_cs @ Q @ A_cs
        L = M_val - 0.5 * Q_val
        
        if L > best_L:
            best_L = L
            best_tau = tau_test
            best_A_cs = A_cs.copy()
    
    N = 1
    A_est = np.sqrt(best_A_cs[0]**2 + best_A_cs[1]**2)
    psi_est = np.arctan2(best_A_cs[1], best_A_cs[0])
    
    return best_tau, A_est, psi_est, best_L

def subtract_ray_fixed(t, A, tau, psi):
    s = np.zeros_like(t)
    for m in range(M):
        s += A * a_m[m] * np.cos(omega_m[m] * (t - tau) - phi_m[m] - psi)
    return s

Z_prime,_,_ = compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

# Shag 1: ocenka 1-go lucha
print("Shag 1: ocenka 1-go lucha...")
tau_search = np.linspace(400e-6, 800e-6, 201)
tau1, A1, psi1, L1 = estimate_1ray_fixed(xi, t, Z_prime, omega_m, a_m, phi_m, tau_search, N0, dt, T_dur)
print(f"  tau1={tau1*1e6:.1f} us (true={tau_true[0]*1e6:.0f}), A1={A1:.3f} (true={A_true[0]:.1f}), psi1={psi1:.3f} (true={psi_true[0]:.2f}), L={L1:.10f}")

# Shag 2: vychitaem 1-y luch
s1 = subtract_ray_fixed(t, A1, tau1, psi1)
xi_res = xi - s1

# Shag 3: ocenka 2-go lucha
print("Shag 2: ocenka 2-go lucha...")
Z_prime2,_,_ = compute_sufficient_statistics(xi_res,t,omega_m,a_m,N0,fs)
tau_search2 = np.linspace(1000e-6, 1400e-6, 201)
tau2, A2, psi2, L2 = estimate_1ray_fixed(xi_res, t, Z_prime2, omega_m, a_m, phi_m, tau_search2, N0, dt, T_dur)
print(f"  tau2={tau2*1e6:.1f} us (true={tau_true[1]*1e6:.0f}), A2={A2:.3f} (true={A_true[1]:.1f}), psi2={psi2:.3f} (true={psi_true[1]:.2f}), L={L2:.10f}")

print(f"\nOshibki:")
print(f"  tau: {(tau1-tau_true[0])*1e6:.2f}, {(tau2-tau_true[1])*1e6:.2f} us")
print(f"  A:   {(A1-A_true[0])/A_true[0]*100:.2f}, {(A2-A_true[1])/A_true[1]*100:.2f} %")
print(f"  psi: {psi1-psi_true[0]:.4f}, {psi2-psi_true[1]:.4f} rad")