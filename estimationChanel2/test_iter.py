"""
Iterativnoe utochnenie: posle ocenki 2-go lucha, pereocenit 1-y s vychitaniem 2-go
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; dt=1/fs; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

def estimate_1ray_from_signal(xi, t, tau_search_range, omega_m, a_m, phi_m, N0, dt, T_dur):
    """Ocenka 1 lucha po signalu cherez pryamuyu FOP"""
    Z_prime,_,_ = compute_sufficient_statistics(xi, t, omega_m, a_m, N0, fs)
    best_L = -np.inf
    best_tau = None
    best_A = None
    best_psi = None
    
    for tau_test in tau_search_range:
        tau = np.array([tau_test])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        Q_inv = np.linalg.inv(Q)
        A_cs = Z_prime @ H @ Q_inv
        A_cs = A_cs * (2.0 / T_dur)
        
        A_est = np.sqrt(A_cs[0]**2 + A_cs[1]**2)
        psi_est = np.arctan2(A_cs[1], A_cs[0])
        
        # Pryamaya FOP
        s = np.zeros_like(t)
        for m in range(M):
            s += A_est * a_m[m] * np.cos(omega_m[m] * (t - tau_test) - phi_m[m] - psi_est)
        
        int_xi_s = np.sum(xi * s) * dt
        int_s2 = np.sum(s**2) * dt
        L = (2.0/N0) * int_xi_s - (1.0/N0) * int_s2
        
        if L > best_L:
            best_L = L
            best_tau = tau_test
            best_A = A_est
            best_psi = psi_est
    
    return best_tau, best_A, best_psi, best_L

def subtract_ray(t, A, tau, psi):
    s = np.zeros_like(t)
    for m in range(M):
        s += A * a_m[m] * np.cos(omega_m[m] * (t - tau) - phi_m[m] - psi)
    return s

# Iterativnaya ocenka
print("Iterativnaya ocenka:")
print("=" * 60)

# Initialnye ocenki
tau_search_1 = np.linspace(400e-6, 800e-6, 201)
tau_search_2 = np.linspace(1000e-6, 1400e-6, 201)

for iteration in range(5):
    print(f"\nIteraciya {iteration + 1}:")
    
    if iteration == 0:
        xi_curr = xi.copy()
    else:
        # Vychitaem 2-y luch dlya ocenki 1-go
        xi_curr = xi - subtract_ray(t, A2, tau2, psi2)
    
    # Ocenka 1-go lucha
    tau1, A1, psi1, L1 = estimate_1ray_from_signal(xi_curr, t, tau_search_1, omega_m, a_m, phi_m, N0, dt, T_dur)
    print(f"  Luch 1: tau={tau1*1e6:.1f} mks, A={A1:.4f}, psi={psi1:.4f}, L={L1:.10f}")
    
    # Vychitaem 1-y luch
    xi_res = xi - subtract_ray(t, A1, tau1, psi1)
    
    # Ocenka 2-go lucha
    tau2, A2, psi2, L2 = estimate_1ray_from_signal(xi_res, t, tau_search_2, omega_m, a_m, phi_m, N0, dt, T_dur)
    print(f"  Luch 2: tau={tau2*1e6:.1f} mks, A={A2:.4f}, psi={psi2:.4f}, L={L2:.10f}")

print(f"\n" + "=" * 60)
print(f"Konechnyy rezultat:")
print(f"  Luch 1: tau={tau1*1e6:.1f} mks (true={tau_true[0]*1e6:.0f}), A={A1:.3f} (true={A_true[0]:.1f}), psi={psi1:.3f} (true={psi_true[0]:.2f})")
print(f"  Luch 2: tau={tau2*1e6:.1f} mks (true={tau_true[1]*1e6:.0f}), A={A2:.3f} (true={A_true[1]:.1f}), psi={psi2:.3f} (true={psi_true[1]:.2f})")
print(f"\nOshibki:")
print(f"  tau: {(tau1-tau_true[0])*1e6:.2f}, {(tau2-tau_true[1])*1e6:.2f} mks")
print(f"  A:   {(A1-A_true[0])/A_true[0]*100:.2f}, {(A2-A_true[1])/A_true[1]*100:.2f} %")
print(f"  psi: {psi1-psi_true[0]:.4f}, {psi2-psi_true[1]:.4f} rad")