"""
Sovmestnaya ocenka cherez posledovatelnoe vychitanie:
dlya kazhdogo kandidata tau2: ocenit 2-y luch, vychest, ocenit 1-y, poschitat FOP
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; dt=1/fs; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
t=data['t']
xi=data['xi_observed'].copy()

def estimate_1ray(xi, t, tau_search, omega_m, a_m, phi_m, N0, dt, T_dur):
    """Ocenka 1 lucha po signalu cherez (99) + pryamaya FOP"""
    Z_prime,_,_ = compute_sufficient_statistics(xi, t, omega_m, a_m, N0, fs)
    best_L = -np.inf
    best_tau = None
    best_A = None
    best_psi = None
    
    for tau_test in tau_search:
        tau = np.array([tau_test])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        Q_inv = np.linalg.inv(Q)
        A_cs = Z_prime @ H @ Q_inv
        A_cs = A_cs * (2.0 / T_dur)
        
        A_est = np.sqrt(A_cs[0]**2 + A_cs[1]**2)
        psi_est = np.arctan2(A_cs[1], A_cs[0])
        
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

# Algoritm: dlya kazhdogo kandidata tau2, vychitaem 2-y luch i ocenivaem 1-y
print("Sovmestnaya ocenka cherez posledovatelnoe vychitanie:")
print("=" * 70)

tau2_search = np.linspace(800e-6, 1600e-6, 81)
tau1_search = np.linspace(200e-6, 1000e-6, 81)

best_total_L = -np.inf
best_result = None

for tau2_test in tau2_search:
    # Ocenivaem 2-y luch iz polnogo signala
    tau2_arr = np.array([tau2_test])
    C2, S2 = compute_matrices_C_S(omega_m, phi_m, tau2_arr)
    H2 = compute_matrix_H(C2, S2)
    Q2 = compute_matrix_Q(omega_m, a_m, tau2_arr, N0)
    
    Z_prime_all,_,_ = compute_sufficient_statistics(xi, t, omega_m, a_m, N0, fs)
    Q2_inv = np.linalg.inv(Q2)
    A_cs2 = Z_prime_all @ H2 @ Q2_inv
    A_cs2 = A_cs2 * (2.0 / T_dur)
    A2_test = np.sqrt(A_cs2[0]**2 + A_cs2[1]**2)
    psi2_test = np.arctan2(A_cs2[1], A_cs2[0])
    
    # Vychitaem 2-y luch
    s2 = subtract_ray(t, A2_test, tau2_test, psi2_test)
    xi_res = xi - s2
    
    # Ocenivaem 1-y luch iz ostatka
    tau1_test, A1_test, psi1_test, L1 = estimate_1ray(xi_res, t, tau1_search, omega_m, a_m, phi_m, N0, dt, T_dur)
    
    # Polnaya FOP
    s_total = subtract_ray(t, A1_test, tau1_test, psi1_test) + s2
    int_xi_s = np.sum(xi * s_total) * dt
    int_s2 = np.sum(s_total**2) * dt
    L_total = (2.0/N0) * int_xi_s - (1.0/N0) * int_s2
    
    if L_total > best_total_L:
        best_total_L = L_total
        best_result = (tau1_test, tau2_test, A1_test, A2_test, psi1_test, psi2_test)
    
    if abs(tau2_test*1e6 - 1200) < 5:
        print(f"tau2={tau2_test*1e6:.0f}: tau1={tau1_test*1e6:.1f}, A1={A1_test:.3f}, A2={A2_test:.3f}, L={L_total:.8f}")

tau1_opt, tau2_opt, A1_opt, A2_opt, psi1_opt, psi2_opt = best_result
print(f"\nOptimalnye parametry:")
print(f"  Luch 1: tau={tau1_opt*1e6:.1f} mks (true={tau_true[0]*1e6:.0f}), A={A1_opt:.3f} (true={A_true[0]:.1f}), psi={psi1_opt:.3f} (true={psi_true[0]:.2f})")
print(f"  Luch 2: tau={tau2_opt*1e6:.1f} mks (true={tau_true[1]*1e6:.0f}), A={A2_opt:.3f} (true={A_true[1]:.1f}), psi={psi2_opt:.3f} (true={psi_true[1]:.2f})")
print(f"  L_total = {best_total_L:.10f}")

print(f"\nOshibki:")
print(f"  tau: {(tau1_opt-tau_true[0])*1e6:.2f}, {(tau2_opt-tau_true[1])*1e6:.2f} mks")
print(f"  A:   {(A1_opt-A_true[0])/A_true[0]*100:.2f}, {(A2_opt-A_true[1])/A_true[1]*100:.2f} %")
print(f"  psi: {psi1_opt-psi_true[0]:.4f}, {psi2_opt-psi_true[1]:.4f} rad")