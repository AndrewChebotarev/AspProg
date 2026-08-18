"""
Proverka: nuzhen li mnozhitel T_dur/2 v Q-matricu?
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo import compute_matrices_C_S, compute_matrix_H, compute_sufficient_statistics

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N0=100.0; fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,2,np.array([0.7,0.5]),np.array([600e-6,1200e-6]),np.array([0.52,-0.53]),fs,T1,T2,40)
Z,_,_=compute_sufficient_statistics(data['xi_observed'],data['t'],omega_m,a_m,N0,fs)

A_cs_true = np.array([0.7*np.cos(0.52), 0.5*np.cos(-0.53), 0.7*np.sin(0.52), 0.5*np.sin(-0.53)])

def compute_Q_with_factor(omega_m, a_m, tau, N0, factor=1.0):
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
                sum_sin += weight * np.sin(omega_m[m] * (tau[k] - tau[i]))
            Qc[i, k] = factor * sum_cos
            Qcs[i, k] = factor * sum_sin
    Q = np.vstack([np.hstack([Qc, Qcs]), np.hstack([Qcs.T, Qc])])
    return Q

print("=== Vliyanie mnozhitelya T_dur/2 na L(T) ===")
for factor_name, factor in [("bez (1.0)", 1.0), ("T_dur/2 (0.5)", T_dur/2), ("T_dur (1.0)", T_dur)]:
    print(f"\n--- Mnozhitel Q = {factor_name} ---")
    for t1 in [500e-6, 600e-6, 700e-6]:
        for t2 in [1100e-6, 1200e-6, 1300e-6]:
            tau=np.array([t1,t2])
            C,S=compute_matrices_C_S(omega_m,phi_m,tau)
            H=compute_matrix_H(C,S)
            Q=compute_Q_with_factor(omega_m,a_m,tau,N0,factor)
            M_val = Z @ H @ A_cs_true
            Q_val = A_cs_true @ Q @ A_cs_true
            L_val = M_val - 0.5*Q_val
            print(f"  tau=({t1*1e6:.0f},{t2*1e6:.0f}): M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}")

# Proverim: chto esli v Z_prime dobavit T_dur?
print("\n\n=== S Z_prime * T_dur ===")
Z_scaled = Z * T_dur
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau=np.array([t1,t2])
        C,S=compute_matrices_C_S(omega_m,phi_m,tau)
        H=compute_matrix_H(C,S)
        Q=compute_Q_with_factor(omega_m,a_m,tau,N0,1.0)
        M_val = Z_scaled @ H @ A_cs_true
        Q_val = A_cs_true @ Q @ A_cs_true
        L_val = M_val - 0.5*Q_val
        print(f"  tau=({t1*1e6:.0f},{t2*1e6:.0f}): M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}")

# Proverim: chto esli v Z_prime dobavit 2/T_dur?
print("\n\n=== S Z_prime * 2/T_dur ===")
Z_scaled2 = Z * (2/T_dur)
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau=np.array([t1,t2])
        C,S=compute_matrices_C_S(omega_m,phi_m,tau)
        H=compute_matrix_H(C,S)
        Q=compute_Q_with_factor(omega_m,a_m,tau,N0,1.0)
        M_val = Z_scaled2 @ H @ A_cs_true
        Q_val = A_cs_true @ Q @ A_cs_true
        L_val = M_val - 0.5*Q_val
        print(f"  tau=({t1*1e6:.0f},{t2*1e6:.0f}): M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}")