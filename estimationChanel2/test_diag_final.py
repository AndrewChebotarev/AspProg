"""
Finalnaya diagnostika: sravnenie L(T) dlya N=1 i N=2 s istinnymi A_cs
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo import compute_matrices_C_S, compute_matrix_H, compute_matrix_Q, compute_sufficient_statistics

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N0=100.0; fs=20000; T1=0.0; T2=1.0

# ===== TEST 1: Signal s 1 luchom, N=1 model =====
print("=== TEST 1: 1-ray signal, N=1 model ===")
A1=0.7; tau1=600e-6; psi1=0.52
data1=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,1,np.array([A1]),np.array([tau1]),np.array([psi1]),fs,T1,T2,40)
Z1,_,_=compute_sufficient_statistics(data1['xi_observed'],data1['t'],omega_m,a_m,N0,fs)

A_cs_true_1 = np.array([A1*np.cos(psi1), A1*np.sin(psi1)])

for t in [500e-6, 600e-6, 700e-6]:
    tau=np.array([t])
    C,S=compute_matrices_C_S(omega_m,phi_m,tau)
    H=compute_matrix_H(C,S)
    Q=compute_matrix_Q(omega_m,a_m,tau,N0)
    M_val = Z1 @ H @ A_cs_true_1
    Q_val = A_cs_true_1 @ Q @ A_cs_true_1
    L_val = M_val - 0.5*Q_val
    print(f"  tau={t*1e6:.0f} us: M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}")

# ===== TEST 2: Signal s 2 luchami, N=1 model =====
print("\n=== TEST 2: 2-ray signal, N=1 model (ocenka 1-go lucha) ===")
data2=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,2,np.array([0.7,0.5]),np.array([600e-6,1200e-6]),np.array([0.52,-0.53]),fs,T1,T2,40)
Z2,_,_=compute_sufficient_statistics(data2['xi_observed'],data2['t'],omega_m,a_m,N0,fs)

# Istinnye A_cs dlya 1-go lucha
A_cs_ray1 = np.array([0.7*np.cos(0.52), 0.7*np.sin(0.52)])

for t in [500e-6, 600e-6, 700e-6]:
    tau=np.array([t])
    C,S=compute_matrices_C_S(omega_m,phi_m,tau)
    H=compute_matrix_H(C,S)
    Q=compute_matrix_Q(omega_m,a_m,tau,N0)
    M_val = Z2 @ H @ A_cs_ray1
    Q_val = A_cs_ray1 @ Q @ A_cs_ray1
    L_val = M_val - 0.5*Q_val
    print(f"  tau={t*1e6:.0f} us: M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}")

# ===== TEST 3: Signal s 2 luchami, N=2 model =====
print("\n=== TEST 3: 2-ray signal, N=2 model ===")
A_cs_true_2 = np.array([0.7*np.cos(0.52), 0.5*np.cos(-0.53), 0.7*np.sin(0.52), 0.5*np.sin(-0.53)])
print(f"A_cs_true = {A_cs_true_2}")

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau=np.array([t1,t2])
        C,S=compute_matrices_C_S(omega_m,phi_m,tau)
        H=compute_matrix_H(C,S)
        Q=compute_matrix_Q(omega_m,a_m,tau,N0)
        M_val = Z2 @ H @ A_cs_true_2
        Q_val = A_cs_true_2 @ Q @ A_cs_true_2
        L_val = M_val - 0.5*Q_val
        print(f"  tau=({t1*1e6:.0f},{t2*1e6:.0f}): M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}")

# ===== TEST 4: Chto esli Z_prime poschitan po drugoy formule? =====
print("\n=== TEST 4: 2-ray signal, N=2 model, Z iz quadrature correlators ===")
# Ispolzuem metod 1: X(tau), Y(tau) cherez sdyig i interpolyaciyu
t = data2['t']
xi = data2['xi_observed']
dt = 1/fs

beta_c = np.zeros_like(t)
beta_s = np.zeros_like(t)
for m in range(M):
    beta_c += a_m[m] * np.cos(omega_m[m] * t - phi_m[m])
    beta_s += a_m[m] * np.sin(omega_m[m] * t - phi_m[m])

def compute_XY(tau_val):
    t_shifted = t - tau_val
    beta_c_shift = np.interp(t, t_shifted, beta_c, left=0, right=0)
    beta_s_shift = np.interp(t, t_shifted, beta_s, left=0, right=0)
    X = (2.0/N0) * np.sum(xi * beta_c_shift) * dt
    Y = (2.0/N0) * np.sum(xi * beta_s_shift) * dt
    return X, Y

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        X1,Y1 = compute_XY(t1)
        X2,Y2 = compute_XY(t2)
        Z = np.array([X1, X2, Y1, Y2])
        tau=np.array([t1,t2])
        Q=compute_matrix_Q(omega_m,a_m,tau,N0)
        M_val = Z @ A_cs_true_2
        Q_val = A_cs_true_2 @ Q @ A_cs_true_2
        L_val = M_val - 0.5*Q_val
        print(f"  tau=({t1*1e6:.0f},{t2*1e6:.0f}): M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}")