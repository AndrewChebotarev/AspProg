"""
Proverka: Z' @ H = A_cs @ Q dlya N=1, no ne dlya N=2?
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo import compute_matrices_C_S, compute_matrix_H, compute_matrix_Q, compute_sufficient_statistics

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N0=100.0; fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1

# Signal s 1 luchom
data1=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,1,np.array([0.7]),np.array([600e-6]),np.array([0.52]),fs,T1,T2,40)
Z1,_,_=compute_sufficient_statistics(data1['xi_observed'],data1['t'],omega_m,a_m,N0,fs)

tau1=np.array([600e-6])
C1,S1=compute_matrices_C_S(omega_m,phi_m,tau1)
H1=compute_matrix_H(C1,S1)
Q1=compute_matrix_Q(omega_m,a_m,tau1,N0)

A_cs_1 = np.array([0.7*np.cos(0.52), 0.7*np.sin(0.52)])

print("=== N=1: proverka Z' @ H == A_cs @ Q ===")
print(f"Z' @ H = {Z1 @ H1}")
print(f"A_cs @ Q = {A_cs_1 @ Q1}")
print(f"Raznica = {np.max(np.abs(Z1 @ H1 - A_cs_1 @ Q1)):.15f}")

# Signal s 2 luchami
data2=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,2,np.array([0.7,0.5]),np.array([600e-6,1200e-6]),np.array([0.52,-0.53]),fs,T1,T2,40)
Z2,_,_=compute_sufficient_statistics(data2['xi_observed'],data2['t'],omega_m,a_m,N0,fs)

tau2=np.array([600e-6, 1200e-6])
C2,S2=compute_matrices_C_S(omega_m,phi_m,tau2)
H2=compute_matrix_H(C2,S2)
Q2=compute_matrix_Q(omega_m,a_m,tau2,N0)

A_cs_2 = np.array([0.7*np.cos(0.52), 0.5*np.cos(-0.53), 0.7*np.sin(0.52), 0.5*np.sin(-0.53)])

print("\n=== N=2: proverka Z' @ H == A_cs @ Q ===")
print(f"Z' @ H = {Z2 @ H2}")
print(f"A_cs @ Q = {A_cs_2 @ Q2}")
print(f"Raznica = {np.max(np.abs(Z2 @ H2 - A_cs_2 @ Q2)):.15f}")

# Chto esli Z' poschitat po drugomu?
print("\n=== N=2: Z' @ H po komponentam ===")
print(f"Z' = {Z2}")
print(f"H = \n{H2}")
print(f"Z' @ H = {Z2 @ H2}")

# Proverim: Z' @ H dolzhno byt ravno A_cs @ Q
# A_cs @ Q = 
print(f"\nA_cs = {A_cs_2}")
print(f"Q = \n{Q2}")
print(f"A_cs @ Q = {A_cs_2 @ Q2}")

# Esli by Z' byl tochnym, to Z' = A_cs @ Q @ H^{-1}
# No H ne kvadratnaya (4x4), tak chto nuzhno H^+
print("\n=== Proverka: Z' cherez psevdoobratenie H ===")
H_pinv = np.linalg.pinv(H2)
Z_from_A = A_cs_2 @ Q2 @ H_pinv
print(f"Z' iz A_cs = {Z_from_A}")
print(f"Z' iz signala = {Z2}")
print(f"Raznica = {np.max(np.abs(Z_from_A - Z2)):.15f}")

# Esche odna proverka: formula (97) M = Z' @ H @ A_cs^T
# i formula (49) Q_val = A_cs @ Q @ A_cs^T
print("\n=== Proverka formul (97) i (49) ===")
M_val = Z2 @ H2 @ A_cs_2
Q_val = A_cs_2 @ Q2 @ A_cs_2
print(f"M = Z' @ H @ A_cs^T = {M_val:.10f}")
print(f"Q_val = A_cs @ Q @ A_cs^T = {Q_val:.10f}")
print(f"L = M - 0.5*Q = {M_val - 0.5*Q_val:.10f}")

# Dlya N=1
M1 = Z1 @ H1 @ A_cs_1
Q1_val = A_cs_1 @ Q1 @ A_cs_1
print(f"\nN=1: M = {M1:.10f}, Q = {Q1_val:.10f}, L = {M1 - 0.5*Q1_val:.10f}")

# Teoreticheski: dlya N=1, Z' @ H = A_cs @ Q, poetomu M = A_cs @ Q @ A_cs^T = Q_val
# i L = Q_val - 0.5*Q_val = 0.5*Q_val
# No my vidim L=0, znachit chto-to ne tak
print(f"\nN=1: 0.5*Q = {0.5*Q1_val:.10f}")
print(f"N=1: M - 0.5*Q = {M1 - 0.5*Q1_val:.10f}")