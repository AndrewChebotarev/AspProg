"""
Diagnostika: proverka Q-matricy i L(T) dlya N=2
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)
print(f"Z_prime = {Z_prime}")

# Proverim Q-matricu pri istinnyh znacheniyah
tau_true_arr = np.array([600e-6, 1200e-6])
C, S = compute_matrices_C_S(omega_m, phi_m, tau_true_arr)
H = compute_matrix_H(C, S)
Q = compute_matrix_Q(omega_m, a_m, tau_true_arr, N0)

print(f"\n=== Pri istinnyh znacheniyah (tau1=600, tau2=1200 us) ===")
print(f"C = \n{C}")
print(f"S = \n{S}")
print(f"H = \n{H}")
print(f"Q = \n{Q}")
print(f"det(Q) = {np.linalg.det(Q):.6f}")
print(f"cond(Q) = {np.linalg.cond(Q):.6f}")

Q_inv = np.linalg.inv(Q)
L_true = 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
print(f"L(true) = {L_true:.6f}")

# Proverim pri (40, 60) us
tau_test = np.array([40e-6, 60e-6])
C2, S2 = compute_matrices_C_S(omega_m, phi_m, tau_test)
H2 = compute_matrix_H(C2, S2)
Q2 = compute_matrix_Q(omega_m, a_m, tau_test, N0)
print(f"\n=== Pri (tau1=40, tau2=60 us) ===")
print(f"Q = \n{Q2}")
print(f"det(Q) = {np.linalg.det(Q2):.6f}")
print(f"cond(Q) = {np.linalg.cond(Q2):.6f}")

Q2_inv = np.linalg.inv(Q2)
L_test = 0.5 * Z_prime @ H2 @ Q2_inv @ H2.T @ Z_prime.T
print(f"L(test) = {L_test:.6f}")

# Chto esli N0=1.0?
print("\n\n=== S N0=1.0 ===")
Z_prime1,_,_ = compute_sufficient_statistics(xi,t,omega_m,a_m,1.0,fs)
print(f"Z_prime = {Z_prime1}")

Q_true1 = compute_matrix_Q(omega_m, a_m, tau_true_arr, 1.0)
Q_test1 = compute_matrix_Q(omega_m, a_m, tau_test, 1.0)
L_true1 = 0.5 * Z_prime1 @ H @ np.linalg.inv(Q_true1) @ H.T @ Z_prime1.T
L_test1 = 0.5 * Z_prime1 @ H2 @ np.linalg.inv(Q_test1) @ H2.T @ Z_prime1.T
print(f"L(true) = {L_true1:.6f}")
print(f"L(test) = {L_test1:.6f}")

# Chto esli N0=0.1?
print("\n\n=== S N0=0.1 ===")
Z_prime01,_,_ = compute_sufficient_statistics(xi,t,omega_m,a_m,0.1,fs)
print(f"Z_prime = {Z_prime01}")

Q_true01 = compute_matrix_Q(omega_m, a_m, tau_true_arr, 0.1)
Q_test01 = compute_matrix_Q(omega_m, a_m, tau_test, 0.1)
L_true01 = 0.5 * Z_prime01 @ H @ np.linalg.inv(Q_true01) @ H.T @ Z_prime01.T
L_test01 = 0.5 * Z_prime01 @ H2 @ np.linalg.inv(Q_test01) @ H2.T @ Z_prime01.T
print(f"L(true) = {L_true01:.6f}")
print(f"L(test) = {L_test01:.6f}")

# Proverim: chto esli ispolzovat Z_prime iz MonteCarlo.py?
print("\n\n=== Proverka s MonteCarlo.py (N=1) ===")
from MonteCarlo import compute_sufficient_statistics as mc_stats
Z_prime_mc,_,_ = mc_stats(xi, t, omega_m, a_m, N0, fs)
print(f"Z_prime (MC) = {Z_prime_mc}")

# Dlya N=1
tau_1 = np.array([600e-6])
C1, S1 = compute_matrices_C_S(omega_m, phi_m, tau_1)
H1 = compute_matrix_H(C1, S1)
Q1 = compute_matrix_Q(omega_m, a_m, tau_1, N0)
L_1 = 0.5 * Z_prime_mc @ H1 @ np.linalg.inv(Q1) @ H1.T @ Z_prime_mc.T
print(f"L(N=1, tau=600) = {L_1:.6f}")

tau_1_test = np.array([56e-6])
C1t, S1t = compute_matrices_C_S(omega_m, phi_m, tau_1_test)
H1t = compute_matrix_H(C1t, S1t)
Q1t = compute_matrix_Q(omega_m, a_m, tau_1_test, N0)
L_1t = 0.5 * Z_prime_mc @ H1t @ np.linalg.inv(Q1t) @ H1t.T @ Z_prime_mc.T
print(f"L(N=1, tau=56) = {L_1t:.6f}")

# KLUCH: sravnenie MonteCarlo.py i MonteCarlo2.py
print("\n\n=== Sravnenie funkciy MonteCarlo.py vs MonteCarlo2.py ===")
from MonteCarlo import compute_likelihood as mc_likelihood
from MonteCarlo import compute_matrix_Q as mc_Q
from MonteCarlo import compute_matrices_C_S as mc_CS
from MonteCarlo import compute_matrix_H as mc_H

C_mc, S_mc = mc_CS(omega_m, phi_m, tau_1)
H_mc = mc_H(C_mc, S_mc)
Q_mc = mc_Q(omega_m, a_m, tau_1, N0)
L_mc = mc_likelihood(Z_prime_mc, H_mc, Q_mc)
print(f"MonteCarlo.py: L(N=1, tau=600) = {L_mc:.6f}")

# Te zhe funkcii iz MonteCarlo2.py
C2_mc, S2_mc = compute_matrices_C_S(omega_m, phi_m, tau_1)
H2_mc = compute_matrix_H(C2_mc, S2_mc)
Q2_mc = compute_matrix_Q(omega_m, a_m, tau_1, N0)
L2_mc = compute_likelihood(Z_prime_mc, H2_mc, Q2_mc)
print(f"MonteCarlo2.py: L(N=1, tau=600) = {L2_mc:.6f}")

# Sravnivaem Q-matricy
print(f"\nQ (MC.py):\n{Q_mc}")
print(f"Q (MC2.py):\n{Q2_mc}")
print(f"Q diff: {np.max(np.abs(Q_mc - Q2_mc)):.15f}")