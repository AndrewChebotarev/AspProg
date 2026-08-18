"""
Diagnostika: dobavlyaem regulyarizaciyu Q-matricy
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

def compute_likelihood_reg(Z_prime, H, Q, reg=1e-6):
    try:
        Q_reg = Q + reg * np.eye(len(Q))
        Q_inv = np.linalg.inv(Q_reg)
        return 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    except:
        return -np.inf

# Sravnenie L(T) s regulyarizaciey i bez
print("=== Sravnenie L(T) dlya N=2 ===")
for reg in [0, 1e-10, 1e-8, 1e-6, 1e-4]:
    tau_range = np.linspace(0, 2000e-6, 201)
    L_vals = []
    for t1 in tau_range:
        for t2 in tau_range:
            if t2 <= t1:
                continue
            tau = np.array([t1, t2])
            C, S = compute_matrices_C_S(omega_m, phi_m, tau)
            H = compute_matrix_H(C, S)
            Q = compute_matrix_Q(omega_m, a_m, tau, N0)
            L_vals.append(compute_likelihood_reg(Z_prime, H, Q, reg))
    L_vals = np.array(L_vals)
    print(f"reg={reg}: max L = {np.max(L_vals):.6f}")

# Detalno: L(T) s regulyarizaciey v istinnyh i lozhnyh tochkah
print("\n=== Detalno s reg=1e-6 ===")
reg = 1e-6

# Istinnye
tau_true_arr = np.array([600e-6, 1200e-6])
C, S = compute_matrices_C_S(omega_m, phi_m, tau_true_arr)
H = compute_matrix_H(C, S)
Q = compute_matrix_Q(omega_m, a_m, tau_true_arr, N0)
L_true = compute_likelihood_reg(Z_prime, H, Q, reg)
print(f"L(true) = {L_true:.6f}")

# Lozhnye (40, 60)
tau_test = np.array([40e-6, 60e-6])
C2, S2 = compute_matrices_C_S(omega_m, phi_m, tau_test)
H2 = compute_matrix_H(C2, S2)
Q2 = compute_matrix_Q(omega_m, a_m, tau_test, N0)
L_test = compute_likelihood_reg(Z_prime, H2, Q2, reg)
print(f"L(40,60) = {L_test:.6f}")

# Poisk maksimuma s regulyarizaciey
print("\n=== Poisk maksimuma L(T1,T2) s reg=1e-6 ===")
N_grid = 101
tau1_range = np.linspace(0, 2000e-6, N_grid)
tau2_range = np.linspace(0, 2000e-6, N_grid)
L_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            L_grid[i, j] = -np.inf
            continue
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L_grid[i, j] = compute_likelihood_reg(Z_prime, H, Q, reg)

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
t1_est = tau1_range[idx_max[0]]
t2_est = tau2_range[idx_max[1]]
print(f"Max at: tau1={t1_est*1e6:.1f} us, tau2={t2_est*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.6f}")

# Proverim: chto esli ispolzovat tolko diagonal Q?
print("\n=== Ispolzuem tolko diagonal Q ===")
def compute_likelihood_diagQ(Z_prime, H, Q):
    try:
        Q_diag = np.diag(np.diag(Q))
        Q_inv = np.linalg.inv(Q_diag)
        return 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    except:
        return -np.inf

L_true_diag = compute_likelihood_diagQ(Z_prime, H, Q)
print(f"L(true) diagQ = {L_true_diag:.6f}")

C2, S2 = compute_matrices_C_S(omega_m, phi_m, tau_test)
H2 = compute_matrix_H(C2, S2)
Q2 = compute_matrix_Q(omega_m, a_m, tau_test, N0)
L_test_diag = compute_likelihood_diagQ(Z_prime, H2, Q2)
print(f"L(40,60) diagQ = {L_test_diag:.6f}")

# Poisk s diagQ
L_grid_diag = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            L_grid_diag[i, j] = -np.inf
            continue
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L_grid_diag[i, j] = compute_likelihood_diagQ(Z_prime, H, Q)

idx_max_diag = np.unravel_index(np.argmax(L_grid_diag), L_grid_diag.shape)
print(f"Max diagQ at: tau1={tau1_range[idx_max_diag[0]]*1e6:.1f} us, tau2={tau2_range[idx_max_diag[1]]*1e6:.1f} us")
print(f"L_max diagQ = {L_grid_diag[idx_max_diag]:.6f}")