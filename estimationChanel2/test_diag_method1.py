"""
Proverka: metod 1 (kvadraturnye korrelyatory) dlya N=2
Ispolzuem formulu (98): L = Z @ A_cs^T - 0.5 * A_cs @ Q @ A_cs^T
gde Z = [X(tau1), X(tau2), Y(tau1), Y(tau2)] - kvadraturnye korrelyatory
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo import compute_matrix_Q

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N0=100.0; fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,2,np.array([0.7,0.5]),np.array([600e-6,1200e-6]),np.array([0.52,-0.53]),fs,T1,T2,40)
t=data['t']
xi=data['xi_observed']
dt=1/fs

# Opornye signaly (formuly 29-30)
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

A_cs_true = np.array([0.7*np.cos(0.52), 0.5*np.cos(-0.53), 0.7*np.sin(0.52), 0.5*np.sin(-0.53)])

print("=== Metod 1: L(T) cherez kvadraturnye korrelyatory ===")
print("Formula: L = Z @ A_cs^T - 0.5 * A_cs @ Q @ A_cs^T")
print()

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        X1,Y1 = compute_XY(t1)
        X2,Y2 = compute_XY(t2)
        Z = np.array([X1, X2, Y1, Y2])
        tau=np.array([t1,t2])
        Q=compute_matrix_Q(omega_m,a_m,tau,N0)
        M_val = Z @ A_cs_true
        Q_val = A_cs_true @ Q @ A_cs_true
        L_val = M_val - 0.5*Q_val
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): X1={X1:.6f}, Y1={Y1:.6f}, X2={X2:.6f}, Y2={Y2:.6f}, L={L_val:.6f}")

# Teper ispolzuem formulu (100): L = 0.5 * Z @ Q^{-1} @ Z^T
print("\n=== Metod 1: L(T) cherez formulu (100) ===")
print("Formula: L = 0.5 * Z @ Q^{-1} @ Z^T")
print()

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        X1,Y1 = compute_XY(t1)
        X2,Y2 = compute_XY(t2)
        Z = np.array([X1, X2, Y1, Y2])
        tau=np.array([t1,t2])
        Q=compute_matrix_Q(omega_m,a_m,tau,N0)
        try:
            Q_inv = np.linalg.inv(Q)
            L_val = 0.5 * Z @ Q_inv @ Z.T
        except:
            L_val = -np.inf
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L_val:.6f}")

# Poisk maksimuma L(T1,T2) po setke
print("\n=== Poisk maksimuma L(T1,T2) metod 1, formula (100) ===")
N_grid = 51
tau1_range = np.linspace(0, 2000e-6, N_grid)
tau2_range = np.linspace(0, 2000e-6, N_grid)
L_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            L_grid[i,j] = -np.inf
            continue
        X1,Y1 = compute_XY(t1)
        X2,Y2 = compute_XY(t2)
        Z = np.array([X1, X2, Y1, Y2])
        tau=np.array([t1,t2])
        Q=compute_matrix_Q(omega_m,a_m,tau,N0)
        try:
            Q_inv = np.linalg.inv(Q)
            L_grid[i,j] = 0.5 * Z @ Q_inv @ Z.T
        except:
            L_grid[i,j] = -np.inf

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
print(f"Max at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.6f}")

# L pri istinnyh
i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"L(true) = {L_grid[i_true, j_true]:.6f}")