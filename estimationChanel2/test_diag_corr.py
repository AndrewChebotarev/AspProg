"""
Ispolzuem korrelyaciyu mezhdu xi(t) i s(t) kak kriteriy dlya N=2
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

Z_prime,_,_ = compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

def compute_criteria(tau, xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur):
    """Vychislenie raznyh kriteriev dlya kandidata (tau1, tau2)"""
    N = len(tau)
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    
    A_c = A_cs[:N]
    A_s = A_cs[N:]
    A = np.sqrt(A_c**2 + A_s**2)
    psi = np.arctan2(A_s, A_c)
    
    s = np.zeros_like(t)
    for k in range(N):
        for m in range(M):
            s += A[k] * a_m[m] * np.cos(omega_m[m] * (t - tau[k]) - phi_m[m] - psi[k])
    
    # Kriteriy 1: Pryamaya FOP
    int_xi_s = np.sum(xi * s) * dt
    int_s2 = np.sum(s**2) * dt
    L_fop = (2.0/N0) * int_xi_s - (1.0/N0) * int_s2
    
    # Kriteriy 2: Korrelyaciya mezhdu xi i s
    xi_mean = xi - np.mean(xi)
    s_mean = s - np.mean(s)
    corr = np.sum(xi_mean * s_mean) / np.sqrt(np.sum(xi_mean**2) * np.sum(s_mean**2))
    
    # Kriteriy 3: Otnoshenie pravdepodobiya (chi-square)
    residual = xi - s
    chi2 = np.sum(residual**2) * dt
    
    # Kriteriy 4: MSE vosstanovleniya
    mse = np.mean((xi - s)**2)
    
    return L_fop, corr, chi2, mse, A, psi

print("=== Sravnenie kriteriev dlya raznyh (tau1, tau2) ===")
print(f"{'tau1':>6} {'tau2':>6} | {'L_FOP':>12} {'Corr':>8} {'chi2':>12} {'MSE':>12} | {'A1':>6} {'A2':>6}")
print("-" * 85)

for t1 in [400e-6, 500e-6, 600e-6, 700e-6, 800e-6]:
    for t2 in [1000e-6, 1100e-6, 1200e-6, 1300e-6, 1400e-6]:
        if t2 <= t1:
            continue
        tau = np.array([t1, t2])
        L_fop, corr, chi2, mse, A, psi = compute_criteria(tau, xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur)
        print(f"{t1*1e6:6.0f} {t2*1e6:6.0f} | {L_fop:12.8f} {corr:8.6f} {chi2:12.8f} {mse:12.8f} | {A[0]:6.3f} {A[1]:6.3f}")

# Poisk maksimuma po korrelyacii
print("\n=== Poisk maksimuma po korrelyacii ===")
N_grid = 51
tau1_range = np.linspace(0, 2000e-6, N_grid)
tau2_range = np.linspace(0, 2000e-6, N_grid)
corr_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            corr_grid[i,j] = -1
            continue
        tau = np.array([t1, t2])
        _, corr, _, _, _, _ = compute_criteria(tau, xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur)
        corr_grid[i,j] = corr

idx_max = np.unravel_index(np.argmax(corr_grid), corr_grid.shape)
print(f"Max corr at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"corr_max = {corr_grid[idx_max]:.6f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"corr(true) = {corr_grid[i_true, j_true]:.6f}")

# Poisk minimuma po MSE
print("\n=== Poisk minimuma po MSE ===")
mse_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            mse_grid[i,j] = np.inf
            continue
        tau = np.array([t1, t2])
        _, _, _, mse, _, _ = compute_criteria(tau, xi, t, Z_prime, omega_m, a_m, phi_m, N0, dt, T_dur)
        mse_grid[i,j] = mse

idx_min = np.unravel_index(np.argmin(mse_grid), mse_grid.shape)
print(f"Min MSE at: tau1={tau1_range[idx_min[0]]*1e6:.1f} us, tau2={tau2_range[idx_min[1]]*1e6:.1f} us")
print(f"MSE_min = {mse_grid[idx_min]:.10f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"MSE(true) = {mse_grid[i_true, j_true]:.10f}")