"""
Polnaya proverka: L(T) + A_cs otsenka dlya M=3, N=2
"""
import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import compute_matrices_C_S, compute_matrix_H, compute_matrix_Q, compute_sufficient_statistics

def compute_likelihood(Z_prime, H, Q):
    Q_inv = np.linalg.inv(Q)
    L = 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    return L

def estimate_A_cs(Z_prime, H, Q, T_dur):
    """Ocenka A_cs cherez (99): A_cs = (2/T_dur) * Z' @ H @ Q^{-1}"""
    Q_inv = np.linalg.inv(Q)
    A_cs = Z_prime @ H @ Q_inv
    A_cs = A_cs * (2.0 / T_dur)
    return A_cs

M=3; a_m=np.array([1.0,1.0,1.0]); omega_m=2*np.pi*np.array([700, 900, 1100]); phi_m=np.array([0,0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6, 1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

# Proverim ocenku A_cs v istinnyh tau
tau = np.array([600e-6, 1200e-6])
C, S = compute_matrices_C_S(omega_m, phi_m, tau)
H = compute_matrix_H(C, S)
Q = compute_matrix_Q(omega_m, a_m, tau, N0)
A_cs = estimate_A_cs(Z_prime, H, Q, T_dur)

print("=== Ocenka A_cs v istinnyh tau (M=3) ===")
print(f"A_cs = {A_cs}")
A1_est = np.sqrt(A_cs[0]**2 + A_cs[1]**2)
psi1_est = np.arctan2(A_cs[1], A_cs[0])
A2_est = np.sqrt(A_cs[2]**2 + A_cs[3]**2)
psi2_est = np.arctan2(A_cs[3], A_cs[2])
print(f"Luch 1: A={A1_est:.4f}, psi={psi1_est:.4f} (istina: A=0.7, psi=0.52)")
print(f"Luch 2: A={A2_est:.4f}, psi={psi2_est:.4f} (istina: A=0.5, psi=-0.53)")

# Proverim H^T @ H = (N0/2) * Q
print("\n=== Proverka H^T @ H = (N0/2) * Q ===")
HtH = H.T @ H
Q_N0 = (N0/2) * Q
print(f"Max|H^T H - (N0/2)Q| = {np.max(np.abs(HtH - Q_N0)):.2e}")

# 2D scan L(T) s bolee melkoj setkoj
print("\n=== 2D poisk s setkoj 101x101 ===")
N_grid = 101
tau1_range = np.linspace(100e-6, 2000e-6, N_grid)
tau2_range = np.linspace(100e-6, 2000e-6, N_grid)
L_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            L_grid[i,j] = -np.inf
            continue
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L_grid[i,j] = compute_likelihood(Z_prime, H, Q)

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
print(f"Max L(T) pri: tau1={tau1_range[idx_max[0]]*1e6:.1f} mks, tau2={tau2_range[idx_max[1]]*1e6:.1f} mks")
print(f"Istinnye: tau1=600 mks, tau2=1200 mks")

# Ocenka A_cs v maksimume L(T)
tau_max = np.array([tau1_range[idx_max[0]], tau2_range[idx_max[1]]])
C, S = compute_matrices_C_S(omega_m, phi_m, tau_max)
H = compute_matrix_H(C, S)
Q = compute_matrix_Q(omega_m, a_m, tau_max, N0)
A_cs_max = estimate_A_cs(Z_prime, H, Q, T_dur)

print(f"\n=== Ocenka A_cs v maksimume L(T) ===")
A1_max = np.sqrt(A_cs_max[0]**2 + A_cs_max[1]**2)
psi1_max = np.arctan2(A_cs_max[1], A_cs_max[0])
A2_max = np.sqrt(A_cs_max[2]**2 + A_cs_max[3]**2)
psi2_max = np.arctan2(A_cs_max[3], A_cs_max[2])
print(f"Luch 1: A={A1_max:.4f}, psi={psi1_max:.4f} (istina: A=0.7, psi=0.52)")
print(f"Luch 2: A={A2_max:.4f}, psi={psi2_max:.4f} (istina: A=0.5, psi=-0.53)")

# Grafik
plt.figure(figsize=(15, 5))

plt.subplot(131)
plt.contourf(tau1_range*1e6, tau2_range*1e6, L_grid.T, levels=50)
plt.colorbar(label='L(T)')
plt.plot(600, 1200, 'r*', markersize=15, label='Istina')
plt.plot(tau1_range[idx_max[0]]*1e6, tau2_range[idx_max[1]]*1e6, 'wo', markersize=10, 
         label=f'Max={tau1_range[idx_max[0]]*1e6:.0f},{tau2_range[idx_max[1]]*1e6:.0f}')
plt.xlabel('tau1, mks')
plt.ylabel('tau2, mks')
plt.title('L(T1,T2) contour (M=3, N=2)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(132)
j_fixed = np.argmin(np.abs(tau2_range - 1200e-6))
L_slice = L_grid[:, j_fixed]
plt.plot(tau1_range*1e6, L_slice, 'b-', linewidth=1)
plt.axvline(600, color='g', linestyle='--', label='tau1=600')
plt.axvline(tau1_range[idx_max[0]]*1e6, color='r', linestyle=':', label=f'max={tau1_range[idx_max[0]]*1e6:.0f}')
plt.xlabel('tau1, mks')
plt.ylabel('L(T)')
plt.title('L(T1) pri tau2=1200 mks')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(133)
i_fixed = np.argmin(np.abs(tau1_range - 600e-6))
L_slice2 = L_grid[i_fixed, :]
plt.plot(tau2_range*1e6, L_slice2, 'b-', linewidth=1)
plt.axvline(1200, color='orange', linestyle='--', label='tau2=1200')
plt.axvline(tau2_range[idx_max[1]]*1e6, color='r', linestyle=':', label=f'max={tau2_range[idx_max[1]]*1e6:.0f}')
plt.xlabel('tau2, mks')
plt.ylabel('L(T)')
plt.title('L(T2) pri tau1=600 mks')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()