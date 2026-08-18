"""
Proverka L(T) dlya M=3, N=2 s originalnymi parametrami (tau1=600, tau2=1200 mks)
"""
import numpy as np
import matplotlib.pyplot as plt
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import compute_matrices_C_S, compute_matrix_H, compute_matrix_Q, compute_sufficient_statistics

def compute_likelihood(Z_prime, H, Q):
    """L(T) = 0.5 * Z' @ H @ Q^{-1} @ H^T @ Z'^T (formula 100)"""
    Q_inv = np.linalg.inv(Q)
    L = 0.5 * Z_prime @ H @ Q_inv @ H.T @ Z_prime.T
    return L

# M=3, N=2, originalnye parametry
M=3; a_m=np.array([1.0,1.0,1.0]); omega_m=2*np.pi*np.array([700, 900, 1100]); phi_m=np.array([0,0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6, 1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)
print(f"Z_prime (M=3) = {Z_prime}")
print(f"Razmer Z_prime: {Z_prime.shape}")

# Proverim v tochnyh tochkah
print("\n=== L(T) v tochnyh tochkah ===")
for t1 in [500e-6, 550e-6, 600e-6, 650e-6, 700e-6]:
    for t2 in [1100e-6, 1150e-6, 1200e-6, 1250e-6, 1300e-6]:
        if t2 <= t1:
            continue
        tau = np.array([t1, t2])
        C, S = compute_matrices_C_S(omega_m, phi_m, tau)
        H = compute_matrix_H(C, S)
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        L = compute_likelihood(Z_prime, H, Q)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}): L={L:.10f}")

# 2D scan po setke
print("\n=== 2D poisk maksimuma L(T) ===")
N_grid = 51
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
print(f"L_max = {L_grid[idx_max]:.10f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"Istinnye: tau1=600 mks, tau2=1200 mks")
print(f"L(true) = {L_grid[i_true, j_true]:.10f}")

# Sravnenie: est li maksimum v istinnyh tochkah?
L_true = L_grid[i_true, j_true]
L_max = L_grid[idx_max]
print(f"\nL_max / L_true = {L_max / L_true:.6f}")
print(f"Raznica L_max - L_true = {L_max - L_true:.10f}")

# Grafik L(T) v razreze pri fiksirovannom tau2=1200 mks
print("\n=== L(T1) pri tau2=1200 mks ===")
j_fixed = np.argmin(np.abs(tau2_range - 1200e-6))
L_slice = L_grid[:, j_fixed]
idx_slice = np.argmax(L_slice)
print(f"Max L(T1) pri tau2=1200: tau1={tau1_range[idx_slice]*1e6:.1f} mks, L={L_slice[idx_slice]:.10f}")

plt.figure(figsize=(12, 5))
plt.subplot(121)
plt.plot(tau1_range*1e6, L_slice, 'b-', linewidth=1)
plt.axvline(600, color='g', linestyle='--', label='tau1=600 mks')
plt.axvline(tau1_range[idx_slice]*1e6, color='r', linestyle=':', label=f'max={tau1_range[idx_slice]*1e6:.1f}')
plt.xlabel('tau1, mks')
plt.ylabel('L(T)')
plt.title(f'L(T1) pri tau2=1200 mks (M=3)')
plt.legend()
plt.grid(True, alpha=0.3)

# Grafik L(T) v razreze pri fiksirovannom tau1=600 mks
print(f"\n=== L(T2) pri tau1=600 mks ===")
i_fixed = np.argmin(np.abs(tau1_range - 600e-6))
L_slice2 = L_grid[i_fixed, :]
idx_slice2 = np.argmax(L_slice2)
print(f"Max L(T2) pri tau1=600: tau2={tau2_range[idx_slice2]*1e6:.1f} mks, L={L_slice2[idx_slice2]:.10f}")

plt.subplot(122)
plt.plot(tau2_range*1e6, L_slice2, 'b-', linewidth=1)
plt.axvline(1200, color='orange', linestyle='--', label='tau2=1200 mks')
plt.axvline(tau2_range[idx_slice2]*1e6, color='r', linestyle=':', label=f'max={tau2_range[idx_slice2]*1e6:.1f}')
plt.xlabel('tau2, mks')
plt.ylabel('L(T)')
plt.title(f'L(T2) pri tau1=600 mks (M=3)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# 2D contour plot
plt.figure(figsize=(10, 8))
plt.contourf(tau1_range*1e6, tau2_range*1e6, L_grid.T, levels=50)
plt.colorbar(label='L(T)')
plt.plot(600, 1200, 'r*', markersize=15, label='Istina')
plt.plot(tau1_range[idx_max[0]]*1e6, tau2_range[idx_max[1]]*1e6, 'wo', markersize=10, label=f'Max={tau1_range[idx_max[0]]*1e6:.0f},{tau2_range[idx_max[1]]*1e6:.0f}')
plt.xlabel('tau1, mks')
plt.ylabel('tau2, mks')
plt.title('L(T1,T2) contour (M=3, N=2)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()