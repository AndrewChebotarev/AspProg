"""
Pryamaya maksimizaciya FOP po vsem 6 parametram
Formula (5): L(theta) = (2/N0)*int(xi(t)*s(t,theta)dt) - (1/N0)*int(s^2(t,theta)dt)
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
t=data['t']
xi=data['xi_observed']
dt=1/fs

def compute_signal(t, A, tau, psi):
    """Vychislenie s(t,theta) po formule (3)-(4)"""
    s = np.zeros_like(t)
    for k in range(len(A)):
        for m in range(M):
            s += A[k] * a_m[m] * np.cos(omega_m[m] * (t - tau[k]) - phi_m[m] - psi[k])
    return s

def compute_likelihood_direct(theta, t, xi, N0, dt):
    """Pryamoe vychislenie FOP"""
    tau1, tau2, A1, A2, psi1, psi2 = theta
    A = np.array([A1, A2])
    tau = np.array([tau1, tau2])
    psi = np.array([psi1, psi2])
    
    s = compute_signal(t, A, tau, psi)
    
    # Formula (5): L = (2/N0)*int(xi*s)dt - (1/N0)*int(s^2)dt
    int_xi_s = np.sum(xi * s) * dt
    int_s2 = np.sum(s**2) * dt
    
    L = (2.0/N0) * int_xi_s - (1.0/N0) * int_s2
    return L

# Proverim dlya neskolkih kombinaciy
print("=== Pryamaya FOP dlya raznyh parametrov ===")
print()

# Fiksirovannye A i psi, menyaem tau
A_fixed = np.array([0.7, 0.5])
psi_fixed = np.array([0.52, -0.53])

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        theta = [t1, t2, A_fixed[0], A_fixed[1], psi_fixed[0], psi_fixed[1]]
        L = compute_likelihood_direct(theta, t, xi, N0, dt)
        print(f"tau=({t1*1e6:.0f},{t2*1e6:.0f}), A=({A_fixed[0]},{A_fixed[1]}), psi=({psi_fixed[0]:.2f},{psi_fixed[1]:.2f}): L={L:.10f}")

# Menyaem A i psi, fiksirovannye tau
print("\n=== Fiksirovannye tau, menyaem A i psi ===")
tau_fixed = np.array([600e-6, 1200e-6])

for A1_test in [0.5, 0.7, 0.9]:
    for A2_test in [0.3, 0.5, 0.7]:
        theta = [tau_fixed[0], tau_fixed[1], A1_test, A2_test, psi_fixed[0], psi_fixed[1]]
        L = compute_likelihood_direct(theta, t, xi, N0, dt)
        print(f"A=({A1_test},{A2_test}), tau=({tau_fixed[0]*1e6:.0f},{tau_fixed[1]*1e6:.0f}): L={L:.10f}")

# Poisk maksimuma po setke tau1,tau2 s istinnymi A,psi
print("\n=== Poisk maksimuma po tau1,tau2 (istinnye A,psi) ===")
N_grid = 51
tau1_range = np.linspace(0, 2000e-6, N_grid)
tau2_range = np.linspace(0, 2000e-6, N_grid)
L_grid = np.zeros((N_grid, N_grid))
for i, t1 in enumerate(tau1_range):
    for j, t2 in enumerate(tau2_range):
        if t2 <= t1:
            L_grid[i,j] = -np.inf
            continue
        theta = [t1, t2, A_fixed[0], A_fixed[1], psi_fixed[0], psi_fixed[1]]
        L_grid[i,j] = compute_likelihood_direct(theta, t, xi, N0, dt)

idx_max = np.unravel_index(np.argmax(L_grid), L_grid.shape)
print(f"Max at: tau1={tau1_range[idx_max[0]]*1e6:.1f} us, tau2={tau2_range[idx_max[1]]*1e6:.1f} us")
print(f"L_max = {L_grid[idx_max]:.10f}")

i_true = np.argmin(np.abs(tau1_range - 600e-6))
j_true = np.argmin(np.abs(tau2_range - 1200e-6))
print(f"L(true) = {L_grid[i_true, j_true]:.10f}")