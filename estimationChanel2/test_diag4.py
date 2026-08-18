import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo import compute_matrix_Q

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
t = data['t']
xi = data['xi_observed']
dt = 1/fs

# Опорные сигналы (формулы 29-30)
beta_c = np.zeros_like(t)
beta_s = np.zeros_like(t)
for m in range(M):
    beta_c += a_m[m] * np.cos(omega_m[m] * t - phi_m[m])
    beta_s += a_m[m] * np.sin(omega_m[m] * t - phi_m[m])

# Функция для вычисления X(tau), Y(tau) через сдвиг и интегрирование
def compute_XY(tau_val):
    t_shifted = t - tau_val
    beta_c_shift = np.interp(t, t_shifted, beta_c, left=0, right=0)
    beta_s_shift = np.interp(t, t_shifted, beta_s, left=0, right=0)
    X = (2.0/N0) * np.sum(xi * beta_c_shift) * dt
    Y = (2.0/N0) * np.sum(xi * beta_s_shift) * dt
    return X, Y

print("Computing L via method 1 (quadrature correlators)...")
for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        X1, Y1 = compute_XY(t1)
        X2, Y2 = compute_XY(t2)
        Z = np.array([X1, X2, Y1, Y2])
        tau = np.array([t1, t2])
        Q = compute_matrix_Q(omega_m, a_m, tau, N0)
        try:
            Q_inv = np.linalg.inv(Q)
            L = 0.5 * Z @ Q_inv @ Z.T
        except:
            L = -np.inf
        print(f'tau=({t1*1e6:.0f},{t2*1e6:.0f}): X1={X1:.4f}, Y1={Y1:.4f}, X2={X2:.4f}, Y2={Y2:.4f}, L={L:.6f}')
print("Done")