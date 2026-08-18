import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo import compute_matrices_C_S, compute_matrix_H, compute_matrix_Q, compute_sufficient_statistics

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
Z_prime,x_m,y_m=compute_sufficient_statistics(data['xi_observed'],data['t'],omega_m,a_m,N0,fs)

# Вычислим истинный A_cs из A_true и psi_true
A_c_true = A_true * np.cos(psi_true)  # A_ck = A_k * cos(psi_k)
A_s_true = A_true * np.sin(psi_true)  # A_sk = A_k * sin(psi_k)
A_cs_true = np.concatenate([A_c_true, A_s_true])
print("A_cs_true:", A_cs_true)

# Формула (98): L = Z' * H * A_cs^T - 0.5 * A_cs * Q * A_cs^T
# Формула (97): M = Z' * H * A_cs^T
# Формула (49): Q_val = A_cs * Q * A_cs^T

for t1 in [500e-6, 600e-6, 700e-6]:
    for t2 in [1100e-6, 1200e-6, 1300e-6]:
        tau=np.array([t1,t2])
        C,S=compute_matrices_C_S(omega_m,phi_m,tau)
        H=compute_matrix_H(C,S)
        Q=compute_matrix_Q(omega_m,a_m,tau,N0)
        
        # Формула (97): M = Z' * H * A_cs^T
        M_val = Z_prime @ H @ A_cs_true
        
        # Формула (49): Q_val = A_cs * Q * A_cs^T
        Q_val = A_cs_true @ Q @ A_cs_true
        
        # Формула (98): L = M - 0.5 * Q
        L_val = M_val - 0.5 * Q_val
        
        print(f'tau=({t1*1e6:.0f},{t2*1e6:.0f}): M={M_val:.6f}, Q={Q_val:.6f}, L={L_val:.6f}')
print("Done")