"""
Test: posledovatelnaya ocenka s novym algoritmom (pryamaya FOP + ispravlennaya Q)
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo2 import *

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; dt=1/fs; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

# Dostatochnye statistiki
Z_prime,x_m,y_m=compute_sufficient_statistics(xi,t,omega_m,a_m,N0,fs)

# 2D poisk zaderzhek
print("2D poisk zaderzhek...")
tau1_est, tau2_est, L_max = estimate_delays_2d(xi, t, Z_prime, omega_m, a_m, phi_m, tau_true, N0, dt, T_dur, n_grid=31)
print(f"tau1_est={tau1_est*1e6:.1f} mks (true={tau_true[0]*1e6:.0f})")
print(f"tau2_est={tau2_est*1e6:.1f} mks (true={tau_true[1]*1e6:.0f})")
print(f"L_max={L_max:.10f}")

# Ocenka A i psi
tau_est = np.array([tau1_est, tau2_est])
_, A_est, psi_est = estimate_parameters_for_tau(xi, t, Z_prime, omega_m, a_m, phi_m, tau_est, N0, dt, T_dur)
print(f"\nA_est={A_est} (true={A_true})")
print(f"psi_est={psi_est} (true={psi_true})")

print(f"\nOshibki:")
print(f"  tau: {(tau1_est-tau_true[0])*1e6:.2f}, {(tau2_est-tau_true[1])*1e6:.2f} mks")
print(f"  A:   {(A_est[0]-A_true[0])/A_true[0]*100:.2f}, {(A_est[1]-A_true[1])/A_true[1]*100:.2f} %")
print(f"  psi: {psi_est[0]-psi_true[0]:.4f}, {psi_est[1]-psi_true[1]:.4f} rad")