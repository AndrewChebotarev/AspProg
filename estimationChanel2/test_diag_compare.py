"""
Proverka: kak rabotaet MonteCarlo.py (N=1) na signale s 2 luchami
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response
from MonteCarlo import estimate_delay, compute_sufficient_statistics as mc_stats
from MonteCarlo import compute_matrices_C_S, compute_matrix_H, compute_matrix_Q
from MonteCarlo import estimate_amplitudes_and_phases, compute_likelihood

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

# Signal with 2 rays
data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']

Z_prime, x_m, y_m = mc_stats(xi, t, omega_m, a_m, N0, fs)
print(f"Z_prime (from MonteCarlo.py) = {Z_prime}")

tau_est = estimate_delay(Z_prime, omega_m, a_m, phi_m, tau_true[0], T_dur, N0)
print(f"Delay est (N=1) on 2-ray signal: tau={tau_est*1e6:.1f} us (true tau1={tau_true[0]*1e6:.0f} us)")

tau_arr = np.array([tau_est])
C, S = compute_matrices_C_S(omega_m, phi_m, tau_arr)
H = compute_matrix_H(C, S)
Q = compute_matrix_Q(omega_m, a_m, tau_arr, N0)
A_est, psi_est = estimate_amplitudes_and_phases(Z_prime, H, Q, T_dur)
print(f"Amplitude: A={A_est[0]:.3f} (true A1={A_true[0]:.1f})")
print(f"Phase: psi={psi_est[0]:.3f} (true psi1={psi_true[0]:.2f})")

# Compare: 1-ray signal
data1=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,1,np.array([A_true[0]]),np.array([tau_true[0]]),np.array([psi_true[0]]),fs,T1,T2,40)
xi1=data1['xi_observed'].copy()
Z_prime1,_,_ = mc_stats(xi1, t, omega_m, a_m, N0, fs)
tau_est1 = estimate_delay(Z_prime1, omega_m, a_m, phi_m, tau_true[0], T_dur, N0)
print(f"Delay est (N=1) on 1-ray signal: tau={tau_est1*1e6:.1f} us (true tau1={tau_true[0]*1e6:.0f} us)")

# Scan L(T) over wide range
tau_range = np.linspace(0, 2000e-6, 2001)
L_2rays = []
L_1ray = []
for tau_test in tau_range:
    tau = np.array([tau_test])
    C, S = compute_matrices_C_S(omega_m, phi_m, tau)
    H = compute_matrix_H(C, S)
    Q = compute_matrix_Q(omega_m, a_m, tau, N0)
    L_2rays.append(compute_likelihood(Z_prime, H, Q))
    L_1ray.append(compute_likelihood(Z_prime1, H, Q))

L_2rays = np.array(L_2rays)
L_1ray = np.array(L_1ray)

print(f"\nL(T) for 2-ray signal:")
print(f"  Max at tau={tau_range[np.argmax(L_2rays)]*1e6:.1f} us")
print(f"  L(600 us)={L_2rays[np.argmin(np.abs(tau_range-600e-6))]:.6f}")
print(f"  L(1200 us)={L_2rays[np.argmin(np.abs(tau_range-1200e-6))]:.6f}")

print(f"\nL(T) for 1-ray signal:")
print(f"  Max at tau={tau_range[np.argmax(L_1ray)]*1e6:.1f} us")
print(f"  L(600 us)={L_1ray[np.argmin(np.abs(tau_range-600e-6))]:.6f}")

# Plot
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 6))
plt.plot(tau_range*1e6, L_2rays, 'b-', label='2-ray signal, N=1 model')
plt.plot(tau_range*1e6, L_1ray, 'r-', label='1-ray signal, N=1 model')
plt.axvline(tau_true[0]*1e6, color='g', linestyle='--', label=f'tau1={tau_true[0]*1e6:.0f} us')
plt.axvline(tau_true[1]*1e6, color='orange', linestyle='--', label=f'tau2={tau_true[1]*1e6:.0f} us')
plt.xlabel('tau, us')
plt.ylabel('L(T)')
plt.title('L(T) for N=1 model on different signals')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()