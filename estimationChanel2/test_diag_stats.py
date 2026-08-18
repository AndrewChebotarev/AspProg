"""
Diagnostika: proverka dostochnyh statistik
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
xi=data['xi_observed'].copy()
t=data['t']
s_useful = data['s_useful']

dt = 1/fs

# Schitaem integraly vruchnuyu
print("=== Ruchnoe vychislenie integralov ===")
for m in range(M):
    cos_int = np.sum(xi * np.cos(omega_m[m] * t)) * dt
    sin_int = np.sum(xi * np.sin(omega_m[m] * t)) * dt
    print(f"m={m}: int(xi*cos) = {cos_int:.10f}, int(xi*sin) = {sin_int:.10f}")
    print(f"     x_m = {(2/N0)*cos_int:.10f}, y_m = {(2/N0)*sin_int:.10f}")

# Sravnenie s poleznym signalom
print("\n=== Integraly ot poleznogo signala (bez shuma) ===")
for m in range(M):
    cos_int = np.sum(s_useful * np.cos(omega_m[m] * t)) * dt
    sin_int = np.sum(s_useful * np.sin(omega_m[m] * t)) * dt
    print(f"m={m}: int(s*cos) = {cos_int:.10f}, int(s*sin) = {sin_int:.10f}")
    print(f"     x_m = {(2/N0)*cos_int:.10f}, y_m = {(2/N0)*sin_int:.10f}")

# Teoreticheskie znacheniya
print("\n=== Teoreticheskie znacheniya (formula 87-88) ===")
print("x_m = (2/N0) * sum_k A_k * a_m * (T_dur/2) * cos(omega_m*tau_k + phi_m + psi_k)")
print("y_m = (2/N0) * sum_k A_k * a_m * (T_dur/2) * sin(omega_m*tau_k + phi_m + psi_k)")

for m in range(M):
    x_theor = 0
    y_theor = 0
    for k in range(N_true):
        phase = omega_m[m]*tau_true[k] + phi_m[m] + psi_true[k]
        x_theor += A_true[k] * a_m[m] * (T_dur/2) * np.cos(phase)
        y_theor += A_true[k] * a_m[m] * (T_dur/2) * np.sin(phase)
    x_theor *= (2/N0)
    y_theor *= (2/N0)
    print(f"m={m}: x_theor = {x_theor:.10f}, y_theor = {y_theor:.10f}")

# Teoreticheskie Z_prime
print("\n=== Teoreticheskiy Z_prime ===")
Z_theor = np.zeros(2*M)
for m in range(M):
    x_theor = 0
    y_theor = 0
    for k in range(N_true):
        phase = omega_m[m]*tau_true[k] + phi_m[m] + psi_true[k]
        x_theor += A_true[k] * a_m[m] * (T_dur/2) * np.cos(phase)
        y_theor += A_true[k] * a_m[m] * (T_dur/2) * np.sin(phase)
    x_theor *= (2/N0)
    y_theor *= (2/N0)
    Z_theor[m] = a_m[m] * x_theor
    Z_theor[M+m] = a_m[m] * y_theor
print(f"Z_prime teor = {Z_theor}")

# Proverim: chto esli N0=1?
print("\n=== S N0=1 ===")
for m in range(M):
    cos_int = np.sum(xi * np.cos(omega_m[m] * t)) * dt
    sin_int = np.sum(xi * np.sin(omega_m[m] * t)) * dt
    x_m = 2 * cos_int  # N0=1
    y_m = 2 * sin_int
    print(f"m={m}: x_m = {x_m:.6f}, y_m = {y_m:.6f}")

# Teoreticheskiy signal
print("\n=== Sravnenie: signal s 1 luchom ===")
data1=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,1,np.array([A_true[0]]),np.array([tau_true[0]]),np.array([psi_true[0]]),fs,T1,T2,40)
xi1=data1['xi_observed'].copy()
s1 = data1['s_useful']

for m in range(M):
    cos_int = np.sum(xi1 * np.cos(omega_m[m] * t)) * dt
    sin_int = np.sum(xi1 * np.sin(omega_m[m] * t)) * dt
    print(f"m={m}: int(xi1*cos) = {cos_int:.10f}, int(xi1*sin) = {sin_int:.10f}")
    print(f"     x_m = {(2/N0)*cos_int:.10f}, y_m = {(2/N0)*sin_int:.10f}")
    
    cos_int_s = np.sum(s1 * np.cos(omega_m[m] * t)) * dt
    sin_int_s = np.sum(s1 * np.sin(omega_m[m] * t)) * dt
    print(f"     int(s1*cos) = {cos_int_s:.10f}, int(s1*sin) = {sin_int_s:.10f}")