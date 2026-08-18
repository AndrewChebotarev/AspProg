"""
Posledovatelnaya ocenka s pryamoy FOP dlya N=1
"""
import numpy as np
from CreateDataForEstimation import generate_test_signal_and_channel_response

M=2; a_m=np.array([1.0,1.0]); omega_m=2*np.pi*np.array([900,1100]); phi_m=np.array([0,0])
N_true=2; A_true=np.array([0.7,0.5]); tau_true=np.array([600e-6,1200e-6]); psi_true=np.array([0.52,-0.53])
fs=20000; T1=0.0; T2=1.0; T_dur=T2-T1; N0=100.0

data=generate_test_signal_and_channel_response(M,a_m,omega_m,phi_m,N_true,A_true,tau_true,psi_true,fs,T1,T2,40)
t=data['t']
xi=data['xi_observed'].copy()
dt=1/fs

def compute_signal_1ray(t, A, tau, psi):
    s = np.zeros_like(t)
    for m in range(M):
        s += A * a_m[m] * np.cos(omega_m[m] * (t - tau) - phi_m[m] - psi)
    return s

def likelihood_1ray(xi, t, A, tau, psi, N0, dt):
    s = compute_signal_1ray(t, A, tau, psi)
    int_xi_s = np.sum(xi * s) * dt
    int_s2 = np.sum(s**2) * dt
    return (2.0/N0) * int_xi_s - (1.0/N0) * int_s2

def estimate_1ray(xi, t, tau_search_range, N0, dt):
    """Ocenka parametrov odnogo lucha cherez pryamuyu FOP"""
    best_L = -np.inf
    best_params = None
    
    for tau in tau_search_range:
        # Dlya kazhdogo tau, ocenivaem A i psi cherez skan
        for A_test in np.linspace(0.1, 1.5, 15):
            for psi_test in np.linspace(-np.pi, np.pi, 19):
                L = likelihood_1ray(xi, t, A_test, tau, psi_test, N0, dt)
                if L > best_L:
                    best_L = L
                    best_params = (tau, A_test, psi_test)
    
    return best_params, best_L

# Shag 1: ocenka 1-go lucha
print("Shag 1: ocenka 1-go lucha...")
tau_search = np.linspace(400e-6, 800e-6, 41)
(tau1, A1, psi1), L1 = estimate_1ray(xi, t, tau_search, N0, dt)
print(f"  tau1={tau1*1e6:.1f} us (true={tau_true[0]*1e6:.0f}), A1={A1:.3f} (true={A_true[0]:.1f}), psi1={psi1:.3f} (true={psi_true[0]:.2f})")

# Shag 2: vychitaem 1-y luch
s1 = compute_signal_1ray(t, A1, tau1, psi1)
xi_res = xi - s1

# Shag 3: ocenka 2-go lucha
print("Shag 2: ocenka 2-go lucha...")
tau_search2 = np.linspace(1000e-6, 1400e-6, 41)
(tau2, A2, psi2), L2 = estimate_1ray(xi_res, t, tau_search2, N0, dt)
print(f"  tau2={tau2*1e6:.1f} us (true={tau_true[1]*1e6:.0f}), A2={A2:.3f} (true={A_true[1]:.1f}), psi2={psi2:.3f} (true={psi_true[1]:.2f})")

print(f"\nOshibki:")
print(f"  tau: {(tau1-tau_true[0])*1e6:.2f}, {(tau2-tau_true[1])*1e6:.2f} us")
print(f"  A:   {(A1-A_true[0])/A_true[0]*100:.2f}, {(A2-A_true[1])/A_true[1]*100:.2f} %")
print(f"  psi: {psi1-psi_true[0]:.4f}, {psi2-psi_true[1]:.4f} rad")

# Bolee tochnyy poisk
print("\n=== Bolee tochnyy poisk ===")
tau_search_fine = np.linspace(tau1 - 50e-6, tau1 + 50e-6, 51)
(tau1f, A1f, psi1f), L1f = estimate_1ray(xi, t, tau_search_fine, N0, dt)
print(f"Luch 1: tau={tau1f*1e6:.2f} us, A={A1f:.4f}, psi={psi1f:.4f}")

s1f = compute_signal_1ray(t, A1f, tau1f, psi1f)
xi_resf = xi - s1f

tau_search_fine2 = np.linspace(tau2 - 50e-6, tau2 + 50e-6, 51)
(tau2f, A2f, psi2f), L2f = estimate_1ray(xi_resf, t, tau_search_fine2, N0, dt)
print(f"Luch 2: tau={tau2f*1e6:.2f} us, A={A2f:.4f}, psi={psi2f:.4f}")

print(f"\nOshibki (tochnyy poisk):")
print(f"  tau: {(tau1f-tau_true[0])*1e6:.2f}, {(tau2f-tau_true[1])*1e6:.2f} us")
print(f"  A:   {(A1f-A_true[0])/A_true[0]*100:.2f}, {(A2f-A_true[1])/A_true[1]*100:.2f} %")
print(f"  psi: {psi1f-psi_true[0]:.4f}, {psi2f-psi_true[1]:.4f} rad")