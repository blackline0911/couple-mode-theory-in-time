from scipy.integrate import *
import numpy as np

Cpad = 6.6e-15
Cox = 34.7e-15
Rsub = 1439.8
Lint = 104e-12
Rint = 1.1
Ccc = 5e-15
Rsi = 23e3
Csi = 7e-15
Rs = 53.9
Cj = 20e-15

t0 = 1e-12

Vpp = 1.0
Vext_w = 2*np.pi*50e9*t0

# Normalization
cj_normalize = Cj
Lint_bar = Lint/t0
Ccc_bar = Ccc/t0
Csi_bar = Csi/t0
Cj_bar = Cj/t0
cj_normalize_bar = cj_normalize/t0


def circuit_Ceqs(t, paras):
    Vj, i2, Vb,  Vpn = paras
    # Q_pround = paras

    dVj_dt = (Vpn  - Vj)/Rs/Cj_bar

    di2_dt = 1/Lint_bar * ( (Vpp/2*np.cos(Vext_w*t)-1) - i2*Rint - Vpn)

    dvpn_dt = 1/Ccc_bar* (i2 - (Vpn  - Vj)/Rs \
                        - (Vpn  - Vb)/Rsi )
    
    dVb_dt = (Vpn - Vb) /Rsi/Csi_bar

    return [dVj_dt, di2_dt, dVb_dt, dvpn_dt]
    # return [dQ_dt]
    
t = np.arange(0,3000,0.01)
sol = solve_ivp(circuit_Ceqs, t_span=[0,max(t)], 
                y0=[0,0,0,0], t_eval=t, 
                atol=1e-15,rtol = 1e-15)

import matplotlib.pyplot as plt

plt.plot(t,sol.y[0],label="Vj")
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()


plt.plot(t,sol.y[1],label="i2")
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

plt.plot(t,sol.y[2],label="Vb")
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

plt.plot(t,sol.y[3],label="Vpn")
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()
    
