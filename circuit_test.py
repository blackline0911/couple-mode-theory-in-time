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

Vext_amp = 1.0
Vext_w = 2*np.pi*1e3

def circuit_Ceqs(t, paras):
    Q, i2, i4, Vpn = paras

    dQ_dt = (Vpn - Q/Cj)/Rs

    di2_dt = 1/Lint * (-Vext_amp*np.cos(Vext_w*t) + i2*Rint + Vpn)

    dvpn_dt = 1/Ccc* (i2 - dQ_dt - i4)
    
    di4_dt = -i4 + Csi/Rsi*dvpn_dt

    return [dQ_dt, di2_dt, di4_dt, dvpn_dt]
    
t = np.arange(0,0.01,0.000001)
sol = solve_ivp(circuit_Ceqs, t_span=[0,0.01], y0=[0, 0, 0, 0], t_eval=t, atol=1e-20,rtol = 1e-15)
import matplotlib.pyplot as plt
print(sol.y[1])
plt.plot(t,sol.y[0],label="Q")
plt.plot(t,sol.y[1],label="i2")
plt.plot(t,sol.y[2],label="i4")
plt.plot(t,sol.y[3],label="Vpn")
plt.legend()
plt.show()
    
