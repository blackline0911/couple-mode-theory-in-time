import numpy as np
from utility import ploting, dB

fmax = 100
f = np.linspace(0,fmax, 10000)
w = 2*np.pi*f

def parallel(z1,z2):
    return 1/(1/z1+1/z2)

# cj = 6.6e-15
# Rs = 68.1
# Cp = 31.6e-15
# Rsi = 2289.6
# Cox = 1.3e-15
# Lp = 45.2e-12
# Rp = 13.1
cj = 34e-15
Rs = 33
Cp = 6e-15
Rsi = 776
Cox = 113e-15

Z0 = 50
ZL = parallel( 1/(1j*w*1e9*Cp), \
              parallel (Rs + 1/(1j*w*1e9*cj), \
                        Rsi + 1/(1j*w*1e9*Cox) ) \
            ) 

s11 = (ZL-Z0)/(Z0+ZL)

ploting(f, dB(abs(s11)) ,x_label="frequency(GHz)",title="s11 (dB)")
ploting(f, ZL ,x_label="frequency(GHz)",title="ZL")


