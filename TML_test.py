from scipy.integrate import *
import numpy as np
from utility import *

def parallel(Z1,Z2):
    return 1/(1/Z1+1/Z2)

Cj = 20e-15
Rs = 59.3
Rsi = 1439.8
Cox = 34.7e-15
Cp = 6.6e-15
f0 = 1e12

Cp_bar = Cp*f0
Cox_bar = Cox*f0
Cj_bar = Cj*f0

f_bar = np.arange(1e3,100e9,1e6)/f0
w_bar = 2*np.pi*f_bar
Z0 = 50 # TML characteristic impedence
ZL = parallel( parallel(1/(1j*w_bar*Cp_bar), (Rsi+1/(1j*w_bar*Cox_bar) ) ) , (Rs+1/(1j*w_bar*Cj_bar) ))

s11 = (ZL-Z0)/(Z0+ZL)
ploting(f_bar*f0/1e9, dB(abs(s11)), x_label='Frequency (GHz)',title="|s11|")
