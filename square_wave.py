import numpy as np
import matplotlib.pyplot as plt
from math import pi as p

def square_wave(t,driver):
    X0=driver.v_bias  # v bias
    signal = X0
    vpp = driver.vpp
    k=1
    i=0  
    f = driver.f_drive
    T=1/f/2 # Time Period
    
    if isinstance(t,np.ndarray):
        N = len(t)
        while i<N:  #Change the value of i to change the number of harmonics added
            X=vpp*(np.power(-1, (k-1)/2))*(2/(p*k))*np.exp(1j*(k*2*p*t*f))
            signal=signal+X
            k+=2
            i+=1
   