import numpy as np
from utility import *

dt = 0.001
t = np.arange(-15,15,dt)
L = len(t)
def func(t):
    func = np.zeros(len(t))
    for i in range(len(t)):
        if abs(t[i])<=5:
            func[i] = 1
        else:
            func[i]=0
    return func
x = func(t)
F = 1/L*np.fft.fft(x)
F  = np.fft.fftshift(F)
f = np.fft.fftfreq(t.shape[-1],d = dt)
f = np.fft.fftshift(f)
ploting(t, func(t),x_label="time",title="Rectangle")
ploting(f, F,x_label="frequency",title="FT of Rectangle")
