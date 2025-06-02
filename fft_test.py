import numpy as np
from utility import *

dt = 0.001
vpp=1
t = np.arange(-150,150,dt)
L = len(t)
def func(t, T,beta, shift,bit,level):
    return np.where( ((t-shift)==T/2/beta)|((t-shift)==-T/2/beta), \
                (bit)/(level-1)*vpp*np.pi/4*sinc(1/2/beta), \
                (bit)/(level-1)*vpp*sinc((t-shift)/T)* np.cos(np.pi*beta*(t-shift)/T)/ (1 - (2*beta*(t-shift)/T)**2))
# def func(t):
#     func = np.zeros(len(t))
#     for i in range(len(t)):
#         if abs(t[i])<=5:
#             func[i] = 1
#         else:
#             func[i]=0
#     return func
bias = 1
x = func(t,5,1,0,1,2)+bias
F = np.fft.fft(x)*dt
F  = np.fft.fftshift(F)
f = np.fft.fftfreq(t.shape[-1],d = dt)
f = np.fft.fftshift(f)
ploting(t, func(t,5,1,0,1,2)+bias,x_label="time",title="Rectangle")
ploting(f, (F),x_label="frequency",title="FT of Rectangle")
