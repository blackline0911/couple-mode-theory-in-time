import scipy as sp
from utility import sinc
import numpy as np
import matplotlib.pyplot as  plt
import serdespy as sd
import time

# sig = np.repeat([0., 1., 1.,0. ,1.], 1)
bit_num = 1000
dt = 1e-3
T = 1.0
tmax = T*(bit_num)
sig = sd.prbs24(1)[:bit_num ]
t = np.arange(0,tmax,dt)
beta = 1.0

pi = np.pi
def rcos(t, T,beta, shift):
    return np.where( ((t-shift)==T/2/beta)|((t-shift)==-T/2/beta), \
               pi/4/T*sinc(1/2/beta), \
                1/T*sinc((t-shift)/T)* np.cos(pi*beta*(t-shift)/T)/ (1 - (2*beta*(t-shift)/T)**2))
rs = rcos(t, T,beta,tmax/2)

plt.plot(sig)
plt.show()
plt.plot(rs)
plt.show()

t1 = time.time()
signal_ideal = np.zeros(len(t))
for i in range(bit_num):
    signal_ideal +=  rcos(t, T,beta,i*T) * sig[i]
t2 = time.time()

print("time of for loop = ",t2-t1)
plt.plot(t,signal_ideal)
plt.show()
samples_per_symbol = int(T/dt)
print(samples_per_symbol)
print(len(signal_ideal))
eye_num = 3
sd.simple_eye(signal_ideal, samples_per_symbol*eye_num, bit_num//eye_num , dt*1e-9, res=100, title="ideal {}Gbps 4-PAM Signal".format(T))
plt.show()


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# CTLE
f = np.fft.fftshift(np.fft.fftfreq(len(t),dt*1e-9))
w = 2*pi*f
wz = 5e10
wp = 1.7e11
k = wp**2/wz
# H(w) = [ (k/wp^2)*s + k*wz/wp^2 ] / [ 1/wp^2*s^2 + 2/wp^2*s + 1 ]

w, H_ctle = sp.signal.freqs([k/wp**2, k*wz/wp**2], [1/wp**2, 2/wp, 1], w)
plt.plot(f)
plt.show()
plt.semilogx(w/2/np.pi/1e9, 20*np.log10(abs(H_ctle)))
plt.xlabel("frequency (GHz)")
plt.title("CTLE filter")
plt.grid()
plt.axvline(x=max(f*1e-9),color = 'grey', label = "Nyquist Frequency")
plt.axvline(x=wz/(2*np.pi)*1e-9,color = 'green', label = "Zero Location")
plt.axvline(x=wp/(2*np.pi)*1e-9,color = 'blue', label = "Pole Location")
plt.legend()
plt.show()

h_ctle , t_ctle = sd.freq2impulse(H_ctle,f) 
signal_ctle = sp.signal.convolve(signal_ideal,h_ctle)
sd.simple_eye(signal_ctle, samples_per_symbol*eye_num, bit_num//eye_num , dt*1e-9,  res=100 ,title="{}Gbps PAM4 signal with CTLE".format(T*1e-9))
plt.show()

h_pulse_ctle = sp.signal.convolve(rs, h_ctle)