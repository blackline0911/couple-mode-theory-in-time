import scipy as sp
from utility import sinc
import numpy as np
import matplotlib.pyplot as  plt
import serdespy as sd
import time

# sig = np.repeat([0., 1., 1.,0. ,1.], 1)
vbias = 0.5
vpp = 1
bit_num = 1000
dt = 1e-13
t0 = 1e-12
f = 100e9
T = 1/f
tmax = T/t0*(bit_num)
sig = sd.prbs24(1)[:bit_num ]
t = np.arange(0,tmax,dt/t0)
beta = 1.0

pi = np.pi
def rcos(t, T,beta, shift,bit,level):
    return np.where( ((t-shift)==T/2/beta)|((t-shift)==-T/2/beta), \
                (bit)/(level-1)*vpp*pi/4*sinc(1/2/beta), \
                (bit)/(level-1)*vpp*sinc((t-shift)/T)* np.cos(pi*beta*(t-shift)/T)/ (1 - (2*beta*(t-shift)/T)**2)) 
rs = rcos(t, T/t0,beta,tmax/2,bit = 1,level = 2)

t1 = time.time()
signal_ideal = np.zeros(len(t))
for i in range(bit_num):
    signal_ideal +=  rcos(t, T/t0,beta,shift=i*T/t0,bit=sig[i],level = 2) 
signal_ideal += (vbias-vpp/2)
t2 = time.time()

print("time of for loop = ",t2-t1)
samples_per_symbol = round(T/dt)
eye_num = 3
sd.simple_eye(signal_ideal, samples_per_symbol*eye_num, bit_num//eye_num , dt, res=100, title="ideal {}Gbps 4-PAM Signal".format(1/T*1e-9))
plt.show()
plt.plot(signal_ideal)
plt.show()

fft_signal_ideal = np.fft.fftshift(np.fft.fft(signal_ideal))*dt/t0
f = np.fft.fftshift(np.fft.fftfreq(len(t),d=dt))
ny_frequency = 1/2/dt
f = np.linspace(0,ny_frequency,len(t)/2)
plt.plot(f/1e9, 20*np.log10(abs(fft_signal_ideal[:len(t)])))
plt.xlabel("frequency (GHz)")
plt.title("Frequency response of raise cosine PRBS")
plt.grid()
plt.show()

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# CTLE

w = 2*pi*f
wz = 5e10
wp = 1.7e11
k = wp**2/wz
# H(w) = [ (k/wp^2)*s + k*wz/wp^2 ] / [ 1/wp^2*s^2 + 2/wp^2*s + 1 ]

w, H_ctle = sp.signal.freqs([k/wp**2, k*wz/wp**2], [1/wp**2, 2/wp, 1], w)
# plt.plot(f[1:]-f[0:-1])
# plt.show()
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
signal_ctle = sp.signal.fftconvolve(signal_ideal,h_ctle,mode="same")
sd.simple_eye(signal_ctle, samples_per_symbol*eye_num, bit_num//eye_num , dt,  res=100 ,title="{}Gbps PAM4 signal with CTLE".format(1/T/1e9))
plt.show()

h_pulse_ctle = sp.signal.convolve(rs, h_ctle)

plt.plot(t,signal_ctle)
plt.title("signal_ctle")
plt.show()
plt.plot(h_ctle)
plt.title("h_ctle")
plt.show()


fft_signal_ideal = fft_signal_ideal[int(len(f)/2):]
fa = H_ctle*fft_signal_ideal
Hd = np.concatenate(  fa, np.conj(np.flip(fa[1:fa.size-1]) )  )
ifft = np.real(np.fft.ifft(Hd))
# plt.plot(ifft)
plt.title("impulse response by self do ifft")
# plt.show()

sd.simple_eye(ifft, samples_per_symbol*eye_num, bit_num//eye_num , dt,  res=100 ,title="impulse response by do ifft")
plt.show()
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# FFE and DFE equalization

FFE_pre = 2
FFE_taps = 7
FFE_post = FFE_taps - FFE_pre - 1
DFE_taps = 2

# h = sd.channel_coefficients(rs[:t.size], t*t0, samples_per_symbol,FFE_pre,FFE_post,res=100)
h = sd.channel_coefficients(h_pulse_ctle[:t.size], t*t0, samples_per_symbol,FFE_pre,FFE_post,res=100)
print("h = ",h)
plt.show()  
channel_main = h.argmax()

main_cursor = 1

voltage_level = np.array([vbias-vpp/2, vbias+vpp/2])
signal_BR = sd.nrz_input_BR(sig,voltage_levels=voltage_level)

signal_rx = sp.signal.fftconvolve(h, signal_BR)[:len(signal_BR)]

signal_rx_cropped = signal_rx[channel_main:]

reference_signal = signal_BR[:1000]

w_ffe_init = np.zeros([FFE_taps,])
w_dfe_init = np.zeros([FFE_pre,])

w_ffe, w_dfe, v_combined_ffe, v_combined_dfe, z_combined, e_combined = \
sd.lms_equalizer(signal_rx_cropped, 1e-3, len(signal_rx_cropped), w_ffe_init, FFE_pre, w_dfe_init, voltage_level, reference=reference_signal[:1000])

# RX = sd.Receiver(signal_ideal, samples_per_symbol,max(f), voltage_level, main_cursor=main_cursor)
RX = sd.Receiver(signal_ctle, samples_per_symbol,max(f), voltage_level, main_cursor=main_cursor)

RX.FFE(w_ffe,FFE_pre)

sd.simple_eye(RX.signal, samples_per_symbol*eye_num, bit_num//eye_num , dt,  res=100 ,title="{}Gbps PAM4 signal with CTLE".format(1/T/1e9))
plt.show()