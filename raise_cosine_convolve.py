import scipy as sp
from utility import sinc
import numpy as np
import matplotlib.pyplot as  plt
import serdespy as sd
import time

# sig = np.repeat([0., 1., 1.,0. ,1.], 1)
vbias = 0.5
vpp = 1
bit_num = 3000
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

fft_signal_ideal = np.fft.fftshift(np.fft.fft(signal_ideal))*dt/t0
ny_frequency = 1/2/dt
f_total = np.linspace(-ny_frequency,ny_frequency,len(t))
f = f_total[int(len(f_total)/2):]
plt.plot(f_total/1e9, 20*np.log10(abs(fft_signal_ideal)))
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
sd.simple_eye(signal_ctle, samples_per_symbol*eye_num, bit_num//eye_num , dt,  res=100 ,title="{}Gbps NRZ signal with CTLE".format(1/T/1e9))
plt.show()

h_pulse_ctle = sp.signal.convolve(rs, h_ctle)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# FFE and DFE equalization

FFE_pre = 2
FFE_taps = 7
FFE_post = FFE_taps - FFE_pre - 1
DFE_taps = 2

# h = sd.channel_coefficients(rs[:t.size], t*t0, samples_per_symbol,FFE_pre,FFE_post,res=100)
h = sd.channel_coefficients(h_pulse_ctle[:t.size], t*t0, samples_per_symbol,FFE_pre,FFE_post,res=100)
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
RX = sd.Receiver(signal_ctle, samples_per_symbol,max(f), voltage_level, main_cursor=main_cursor,shift=False)

signal_temp = RX.signal

RX.FFE(w_ffe,FFE_pre)
signal_ffe_temp = RX.signal
sd.simple_eye(RX.signal, samples_per_symbol*eye_num, bit_num//eye_num , dt,  res=100 ,title="{}Gbps NRZ signal with FFE".format(1/T/1e9))
plt.show()

RX.signal = signal_temp
RX.nrz_DFE(w_dfe)
sd.simple_eye(RX.signal, samples_per_symbol*eye_num, bit_num//eye_num , dt,  res=100 ,title="{}Gbps NRZ signal with DFE".format(1/T/1e9))
plt.show()

RX.signal = signal_ffe_temp
RX.nrz_DFE(w_dfe)
sd.simple_eye(RX.signal, samples_per_symbol*eye_num, bit_num//eye_num , dt,  res=100 ,title="{}Gbps NRZ signal with FFE and DFE".format(1/T/1e9))
plt.show()

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# jitter

#TX signal with jitter
signal_jitter = sd.gaussian_jitter(signal_ideal, T, bit_num, samples_per_symbol, stdev=1000e-15)

#eye diagram of TX with jitter
sd.simple_eye(signal_jitter, samples_per_symbol*eye_num, bit_num//eye_num, dt, "{}Gbps NRZ Signal with jitter".format(1/T/1e9),linewidth=0.05,res=100)
plt.show()

# #signal at receiver with no jitter
# signal_out_ideal = sp.signal.convolve(h_ctle,signal_ideal)
# sd.simple_eye(signal_out_ideal[100*samples_per_symbol:], samples_per_symbol, bit_num//eye_num, dt, "rx eye diagram no jitter",linewidth=0.05)
# plt.show()

# #signal at reciever with tx jitter
# signal_out_jitter_tx = sp.signal.convolve(h,signal_jitter)
# sd.simple_eye(signal_out_jitter_tx[100*samples_per_symbol:], samples_per_symbol, bit_num//eye_num,  dt, "rx eye diagram with tx jitter",linewidth=0.05)
# plt.show()