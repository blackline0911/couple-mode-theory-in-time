# Lossy transmission line modeling
import serdespy as sdp
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#% electrical loss parameters

#All AC parameters are specified at this frequency
w_0 = 2*np.pi*10e9

#loss tangent = theta_0/(2*pi)
#theta_0 = 0.022
theta_0 = 0.01

#skin-effect scaling factor [ohms/m @ w_0]
k_r = 87

#dc resistance [ohm/m]
RDC = 0.0001

# Biuld frequency vector

#fmax=1e12 -> time step = 500fs
fmax=1e12

#frequency vector (Hz)
k = 14
f = np.linspace(0,fmax,2**k+1)

#frequency vector (rad/s)
w = f*2*np.pi

# Constants

#speed of light [m/s]
c = 2.998e8

#Vacuum permittivity [F/m]
eps_0 = 8.85*1e-12

# Transmission line parameters

#Effective relative dielectric constant
eps_r = 4.9

#Propagation velocity of the transmission line [m/s]
v0 = np.sqrt(1/eps_r)*c

#Characteristic impedance [ohm]
Z0 = 50

#Inductance [H/m]
L0 = Z0/v0

#Capacitance [F/m]
C0 = 1/(Z0*v0)

#Conductance [S/m]
G0 = 1e-12

#Resistance
RAC = (k_r*(1+1j)*np.sqrt(w/w_0))

#%% Generate frequency-dependent RLGC for the lossy transmission line
R=np.sqrt(RDC**2 + RAC**2)
L=L0*np.ones(np.size(f))
G=G0*np.ones(np.size(f))
C= C0 * (1j*w/w_0)**(-2*theta_0/np.pi)

if (f[0]==0):
   C[0] = C[1]
   
# transmission line length [m]
d = 0.1

#create transmission ABCD paramaters
tline = sdp.rlgc(R,L,G,C,d,f);

#%% source impedance
r_source = 50
source = sdp.impedance(r_source*np.ones(np.size(f)))

#termination admittance
r_term = 50
termination = sdp.admittance(np.ones(np.size(f))/r_term)

#channel is the series connection of the source, transmission line, and termination
channel = sdp.series(np.array([source,tline,termination]))

#%%
# frequency domain response
Hchannel = 1/channel[:,0,0]

#Hchannel = Hchannel/abs(Hchannel[0])

np.save("./data/Hchannel",Hchannel)

# impulse response
h,t = sdp.freq2impulse(Hchannel,f);
np.save("./data/h.npy",h)
np.save("./data/t.npy",t)

#step response
hstep = sp.signal.convolve(h,np.ones(np.shape(h)))[:np.size(h)]

#100Gbps
data_rate = 100e9

#Pam-4 Signalling
t_symbol = 2/data_rate

#time step between samples of impulse response
t_sample = 1/(2*fmax)

#number of time samples in one PAM-4 symbol
samples_per_symbol = int(t_symbol/t_sample)

#response of transmission line to one UI pulse
hpulse = sp.signal.convolve(h,np.ones(np.array([samples_per_symbol,])))[:np.size(h)]
np.save("./data/hpulse.npy",hpulse)

#%% Plots
dpi = 100
plt.figure(dpi=dpi)
plt.title('Transmission Line Frequency Response')
plt.semilogx(1e-9*f,20*np.log10(np.abs(Hchannel)))
plt.xlim([0.1, 100])
plt.ylim([-40, 2])
plt.xlabel('Frequency [GHz]')
plt.ylabel('Mag Response [dB]')
plt.grid()
#nyquist f for 100G 4-PAM is 25G
plt.axvline(x=25,color = 'grey', label = "Nyquist Frequency")
plt.show()

plt.figure(dpi=dpi)
plt.plot((t*1e9),h)
plt.title('Transmission Line Impulse Response')
plt.ylabel('Impulse Response')
plt.xlabel('Time (ns)')
plt.xlim([0, 5])
#plt.ylim([-0.01, 0.08])
plt.show()

plt.figure(dpi=dpi)
plt.plot((t*1e9),hpulse)
plt.title('Transmission Line Pulse Response')
plt.ylabel('Pulse Response')
plt.xlabel('Time (ns)')
plt.xlim([0, 5])
#plt.ylim([-0.01, 0.08])
plt.show()

# plt.figure(dpi=dpi)
# plt.plot(t*1e9,hstep)
# plt.title('Transmission Line Step Response')
# plt.ylabel('Step Response [V]')
# plt.xlabel('Time (ns)')
# plt.xlim([0, 5])
# plt.show()


#%% Eye Diagram

#generate binary data
data = sdp.prqs10(1)[:10000]

#generate Baud-Rate sampled signal from data
signal_BR = sdp.pam4_input_BR(data)

#oversampled signal
signal_ideal = np.repeat(signal_BR, samples_per_symbol)

#eye diagram of ideal signal

signal_out = sp.signal.convolve(h,signal_ideal)
                                
# plt.plot(signal_ideal)
# plt.show()
# plt.plot(data)
# plt.show()
# plt.plot(signal_BR)
# plt.show()
sdp.simple_eye(signal_out[100*samples_per_symbol:], samples_per_symbol*3, 500, t_sample, res=dpi, title="{}Gbps 4-PAM Signal".format(data_rate/1e9))
plt.show()

#%% save data for next homework assignment
np.save("./data/signal.npy",signal_out)
np.save("./data/f.npy",f)
np.save("./data/w.npy",w)


#%% Generate CTLE frequency response
wz = 5e10
wp = 1.7e11
k = wp**2/wz
# H(w) = [ (k/wp^2)*s + k*wz/wp^2 ] / [ 1/wp^2*s^2 + 2/wp^2*s + 1 ]

w, H_ctle = sp.signal.freqs([k/wp**2, k*wz/wp**2], [1/wp**2, 2/wp, 1], w)
plt.semilogx(w/2/np.pi/1e9, 20*np.log10(abs(H_ctle)))
plt.xlabel("frequency (GHz)")
plt.title("CTLE filter")
plt.grid()
plt.axvline(x=25,color = 'grey', label = "Nyquist Frequency")
plt.axvline(x=wz/(2*np.pi)*1e-9,color = 'green', label = "Zero Location")
plt.axvline(x=wp/(2*np.pi)*1e-9,color = 'blue', label = "Pole Location")
plt.legend()
plt.show()

#%% Convert CTLE frequency response to impulse response
h_ctle , t_ctle = sdp.freq2impulse(H_ctle,f) 
plt.plot(t*1e9, h_ctle)
plt.xlabel("time (ns)")
plt.title("CTLE impulse reponse")
plt.grid()
plt.show()

#%% Process equalizer
signal_ctle = sp.signal.convolve(signal_out,h_ctle)
sdp.simple_eye(signal_ctle[100*samples_per_symbol:], samples_per_symbol*3, 500, 500e-15, res=dpi ,title="{}Gbps PAM4 signal with CTLE".format(data_rate*1e-9))
plt.show()

#%% Set FFE cursor
h_pulse_ctle = sp.signal.convolve(hpulse, h_ctle)

FFE_pre = 2
FFE_taps = 7
FFE_post = FFE_taps - FFE_pre - 1
DFE_taps = 2

sdp.channel_coefficients(hpulse[:t.size], t, samples_per_symbol,FFE_pre,FFE_post,res=dpi)
plt.show()
h = sdp.channel_coefficients(h_pulse_ctle[:t.size],t,samples_per_symbol,FFE_pre,FFE_post,res=dpi)
plt.show()

channel_main = h.argmax()
print("channel_main = ",channel_main)

#main_cursor = h[channel_main]
main_cursor = 1

#generate binary data
data = sdp.prqs10(1)[:10000]

voltage_levels = np.array([-3, -1, 1, 3])

#generate Baud-Rate sampled signal from data
signal_BR = sdp.pam4_input_BR(data)

signal_rx = sp.signal.fftconvolve(h, signal_BR)[:len(signal_BR)]

signal_rx_cropped = signal_rx[channel_main:]

reference_signal = signal_BR[:1000]

w_ffe_init = np.zeros([7,])
w_dfe_init = np.zeros([2,])

w_ffe, w_dfe, v_combined_ffe, v_combined_dfe, z_combined, e_combined = \
sdp.lms_equalizer(signal_rx_cropped, 0.001, len(signal_rx_cropped), w_ffe_init, FFE_pre, w_dfe_init,  voltage_levels, reference=reference_signal[:1000])

#%%

#voltage_levels = np.array([-3,-1,1,3])

nyquist_f = 25e9

RX = sdp.Receiver(signal_ctle, samples_per_symbol, nyquist_f, voltage_levels,main_cursor=main_cursor)

#sdp.simple_eye(RX.signal, samples_per_symbol*3, 800, TX.UI/TX.samples_per_symbol, "Eye Diagram with CTLE")

RX.FFE(w_ffe, FFE_pre)

sdp.simple_eye(RX.signal[int(100.5*samples_per_symbol):], samples_per_symbol*3, 800, 500*1e-15, "Eye Diagram with CTLE and FFE",res=dpi)
plt.show()
RX.pam4_DFE(w_dfe)

sdp.simple_eye(RX.signal[int(100.5*samples_per_symbol):], samples_per_symbol*3, 800, 500*1e-15, "Eye Diagram with CTLE, FFE, and DFE",res=dpi)
plt.show()