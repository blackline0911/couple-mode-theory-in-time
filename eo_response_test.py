from utility import *
import os
os.chdir("./eye_diagram_test_v4/")
signal = np.load("eye_SmallSignal_1000bits_with_TPA_FCA_f44GHz_vpp_-3000mV_vbias_3000mVNRZ.npy")
# vpos = np.load("voltage_eye_44GHz_vpp_3000mV_vbias_-3000mV.npy")
vpos_20ghz = np.load("voltage_eye_20GHz_vpp_3000mV_vbias_-3000mV.npy")
vpos_neg_20ghz = np.load("vpos + vneg SmallSignalNRZ.npy")
vj_20ghz = np.load("vj_SmallSignalNRZ.npy")
# vneg = np.load("vneg_SmallSignalNRZ.npy")

# F_vpos,f = FFT(vpos, 1e-14)
F_vpos_20ghz,f = FFT(vpos_20ghz, 1e-14)
F_vpos_neg_20ghz,f = FFT(vpos_neg_20ghz, 1e-14)
F_vj_20ghz,f = FFT(vj_20ghz, 1e-14)

# F_vneg,f = FFT(vneg, 1e-14)
m = 10000
# ploting(f[:m]/1e9, dB(abs(F_vpos_20ghz[:m])),
#         x_label="Frequency (GHz)", title="Voltage Spectrum (Positive 20GHz)",
#         filename="voltage_spectrum_pos_20GHz")
ploting(f[:m]/1e9, dB(abs(F_vpos_neg_20ghz[:m])),
        x_label="Frequency (GHz)", title="Vpos+Vneg Spectrum (20GHz)",
        filename="voltage_spectrum_pos_neg_20GHz")
ploting(f[:m]/1e9, dB(abs(F_vj_20ghz[:m])),dB(abs(F_vpos_neg_20ghz[:m])),
        x_label="Frequency (GHz)", title="Vj Spectrum (20GHz)",
        filename="vj_spectrum_20GHz")
ploting(f[:m]/1e9, dB(abs(F_vj_20ghz[:m]/F_vpos_neg_20ghz[:m])),
        x_label="Frequency (GHz)", title="V Spectrum (20GHz)",
        filename="v_spectrum_20GHz")
# ploting(f[:m]/1e9, dB(abs(F_vpos[:m])),
#         x_label="Frequency (GHz)", title="Voltage Spectrum (Positive, 44GHz)",
#         filename="voltage_spectrum_pos")
# ploting(f[:m]/1e9, dB(abs(F_vneg[:m])),
#         x_label="Frequency (GHz)", title="Voltage Spectrum (reflective)",
#         filename="voltage_spectrum_neg")
# ploting(f[:m]/1e9, dB(abs(F_vneg[:m]/F_vpos[:m])),
#         x_label="Frequency (GHz)", title="s11 response",
#         filename="s11_response")




# optical_power_spectrum, frequency_axis = FFT(signal , 1e-14)
# ploting(frequency_axis[1:3000]/1e9, dB(optical_power_spectrum[1:3000]),
#         x_label="Frequency (GHz)", title="Optical Power Spectrum",
#         filename="optical_power_spectrum")