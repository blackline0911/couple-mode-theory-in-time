import numpy as np
import matplotlib.pyplot as plt
# Refer to imec paper : 
# Thermal Modelling of Silicon Photonic Ring Modulator with Substrate Undercut
HE = 222 # Heater efficiency (pm/mW)
thermal_optic_coeff = 1.95e-4 # dneff/dT (1/K, at 1550nm)
radius = 15
ng = 4.0815
wl = 1.55
FSR = wl**2 / (ng*2*np.pi*radius)
dT_dlambda = 1.55/(FSR*2*np.pi*radius*thermal_optic_coeff)

print("dT_dlambda = ",dT_dlambda)
print("dlambda_dT = ",1/dT_dlambda)
print("FSR = ",FSR)
print("Heater efficiency * dT/dlambda = ",HE*1e-6*dT_dlambda," oC/mW")


plt.figure()
lambd = np.linspace(1.54,1.61,1000)
ne = -1.0461*lambd+4.0815
plt.plot(lambd,ne)
plt.xlabel("wavelength")
plt.title("Dispersion in Waveguide")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Thermal induced resonant wavelength shift deducing (by PDK)
Rh = 2.99 #oC/mW thermal Resistance
HE = 254.3 #pm/mW
FSR = 0.0195 #um
radius = 5
dT_dlambda = Rh / (HE*1e-6) #1/um
dlambda_dT = 1/dT_dlambda

print("dT_dlambda = ",dT_dlambda)
print("dlambda_dT = ",dlambda_dT)
print("Rh = ",HE*1e-6*dT_dlambda)
print("dT_dlambda = ",1.55/ (FSR*2*np.pi*radius*thermal_optic_coeff))
print("dlambda_dT = ",(1.55/ (FSR*2*np.pi*radius*thermal_optic_coeff))**(-1))
print("Rh = ",HE*1e-6*1.55/ (FSR*2*np.pi*radius*thermal_optic_coeff))
# print("ng = ",( FSR*2*np.pi*radius/1.55**2)**-1)

