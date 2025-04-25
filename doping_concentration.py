import numpy as np
import matplotlib.pyplot as plt

Rn2_sheet = 638.0 
Rp2_sheet = 1020.0

t = 0.2114e-4 # thickness (cm)
# t = 0.3e-4 # thickness (cm)
un = 1400 #cm^2/V/s
up = 450 #cm^2/V/s
q = 1.6e-19
La = 2*np.pi*5e-4*0.1

Nd = 1/(Rn2_sheet*t*q*un)
Na = 1/(Rp2_sheet*t*q*up)
N = 1/(44*t*q*un)

print("N2 doping concentration = ",Nd," 1/cm^3")
print("P2 doping concentration = ",Na," 1/cm^3")
print("NPLUS doping concentration = ",N," 1/cm^3")

# cj_pdk =  [23.6e-15, 20e-15, 18.5e-15]
V = np.linspace(-10,0,1000)
cj_pdk =  3.7675e-14/( (2.5485-V)**0.5 )

ni = 1.07e10
V0 = 0.0259*np.log(Na*Nd/ni**2)
print("V0 = ",V0)
epsilon = 8.854e-14*3.46**2
W = ( 2*epsilon*(V0-V) / (q*Na*Nd*(Nd+Na)) )**0.5 * (Na+Nd)
Cj = epsilon*(t*La)/W
Cj = epsilon*(t*La)/(2*epsilon/q*(1/Na+1/Nd))**0.5 / (V0-V)**0.5
# Cj =  epsilon*(t*La)/(2*epsilon/q*(1/Na+1/Nd))**0.5/ (2.5485-V)**0.5
plt.plot(V,Cj,label="deduced Cj")
plt.plot(V,cj_pdk,label="PDK data")
plt.legend()
plt.xlabel("voltage (V)")
plt.title("Junction capacitance (F)")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

print(epsilon*(t*La)/(2*epsilon/q*(1/Na+1/Nd))**0.5)