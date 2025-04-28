import numpy as np
import matplotlib.pyplot as plt
from scipy.special import *

Rn2_sheet = 638.0 
Rp2_sheet = 1020.0

hrib = 0.2114e-4 # thickness (cm)
# t = 0.3e-4 # thickness (cm)
un = 1400 #cm^2/V/s
up = 450 #cm^2/V/s
q = 1.6e-19
La = 2*np.pi*5e-4
tpb = 0.25e-4
tnb = 0.25e-4

Nd = 1/(Rn2_sheet*hrib*q*un)
Na = 1/(Rp2_sheet*hrib*q*up)
N = 1/(44*hrib*q*un)
# Nd = 3e17
# Na = 5e17
# N = 1/(44*t*q*un)

print("N2 doping concentration = ",Nd," 1/cm^3")
print("P2 doping concentration = ",Na," 1/cm^3")
print("NPLUS doping concentration = ",N," 1/cm^3")

# cj_pdk =  [23.6e-15, 20e-15, 18.5e-15]
V = np.linspace(-2,0,1000)
cj_pdk =  3.7675e-14/( (2.5485-V)**0.5 )

ni = 1.07e10
a = 1e20 #Impurity gradient
V0 = 0.0259*np.log(Na*Nd/ni**2)
# V0 = 2.54

print("V0 = ",V0)
eps_si = 8.854e-14*3.46**2
eps_clad = 8.854e-14*1.44**2
eps_box = eps_clad
W = ( 2*eps_si*(V0-V) / (q*Na*Nd*(Nd+Na)) )**0.5 * (Na+Nd)


# Cj_per_mm = 0.1*epsilon*(t)/W
Cj = eps_si*(hrib*La)/W
C_parallel = eps_si*hrib*La/W
Cf = La*(eps_clad+eps_box)/2/np.pi*np.log(2*np.pi*hrib/W)
k = np.sqrt( W*(W+tpb)*(W+tnb+tpb)/ (W+tnb) )
k_pround = np.sqrt(1-k**2)
C_top = La*eps_clad*ellipk(k_pround) /ellipk(k)
C_bottom = C_top
Cj = C_parallel+Cf+C_top+C_bottom; 

# Cj = t*La*(q*a*epsilon**2/ (12*(V0-V)))**(1/3)

plt.plot(V,Cj,label="deduced Cj")
# plt.plot(V,Cj_per_mm,label="Cj per mm")
# plt.plot(V,cj_pdk,label="PDK data")
plt.scatter([0,-1,-2],[23.6*1e-15,20*1e-15,18.5*1e-15],label="PDK discrete data",marker='o',c='r')
plt.legend()
plt.xlabel("voltage (V)")
plt.title("Junction capacitance (F)")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

plt.plot(V,W*1e4)
plt.ylabel('um')
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

print(eps_si*hrib*La/(2*eps_si/q*(1/Na+1/Nd))**0.5)