import numpy as np
import matplotlib.pyplot as plt

cj = [23.6e-15, 20e-15, 18.5e-15]
# a = [0.95124, 0.95248, 0.95273, 0.95292, 0.95305]
a = [0.95248, 0.95273, 0.95292, 0.95305]
L = 2*np.pi*5*1e-4
LA = L
alpha_pdk = -1/(L)*np.log(a)
print(alpha_pdk,"\n (1/cm)")
a0 = 0.95124
delta_alpha_pdk = alpha_pdk - (-1/(L)*np.log(a0))
eps_si = 8.854e-14*3.46**2
Nd = 3.309987243044366e+17
Na = 6.441134216747775e+17
ni = 1.07e10
q = 1.6e-19
t = 0.22e-4
k=1.38e-23
T=300

def Cd_V(V):
    return 3.7675e-14/( (2.5485-V)**0.5 )/LA
def W(V):
    V0 = k*T/q*np.log(Na*Nd/ni**2)
    return ( 2*eps_si*(V0-V) / (q*Na*Nd*(Nd+Na)) )**0.5 * (Na+Nd)
    # return (Cd_V(V)/t/eps_si)**(-1)
def func(V):
    # return Cd_V(V)*(V)/( W(V)-W(0) )
    return Cd_V(V)*(V)/( - ( W(V)-W(0) ) )

# delta_alpha = A * Cd(V)
# A = alpha_pdk[0]/
print( delta_alpha_pdk[0]/func(-0.5))
print( delta_alpha_pdk[1]/func(-1))
print( delta_alpha_pdk[2]/func(-1.5))
print( delta_alpha_pdk[3]/func(-2))

V = np.array([-0.5,-1,-1.5,-2])
plt.plot(V,delta_alpha_pdk[0]/func(-0.5)* func(V),label='fit by -0.5V')
plt.plot(V,delta_alpha_pdk[1]/func(-1)* func(V),label='fit by -1V')
plt.plot(V,delta_alpha_pdk[2]/func(-1.5)* func(V),label='fit by -1.5V')
plt.plot(V,delta_alpha_pdk[3]/func(-2)* func(V),label='fit by -2V')
plt.scatter(V,delta_alpha_pdk,c='r',marker='o',label='pdk data')
plt.xlabel('Voltage')
plt.ylabel('Delta alpha (1/cm)')
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

plt.plot(V,W(V)*1e4)
plt.show()

alpha0 = (-1/(L)*np.log(a0))
m,b = np.polyfit(V, alpha_pdk, 1)
plt.scatter(V,alpha_pdk,c='r',marker='o',label='pdk data')
plt.plot(V,m*V+b)
plt.show()

me_pdk = [39.5, 37.9]
V = np.array([0,-0.25])
m,b = np.polyfit(V,me_pdk, 1)
plt.scatter(V,me_pdk,c='r',marker='o',label='pdk data')
V = np.array([0,-0.5,-1,-1.5,-2])
plt.plot(V,m*V+b)
plt.show()
