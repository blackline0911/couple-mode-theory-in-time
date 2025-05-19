import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

cj = [23.6e-15, 20e-15, 18.5e-15]
a = [0.95124, 0.95248, 0.95273, 0.95292, 0.95305]
# a = [0.95248, 0.95273, 0.95292, 0.95305]
L = 2*np.pi*5*1e-4
LA = L
alpha_pdk = -2/(L)*np.log(a)
print(alpha_pdk,"\n (1/cm)")
a0 = 0.95124
delta_alpha_pdk = alpha_pdk - (-2/(L)*np.log(a0))
eps_si = 8.854e-14*3.46**2
Nd = 3.309987243044366e+17
Na = 6.441134216747775e+17
ni = 1.07e10
q = 1.6e-19
t = 0.22e-4
k=1.38e-23
T=300

# def Cd_V(V):
#     return 3.7675e-14/( (2.5485-V)**0.5 )/LA
# def W(V):
#     V0 = k*T/q*np.log(Na*Nd/ni**2)
#     return ( 2*eps_si*(V0-V) / (q*Na*Nd*(Nd+Na)) )**0.5 * (Na+Nd)
#     # return (Cd_V(V)/t/eps_si)**(-1)
# def func(V):
#     # return Cd_V(V)*(V)/( W(V)-W(0) )
#     return Cd_V(V)*(V)/( - ( W(V)-W(0) ) )

# # delta_alpha = A * Cd(V)
# # A = alpha_pdk[0]/
# print( delta_alpha_pdk[0]/func(-0.5))
# print( delta_alpha_pdk[1]/func(-1))
# print( delta_alpha_pdk[2]/func(-1.5))
# print( delta_alpha_pdk[3]/func(-2))

# V = np.array([0,-0.5,-1,-1.5,-2])
# plt.plot(V,delta_alpha_pdk[0]/func(-0.5)* func(V),label='fit by -0.5V')
# plt.plot(V,delta_alpha_pdk[1]/func(-1)* func(V),label='fit by -1V')
# plt.plot(V,delta_alpha_pdk[2]/func(-1.5)* func(V),label='fit by -1.5V')
# plt.plot(V,delta_alpha_pdk[3]/func(-2)* func(V),label='fit by -2V')
# plt.scatter(V,delta_alpha_pdk,c='r',marker='o',label='pdk data')
# plt.xlabel('Voltage')
# plt.ylabel('Delta alpha (1/cm)')
# plt.legend()
# plt.grid(color='g',linestyle='--', alpha=0.5)
# plt.show()

# plt.plot(V,W(V)*1e4)
# plt.show()

V = np.array([0,-0.5,-1,-1.5,-2])
alpha0 = -2/(L)*np.log(a0)
def func(v, a, b,c):
    return a*v/(abs(v)+b)**0.5 + c
# para = np.polyfit(V,alpha_pdk, 1)
popt, pcov = curve_fit(func, V, alpha_pdk)
print("popt = ",popt)
print("alpha_pdk = ",alpha_pdk)
plt.scatter(np.array([0,-0.5,-1,-1.5,-2]),alpha_pdk,c='r',marker='o',label='pdk data')
plt.grid(color='g',linestyle='--', alpha=0.5)
v = np.linspace(-3,1,1000)
plt.xlabel("voltage (V)")
plt.title("absorption coefficient (1/cm)")
plt.plot(v,func(v,popt[0],popt[1],popt[2]),label='fitting curve')
plt.legend()
plt.show()


dneff_dV = 7.2244e-05
neff0 = 2.69021788
neff_pdk = [neff0+dneff_dV*(-0.5),neff0, neff0+dneff_dV*0.5, neff0+dneff_dV*1, neff0+dneff_dV*1.5, neff0+dneff_dV*2]
# neff_pdk = [2.51105, 2.5111, 2.51113, 2.51116, 2.51118, 2.5112]
V = np.array([0.5,0,-0.5,-1,-1.5,-2])
# a,b,c = np.polyfit(V,neff_pdk, 2)
# popt, pcov = curve_fit(func, V, neff_pdk)
plt.scatter(V,neff_pdk,c='r',marker='o',label='pdk data')
V = np.linspace(-3,0.5,1000)
# plt.plot(V,m*V+b,label="Linea fitting") 
# plt.plot(V,func(V,popt[0],popt[1],popt[2]),label="Quadra fitting") 
plt.xlabel("voltage")
plt.title("Neff vs Voltage")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.legend()
plt.show()

gamma_pdk = [0.947445, 0.947582, 0.947633, 0.947595, 0.947592, 0.947325]
V = np.array([0.5,0,-0.5,-1,-1.5,-2])
def gamma_func(v, a, b,c):
    # return np.sin((a*v/(abs(v)+b)**0.5 + c)*L)
    return c*np.sin((a*v)-b)
# a,b,c = np.polyfit(V,neff_pdk, 2)
popt, pcov = curve_fit(gamma_func, V, gamma_pdk)
print(popt)
plt.scatter(V,gamma_pdk,c='r',marker='o',label='pdk data')
V = np.linspace(-3,0.5,1000)
# plt.plot(V,m*V+b,label="Linea fitting") 
plt.plot(V,gamma_func(V,popt[0],popt[1],popt[2]),label="func fitting") 
plt.xlabel("voltage")
plt.title("gamma vs Voltage")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.legend()
plt.show()