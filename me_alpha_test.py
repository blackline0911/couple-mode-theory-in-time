import numpy as np
import matplotlib.pyplot as plt

cj = [23.6e-15, 20e-15, 18.5e-15]
a = [0.95124, 0.95248, 0.95273, 0.95292, 0.95305]
# a = [0.95248, 0.95273, 0.95292, 0.95305]
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
alpha0 = -1/(L)*np.log(a0)
# m,b = np.polyfit(V,alpha_pdk-alpha0, 1)
# a,b,c = np.polyfit(V,alpha_pdk, 2)
para = np.polyfit(V,alpha_pdk, 1)
plt.scatter(np.array([0,-0.5,-1,-1.5,-2]),alpha_pdk,c='r',marker='o',label='pdk data')
plt.grid(color='g',linestyle='--', alpha=0.5)
# plt.plot(np.linspace(-2,0,1000),a*(np.linspace(-2,0,1000))**2+b*np.linspace(-2,0,1000)+c)
v = np.linspace(-3,0,1000)
# plt.plot(v,m*v+b+alpha0)
# plt.plot(v,a*v**2+b*v+c)
A = 0
for i in range(len(para)):
    A+=para[i]*v**(len(para)-1-i)
plt.plot(v,A,label='fitting curve')
plt.legend()
plt.show()

# me_pdk = [39.5, 37.9]
# V = np.array([0,-0.25])
# m,b = np.polyfit(V,me_pdk, 1)
# plt.scatter(V,me_pdk,c='r',marker='o',label='pdk data')
# V = np.array([0,-0.5,-1,-1.5,-2])
# plt.plot(V,m*V+b)
# plt.show()

neff_pdk = [2.51105, 2.5111, 2.51113, 2.51116, 2.51118, 2.5112]
V = np.array([0.5,0,-0.5,-1,-1.5,-2])
a,b,c = np.polyfit(V,neff_pdk, 2)
plt.scatter(V,neff_pdk,c='r',marker='o',label='pdk data')
V = np.linspace(-2,0.5,1000)
# plt.plot(V,m*V+b,label="Linea fitting") 
plt.plot(V,a*V**2+b*V+c,label="Quadra fitting") 
plt.xlabel("voltage")
plt.title("Neff vs Voltage")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.legend()
plt.show()

plt.plot(V,1.5486e6/3.905*(m*V+b-b),label="me")  
plt.xlabel("voltage")
plt.title("ME vs Voltage (pm/V)")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.legend()
plt.show()