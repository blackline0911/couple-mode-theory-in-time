import numpy as np
import matplotlib.pyplot as plt

cj = [23.6e-15, 20e-15, 18.5e-15]
a = [0.95124, 0.95248, 0.95273, 0.95292, 0.95305]
L = 2*np.pi*5
LA = L
alpha_pdk = -1/(2*L*1e-4)*np.log(a)
delta_alpha_pdk = alpha_pdk - alpha_pdk[0]
# print(alpha_pdk)

def Cd_V(V):
    return 3.7675e-14/( (2.5485-V)**0.5 )

# delta_alpha = A * Cd(V)
# A = alpha_pdk[0]/
print( delta_alpha_pdk[1]/Cd_V(-0.5)/(-0.5))
print( delta_alpha_pdk[2]/Cd_V(-1)/(-1))
print( delta_alpha_pdk[3]/Cd_V(-1.5)/(-1.5))
print( delta_alpha_pdk[4]/Cd_V(-2)/(-2))



V = np.array([0,-0.5,-1,-1.5,-2])
plt.plot(V,delta_alpha_pdk[1]/Cd_V(-0.5)/(-0.5)* Cd_V(V)*V,label='fit by -0.5V')
plt.plot(V,delta_alpha_pdk[2]/Cd_V(-1)/(-1)* Cd_V(V)*V,label='fit by -1V')
plt.plot(V,delta_alpha_pdk[3]/Cd_V(-1.5)/(-1.5)* Cd_V(V)*V,label='fit by -1.5V')
plt.plot(V,delta_alpha_pdk[4]/Cd_V(-2)/(-2)* Cd_V(V)*V,label='fit by -2V')
plt.scatter(V,alpha_pdk-alpha_pdk[0],c='r',marker='o',label='pdk data')
plt.xlabel('Voltage')
plt.ylabel('Delta alpha')
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)

plt.show()
# def alpha_V(V):
