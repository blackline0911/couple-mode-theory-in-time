import numpy as np

w = 2*np.pi*100e9
C = 1.6e-14
R = 53.9

Vj = 1/(1j*w*R*C+1)
print(180/np.pi*np.angle(Vj))
print(abs(Vj))

import matplotlib.pyplot as plt
t = np.linspace(0,200,100000)
plt.figure()
plt.plot(t,1*np.exp(1j*w/1e12*t))
plt.plot(t,Vj*1*np.exp(1j*w/1e12*t))
plt.show()
