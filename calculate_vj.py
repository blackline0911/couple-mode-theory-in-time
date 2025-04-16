import numpy as np
import matplotlib.pyplot as plt

def Cj_V(vol):
        return 3.7675e-14/( (2.5485-vol)**0.5 )

def Q_V(vol):
        return 3.7675e-14/( (2.5485-vol)**0.5 )*vol

def V_Q(Q):
        return -Q**2/(2*3.7675e-14**2) + Q/(2*3.7675e-14**2)*(Q**2 + 4*2.5485*3.7675e-14**2)**0.5
        # return -Q**2/(2*3.7675e-14**2) - Q/(2*3.7675e-14**2)*(Q**2 + 4*2.5485*3.7675e-14**2)**0.5

v = np.linspace(-2,0,1000)
Q = np.linspace(-23e-15,-18e-15,1000)
plt.plot(v,Q_V(v))
plt.plot(v,Cj_V(v))
plt.show()
plt.plot(Q,V_Q(Q))
plt.show()
