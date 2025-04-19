import numpy as np
import matplotlib.pyplot as plt
from driver import driver

v_bias = -3
vpp = 5
v = driver(f_drive=50,
           v_bias=v_bias,
           vpp=vpp,
           R=53.9,
           raise_cosine=1,
           sine_wave=0,
           PRBS=1)
cj_normalizing = v.Cj_V(v_bias)
def Cj_V(vol):
        return 3.7675e-14/( (2.5485-vol)**0.5 )

def Q_V(vol):
        return 3.7675e-14/( (2.5485-vol)**0.5 )*vol


A = 3.7675e-14**2/(cj_normalizing**2)
B = cj_normalizing**2 / (2*3.7675e-14**2)
def V_Q(Q_bar):
        return  (-Q_bar**2*B + \
                B*Q_bar*(Q_bar**2 + 4*2.5485*A)**0.5 )

v = np.linspace(-20,0,10000)
# Q = np.linspace(-100e-15,180e-15,1000)
Q = np.linspace(-200e-15,-1e-15,10000)
plt.plot(v,Cj_V(v),label='Cj')
plt.xlabel('V')
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

plt.plot(v,Q_V(v),label="Q")
plt.xlabel('V')
plt.legend()
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.show()

plt.plot(v,Cj_V(v)/cj_normalizing,label="Large Signal")
plt.plot(v,np.ones((len(v))),label="Small Signal")
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.legend()
plt.title("Comparing to fixed Cj")
plt.xlabel('V')
plt.show()

plt.plot(Q,V_Q(Q/cj_normalizing),label="Large Signal")
plt.plot(Q,(Q/cj_normalizing),label="Small Signal")
plt.grid(color='g',linestyle='--', alpha=0.5)
name = "Vj at Vbias = "+str((v_bias))
plt.title(name)
plt.legend()
plt.xlabel('Q')
plt.show()

dQ = Q[1]-Q[0]
A = V_Q(Q/cj_normalizing)
dA_dQ = np.zeros(len(A))
dA_dQ[0] = (-3*A[0]+4*A[1]-A[2])/2/dQ
dA_dQ[-1] = (3*A[-1]-4*A[-2]+A[-3])/2/dQ
for i in range(1,len(A)-1):
        dA_dQ[i] =  (A[i+1]-A[i-1])/(2*dQ)
plt.plot(Q,dA_dQ)
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.xlabel('Q')
plt.title('dA_dQ')
plt.show()

dQ = Q[1]-Q[0]
A = dA_dQ
d2A_dQ2 = np.zeros(len(A))
d2A_dQ2 [0] = (-3*A[0]+4*A[1]-A[2])/2/dQ
d2A_dQ2 [-1] = (3*A[-1]-4*A[-2]+A[-3])/2/dQ
for i in range(1,len(A)-1):
        d2A_dQ2 [i] =  (A[i+1]-A[i-1])/(2*dQ)
plt.plot(Q,d2A_dQ2 )
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.xlabel('Q')
plt.title('d2A_dQ2 ')
plt.show()
