import numpy as np
from utility import *

def output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3):
    t = (-q3*t1*t2*np.exp(-1j*phi3) + k1*k2*q1*np.exp(-1j*phi1) + q1*q2*q3*(t1**2+k1**2)*(t2**2+k2**2)*np.exp(-1j*(phi1+phi2+phi3))) / \
        (-1 + q1*q2*t1*t2*np.exp(-1j*(phi1+phi2)) - q2*q3*k1*k2*np.exp(-1j*(phi2+phi3)))
    
    return t


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

radius = 7.5
ng = 3.98

t1 = 0.94898839
t2 = t1

L1 = 2*np.pi*radius/2
L2 = 2*np.pi*radius/2
L3 = 2*np.pi*radius/2 + 60
q1 = 0.94898839**0.5
q2 = 0.94898839**0.5
q3 = 0.9**0.5
print("t1 at critical couple condition = ",((1+q1*q3)/(1+q3/q1))**0.5)
# t1 = ((1+q1*q3)/(1+q3/q1))**0.5
# t2 = t1
k1 = (1-t1**2)**0.5
k2 = (1-t2**2)**0.5
heater_phase = 0

wl = np.linspace(1.5,1.6,100000)
w = 2*np.pi*c/wl
w0 =  2*np.pi*c/1.55
n0 = 2.7
n = n0 + (ng-n0)*(w-w0)/w0
# n = n0

phi1 = 2*np.pi/wl*n*L1
phi2 = 2*np.pi/wl*n*L2
phi3 = 2*np.pi/wl*n*L3 + heater_phase
# phi3 = np.zeros(len(wl))

t = np.zeros(len(wl)) + 1j*np.zeros(len(wl))
t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3)

print("L3 - L1 = ",L3-L1," um")
plt.scatter(wl,phi3,label="phi3")
plt.show()
plt.scatter(wl,phi3%(2*np.pi)*180/np.pi,label="phi3 (in 2pi)")
plt.show()
print("FSR of MZM = ",1550**2/(ng*1000*abs(L3-L1))," nm")
print("FSR of Ring = ",1550**2/(ng*1000*abs(L2+L1))," nm")
ploting(wl*1000,dB(abs(t)**2),x_label="wavelength (nm)",title="transmission of feedback ring")
ploting(wl*1000,dB(abs(np.exp(-1j*2*np.pi/wl*n*L1)+np.exp(-1j*2*np.pi/wl*n*L3))**2),x_label="wavelength (nm)",title="transmission of MZM")

heater_phase = np.linspace(0,2*np.pi,15)
for i in range(len(heater_phase)):
    t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
    plt.plot(wl*1000,dB(abs(t)**2),label=f"heater phase = {(heater_phase[i])*180/np.pi:.1f} deg")
plt.xlabel("wavelength (nm)")
plt.ylabel("transmission (dB)")
plt.title("Transmission of feedback ring with MZM")
plt.legend()
plt.grid()
plt.show()
heater_phase = np.linspace(0,2*np.pi,500)
y = output(t1,k1,t2,k2,q1,q2,q3,-np.pi,-np.pi,heater_phase)
ploting(heater_phase*180/np.pi,dB(abs(y)**2),x_label="heater_phase (degree)",title="transmission of feedback ring with MZM at -pi/2 phase")
# x = (-q3*t1*t2 + q1*q2*q3)*np.exp(-1j*phi3)-k1*k2*q1
# plt.plot()