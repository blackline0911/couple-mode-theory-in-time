import numpy as np
from utility import *

def output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3):
    t = (-q3*t1*t2*np.exp(-1j*phi3) + k1*k2*q1*np.exp(-1j*phi1) + q1*q2*q3*(t1**2+k1**2)*(t2**2+k2**2)*np.exp(-1j*(phi1+phi2+phi3))) / \
        (-1 + q1*q2*t1*t2*np.exp(-1j*(phi1+phi2)) - q2*q3*k1*k2*np.exp(-1j*(phi2+phi3)))
    
    return t


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

radius = 5
ng = 3.98

t1 = 0.95012
t2 = t1
k1 = (1-t1**2)**0.5
k2 = (1-t2**2)**0.5
L1 = 2*np.pi*radius/2
L2 = 2*np.pi*radius/2
L3 = 2*np.pi*radius/2 + 30
q1 = 0.99**0.5
q2 = 0.99**0.5
q3 = 0.9**0.5
heater_phase = 0

wl = np.linspace(1.5,1.6,10000)
w = 2*np.pi*c/wl
w0 =  2*np.pi*c/1.55
n0 = 2.5
n = n0 + (ng-n0)*(w-w0)/w0

phi1 = 2*np.pi/wl*n*L1
phi2 = 2*np.pi/wl*n*L2
phi3 = 2*np.pi/wl*n*L3 + heater_phase

t = np.zeros(len(wl)) + 1j*np.zeros(len(wl))
t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3)

print("L3 - L1 = ",L3-L1," um")
print("FSR of MZM = ",1550**2/(ng*1000*abs(L3-L1))," nm")
print("FSR of Ring = ",1550**2/(ng*1000*abs(L2+L1))," nm")
ploting(wl*1000,dB(abs(t)**2),x_label="wavelength (nm)",title="transmission of feedback ring")
ploting(wl*1000,dB(abs(np.exp(-1j*2*np.pi/wl*n*L1)+np.exp(-1j*2*np.pi/wl*n*L3))**2),x_label="wavelength (nm)",title="transmission of MZM")

