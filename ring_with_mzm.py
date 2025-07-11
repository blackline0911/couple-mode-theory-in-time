import numpy as np
from utility import *
import os

def output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3):
    t = (-q3*t1*t2*np.exp(-1j*phi3) + k1*k2*q1*np.exp(-1j*phi1) + q1*q2*q3*(t1**2+k1**2)*(t2**2+k2**2)*np.exp(-1j*(phi1+phi2+phi3))) / \
        (-1 + q1*q2*t1*t2*np.exp(-1j*(phi1+phi2)) - q2*q3*k1*k2*np.exp(-1j*(phi2+phi3)))
    
    return t


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

radius = 7.5
ng = 3.98
lambda0 = 1.55

t1 = 0.94898839
t2 = t1
a = 0.94898839
L1 = 2*np.pi*radius/2
L2 = 2*np.pi*radius/2
L3 = 2*np.pi*radius/2 +L1 + L2
q1 = a**0.5
q2 = a**0.5
q3 = 0.99**0.5
print("t1 at critical couple condition = ",((1+q1*q3)/(1+q3/q1))**0.5)
# t1 = ((1+q1*q3)/(1+q3/q1))**0.5
# t2 = t1
k1 = (1-t1**2)**0.5
k2 = (1-t2**2)**0.5

wl = np.linspace(1.549,1.551,10000)
w = 2*np.pi*c/wl
w0 =  2*np.pi*c/1.55
n0 = 80/7.5*lambda0/(2*np.pi)
n = n0 + (ng-n0)*(w-w0)/w0
print("n0 = ",n0)

phi1 = 2*np.pi/wl*n*L1
phi2 = 2*np.pi/wl*n*L2
phi3 = 2*np.pi/wl*n*L3

t = np.zeros(len(wl)) + 1j*np.zeros(len(wl))
t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3)

print("L3 - L1 = ",L3-L1," um")
print("FSR of MZM = ",1550**2/(ng*1000*abs(L3-L1))," nm")
print("FSR of Ring = ",1550**2/(ng*1000*abs(L2+L1))," nm")
os.chdir("./feedback_ring")

# ploting(wl*1000,dB(abs(t)**2),x_label="wavelength (nm)",title="transmission of feedback ring heater_phase = "+str(heater_phase*180/np.pi),filename="transmission of feedback ring")

# ploting(wl*1000,dB(abs(np.exp(-1j*2*np.pi/wl*n*L1)+np.exp(-1j*2*np.pi/wl*n*L3))**2),x_label="wavelength (nm)",title="transmission of MZM")

# t1 = 0.94898839
# t2 = t1
# k1 = (1-t1**2)**0.5
# k2 = (1-t2**2)**0.5
a = 0.94898839
q1 = a**0.5
q2 = a**0.5
heater_phase = np.linspace(0,2*np.pi,2000)
Deep1 = np.zeros(len(heater_phase))
for i in range(len(heater_phase)):
    t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
    Deep1[i] = np.min(dB(abs(t)**2))

# t1 = 0.9
# t2 = t1
# k1 = (1-t1**2)**0.5
# k2 = (1-t2**2)**0.5
a = 0.94948669
q1 = a**0.5
q2 = a**0.5
heater_phase = np.linspace(0,2*np.pi,2000)
Deep2 = np.zeros(len(heater_phase))
for i in range(len(heater_phase)):
    t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
    Deep2[i] = np.min(dB(abs(t)**2))

# t1 = 0.97
# t2 = t1
# k1 = (1-t1**2)**0.5
# k2 = (1-t2**2)**0.5
a = 0.9498056
q1 = a**0.5
q2 = a**0.5
heater_phase = np.linspace(0,2*np.pi,2000)
Deep3 = np.zeros(len(heater_phase))
for i in range(len(heater_phase)):
    t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
    Deep3[i] = np.min(dB(abs(t)**2))


# ploting(heater_phase*180/np.pi,Deep1,Deep2,Deep3,x_label="heater phase (degree)",title="deep  of feedback ring  vs heater phase (a = 0.94898839)", \
#     filename="deep of transmission of feedback ring with MZM vs heater phase",leg=['t1=0.94898839','t1=0.9','t1=0.97'])
ploting(heater_phase*180/np.pi,Deep1,Deep2,Deep3,x_label="heater phase (degree)",title="deep  of feedback ring  vs heater phase (t1 = 0.94898839)", \
    filename="deep of transmission of feedback ring with MZM vs heater phase",leg=['a=0.94898839','a=0.94948669(1V)','a=0.9498056(2V)'])