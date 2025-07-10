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

L1 = 2*np.pi*radius/2
L2 = 2*np.pi*radius/2
L3 = 2*np.pi*radius/2 +L2+L1
q1 = 0.9**0.5
q2 = 0.9**0.5
q3 = 0.9**0.5
print("t1 at critical couple condition = ",((1+q1*q3)/(1+q3/q1))**0.5)
# t1 = ((1+q1*q3)/(1+q3/q1))**0.5
# t2 = t1
k1 = (1-t1**2)**0.5
k2 = (1-t2**2)**0.5
heater_phase = np.pi/180*330

wl = np.linspace(1.549,1.551,10000)
w = 2*np.pi*c/wl
w0 =  2*np.pi*c/1.55
n0 = 80/7.5*lambda0/(2*np.pi)
n = n0 + (ng-n0)*(w-w0)/w0
print("n0 = ",n0)

phi1 = 2*np.pi/wl*n*L1
phi2 = 2*np.pi/wl*n*L2
phi3 = 2*np.pi/wl*n*L3 + heater_phase
# phi3 = np.zeros(len(wl))

t = np.zeros(len(wl)) + 1j*np.zeros(len(wl))
t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3)

print("L3 - L1 = ",L3-L1," um")
# plt.scatter(wl,phi3,label="phi3")
# plt.show()
# plt.scatter(wl,phi3%(2*np.pi)*180/np.pi,label="phi3 (in 2pi)")
# plt.show()
print("FSR of MZM = ",1550**2/(ng*1000*abs(L3-L1))," nm")
print("FSR of Ring = ",1550**2/(ng*1000*abs(L2+L1))," nm")
os.chdir("./feedback_ring")
# phi3 = np.zeros(len(wl))
# phi3 = np.zeros(len(wl))
ploting(wl*1000,dB(abs(t)**2),x_label="wavelength (nm)",title="transmission of feedback ring heater_phase = "+str(heater_phase*180/np.pi),filename="transmission of feedback ring")

# ploting(wl*1000,dB(abs(np.exp(-1j*2*np.pi/wl*n*L1)+np.exp(-1j*2*np.pi/wl*n*L3))**2),x_label="wavelength (nm)",title="transmission of MZM")

# heater_phase = np.linspace(0,2*np.pi,15)
# for i in range(len(heater_phase)):
#     t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
#     plt.plot(wl*1000,dB(abs(t)**2),label=f"heater phase = {(heater_phase[i])*180/np.pi:.1f} deg")
# plt.xlabel("wavelength (nm)")
# plt.ylabel("transmission (dB)")
# plt.title("Transmission of feedback ring with MZM")
# plt.legend()
# plt.grid()
# plt.show()
# heater_phase = np.linspace(0,2*np.pi,500)
# y = output(t1,k1,t2,k2,q1,q2,q3,-np.pi,-np.pi,heater_phase)
# ploting(heater_phase*180/np.pi,dB(abs(y)**2),x_label="heater_phase (degree)",title="transmission of feedback ring with MZM at pi phase")
# x = (-q3*t1*t2 + q1*q2*q3)*np.exp(-1j*phi3)-k1*k2*q1
# plt.plot()

x = np.array([0,45,70,90,100,110,120,130,140,160,180,200,220,230,240,250,260,270,290,310,330,350])
y = np.array([-4.5,-5.5,-8,-12.5,-18,-60,-17,-10,-7,-2,0,-2,-6.5,-10,-17,-60,-17.5,-12.5,-8,-6,-5,-4.5])
x1 = np.array([0,45,70,90,100,110,120,130,140,150,170,180,200,220,240,260,280,300,320,340])
y1 = np.array([0,-1,-2,-2.5,-3,-5,-6.5,-10,-20,-15,-1,0,-5.5,-20,-6,-3,-2.5,-2,-1.5,-1])
x2 = np.array([0,45,70,80,90,100,120,140,160,180,200,220,240,260,280,290,300,310,320,330,340,350])
y2 = np.array([-9.5,-12.5,-20,-60,-20,-15,-7.5,-3,-1,0,-1,-3,-7.5,-15,-60,-20,-17.5,-13,-12,-10.5,-10,-9.5])
plt.plot(x,y,label="t1=round-trip loss")
plt.plot(x1,y1,label="t1<round-trip loss")
plt.plot(x2,y2,label="t1>round-trip loss")
plt.grid()
plt.legend()
plt.xlabel("heater phase (degree)")
plt.ylabel("deep (dB)")
plt.title("Transmission of feedback ring with MZM vs heater phase")
plt.savefig("Transmission of feedback ring with MZM vs heater phase")
plt.show()