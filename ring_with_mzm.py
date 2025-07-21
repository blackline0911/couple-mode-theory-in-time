import numpy as np
from utility import *
import os


os.chdir("./feedback_ring_analysis/")

def output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3):
    t = (-q3*t1*t2*np.exp(-1j*phi3) + k1*k2*q1*np.exp(-1j*phi1) + q1*q2*q3*(t1**2+k1**2)*(t2**2+k2**2)*np.exp(-1j*(phi1+phi2+phi3))) / \
        (-1 + q1*q2*t1*t2*np.exp(-1j*(phi1+phi2)) - q2*q3*k1*k2*np.exp(-1j*(phi2+phi3)))
    
    return t

def map_angle(theta):
    """ Map angle to the range  [0, 360) degrees.

    Args:
        theta (float or np.ndarray): Angle in degree.
    """
    output = np.where(theta < 0,theta + 360, theta)
    return output

def neff_dispersion(wl, ng, lambda0):
    f_bar = 299792458e6/wl/1e12
    f_res_bar = 299792458e6/lambda0/1e12
    return n0 + (ng-n0)*(f_bar/f_res_bar-1)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

radius = 11
ng = 3.98
lambda0 = 1.55

t1 = 0.95
t2 = t1

L1 = 2*np.pi*radius/2
L2 = 2*np.pi*radius/2
junction_radius = 5
La = 2*np.pi*junction_radius
L3 = 2*np.pi*radius/2 +L1 + L2
junction_loss = 0.95223**(junction_radius/5)
bend_loss = 0.999   # 90 degree bend loss
bend_loss_in_dB = dB(bend_loss)  # convert to dB
q1 = bend_loss**2
q2 = junction_loss*dB_inv(bend_loss_in_dB*4*((L2-La)/(2*np.pi*radius)))
print("q2 = ",q2)
q3 = 0.99**0.5
a = q1*q2
# print("t1 at critical couple condition = ",((1+q1*q3)/(1+q3/q1))**0.5)
# t1 = ((1+q1*q3)/(1+q3/q1))**0.5
# t2 = t1
k1 = (1-t1**2)**0.5
k2 = (1-t2**2)**0.5


w0 =  2*np.pi*c/lambda0
n0 = 110/radius*lambda0/(2*np.pi)
# n0 = 2.680503


tau_o = - (L1+L2)*ng/(c*np.log(a))
tau_e1 = - (L1+L2)*ng/(c*np.log(t1))
tau_e2 = - (L1+L2)*ng/(c*np.log(t2))
Q = np.real( ( (2*(1/tau_e1 + 1/tau_e2) + 2/tau_o  )/(w0) )**(-1) )


wl = np.linspace(lambda0-lambda0/Q,lambda0+lambda0/Q/2,10000)
# wl = np.linspace(1.552,1.556,10000)
# wl = np.linspace(1.545,1.555,100000)
w = 2*np.pi*c/wl
n = n0 + (ng-n0)*(w-w0)/w0

phi1 = 2*np.pi/wl*n*L1
phi2 = 2*np.pi/wl*n*L2
phi3 = 2*np.pi/wl*n*L3

t = np.zeros(len(wl)) + 1j*np.zeros(len(wl))
t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3)

print("FSR of MZM = ",1550**2/(ng*1000*abs(L3-L1))," nm")
print("FSR of Ring = ",1550**2/(ng*1000*abs(L2+L1))," nm")

ploting(wl*1000,dB(abs(t)**2),x_label="wavelength (nm)",title="transmission of feedback ring heater_phase = "+str(0*180/np.pi),filename="transmission of feedback ring")

# ploting(wl*1000,dB(abs(np.exp(-1j*2*np.pi/wl*n*L1)+np.exp(-1j*2*np.pi/wl*n*L3))**2),x_label="wavelength (nm)",title="transmission of MZM")
# da = 0.01
# t1 = 0.94898839
# t2 = t1
# k1 = (1-t1**2)**0.5
# k2 = (1-t2**2)**0.5
# a0 = 0.94898839
# q1 = a0**0.5
# q2 = a0**0.5
# heater_phase = np.linspace(0,2*np.pi,2000)
# Deep1 = np.zeros(len(heater_phase))
# for i in range(len(heater_phase)):
#     t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
#     Deep1[i] = np.min(dB(abs(t)**2))

# t1 = 0.9
# t2 = t1
# k1 = (1-t1**2)**0.5
# k2 = (1-t2**2)**0.5
# a1 = 0.94948669+da
# q1 = a1**0.5
# q2 = a1**0.5
# heater_phase = np.linspace(0,2*np.pi,2000)
# Deep2 = np.zeros(len(heater_phase))
# for i in range(len(heater_phase)):
#     t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
#     Deep2[i] = np.min(dB(abs(t)**2))

# t1 = 0.97
# t2 = t1
# k1 = (1-t1**2)**0.5
# k2 = (1-t2**2)**0.5
# a2 = 0.9498056+2*da
# q1 = a2**0.5
# q2 = a2**0.5
# heater_phase = np.linspace(0,2*np.pi,2000)
# Deep3 = np.zeros(len(heater_phase))
# for i in range(len(heater_phase)):
#     t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
#     Deep3[i] = np.min(dB(abs(t)**2))


# ploting(heater_phase*180/np.pi,Deep1,Deep2,Deep3,x_label="heater phase (degree)",title="deep  of feedback ring  vs heater phase (a = 0.94898839)", \
#     filename="deep of transmission of feedback ring with MZM vs heater phase",leg=['t1=0.94898839','t1=0.9','t1=0.97'])
# ploting(heater_phase*180/np.pi,Deep1,Deep2,Deep3,x_label="heater phase (degree)",title="deep  of feedback ring  vs heater phase (t1 = 0.94898839)", \
#     filename="deep of transmission of feedback ring with MZM vs heater phase")


heater_phase = np.array([0,80,110,150,180,210,249,300])*np.pi/180
plt.figure()
for i in range(len(heater_phase)):
    t = output(t1,k1,t2,k2,q1,q2,q3,phi1,phi2,phi3+heater_phase[i])
    plt.plot(wl*1000,dB(abs(t)**2),label=str(int(heater_phase[i]*180/np.pi))+" degree")
plt.xlabel("wavelength (nm)")
plt.ylabel("transmission (dB)")
plt.title("transmission of feedback ring with MZM vs heater phase ")
plt.legend()
plt.grid()
plt.savefig("transmission of feedback ring with MZM vs heater phase .png")
plt.show()

print("round trip amplitude loss = ",q2*q1)
print("Q = ",Q)
print("linewidth = ",(2*np.pi*c*1e3/(w0)**2)*w0/Q," nm")
print("f_opt (optical bandwidth) = ",c/1.55/Q/1e9," GHz")
heater_phase = 110*np.pi/180
dn = 9.670417578102628e-05
V = np.array([0,0.5,1,1.5,2])
phi2 = np.array([n*2*np.pi/wl*La,     (n+dn/2)*2*np.pi/wl*La,\
                          (n+dn)*2*np.pi/wl*La, (n+dn*1.5)*2*np.pi/wl*La,\
                          (n+dn*2)*2*np.pi/wl*La]) + n*2*np.pi/wl*(L2-La)
junction_loss = np.array([0.95223,0.95248,0.95273,0.95292,0.95305])**(junction_radius/5)
q2 = junction_loss*dB_inv(bend_loss_in_dB*4*((L2-La)/(2*np.pi*radius)))
plt.figure()
for i in range(len(V)):
    t = output(t1,k1,t2,k2,q1,q2[i],q3,phi1,phi2[i],phi3+heater_phase)
    plt.plot(wl*1000,dB(abs(t)**2),label=str(-V[i])+" V")
plt.xlabel("wavelength (nm)")
plt.ylabel("transmission (dB)")
plt.title("Transmission of feedback ring with MZM vs voltage (a = "+str(a)+")")
plt.legend(loc=3)
plt.grid()
plt.savefig("Transmission of feedback ring with MZM vs voltage (a = "+str(a)+").png")
plt.show()