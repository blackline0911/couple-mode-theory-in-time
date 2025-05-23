# refer：Multibistability and self-pulsation in nonlinear high-Q silicon microring resonators considering
from utility import ploting, h
import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt

c=299792458e6 #um/s
R = 50  #um
Pin = 1e-3 #W
n0 = 2.588
ng = n0
# lambda0 = 2*np.pi*R*n0/125  #um
lambda0 = 2*np.pi*R*n0/621  #um
w0 = 2*np.pi*c/lambda0
alpha_ring = 30  #1/cm
alpha_c = alpha_ring #critical couple 1/cm
tau_ph = 1/(c*1e-4*(alpha_ring+alpha_c)/ng) # (energy life time) c*(alpha_ring+alpha_c)/n0
eta_lin = 0.4
Q = ((1/(tau_ph))/w0)**-1
print("Q = ",Q)
print("Linewidth = ",lambda0/Q*1000," nm")
print("lambda0 = ",lambda0)
print("f0 = ",w0/2/np.pi)
wl_in = 1309.0917632973328 # nm
wL = 2*np.pi*c*1e3/wl_in



# //////////////////////////////////////////////////////
n2 = 4.5e-18            #m^2/W
beta2 = 0.75e-11        #m/W
sigma_FCA = 14.5e-22    #m^2
Aeff=0.204e-12          #m^2
Atpa = 0.1289e-12       #m^2
Afca = 0.116e-12        #m^2
Akerr = Atpa         
rho_si = 2.329e6        #g/m^3
csi = 0.713             #J/(g*K)
sigma_r1 = 8.8e-22      #cm^3
sigma_r2 = 8.5e-18      #cm^3
kappa_theta = 1.86e-4   #1/K
tau_car = 10e-9         #s
tau_th = 100e-9         #s

gamma_c = c*1e-4*alpha_c/ng #energy
L = 2*np.pi*R
alpha = alpha_ring+alpha_c
eta_lin =  0.5 #alpha_abs/alpha_ring

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////
# Normalization

# gamma0 = c*1e-4*alpha/(2*n0)   #1/s amplitude
gamma0 = 1e12
sigma = sigma_r1*wL/(ng*gamma0) #cm^3 * rad
beta = (299792458)**2*beta2/(gamma0*2*h*(wL/2/np.pi)*ng**2*(Afca*L*1e-6)**2) #1/J^2/m^3
Kin = (sigma*1e-6*beta)**0.5*gamma_c/(gamma0**2)
P = Kin*Pin
T_front = wL*kappa_theta/(ng*gamma0) # *dT dummy
tau = gamma0*tau_car
sigma_FCD = wL*sigma_r2/(gamma0*ng*(sigma)**0.8) # no unit
Gamma_FCA = sigma_FCA*(299792458)/(2*ng*gamma0*sigma*1e-6)
alpha_TPA = beta2*(299792458**2)/(2*gamma0*ng**2*Atpa*L*1e-6* (sigma*1e-6*beta)**0.5)
nkerr = wL*n2*(299792458)/(gamma0*ng**2*Akerr*L*1e-6*(sigma*1e-6*beta)**0.5)
epsilon_T = wL*kappa_theta/(ng*gamma0*rho_si*csi*Aeff*L*1e-6*(sigma*1e-6*beta)**0.5)
eta_c = 2*alpha_ring/alpha

P = Kin*Pin
tau_theta = tau_th*gamma0



# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////

buffer = 100
wl_start = lambda0*1000 - lambda0/Q*1000/1.5
wl_end = lambda0*1000 + lambda0/Q*1000/2

f_in_bar = wL/gamma0/2/np.pi
f0_bar = w0/gamma0/2/np.pi
f_start_bar = c*1e3/wl_start/gamma0
f_end_bar = c*1e3/wl_end/gamma0
tmax = 10000
d = (c*1e-4*alpha/2/ng)/gamma0

print(1/(c*1e-4*alpha/2/ng))
print("wl_in = ",wl_in," um")
print("ng = ",ng)
def cmt(t, eqs):
    if t<(buffer):
        fres = f0_bar+(f_in_bar-f_start_bar)
    else:
        fres = -(f_end_bar-f_start_bar)/tmax*(t-buffer)+ f0_bar+ (f_in_bar-f_start_bar)
    delta = (fres - f_in_bar)
    # a, n, T = eqs
    a = eqs
    n = tau*abs(a)**4
    T=0
    da_dt = ( 1j*delta - 1j*(nkerr*abs(a)**2 - 0*(n + sigma_FCD*n**0.8) + 0*T) \
            \
            - ( d + 0*alpha_TPA*abs(a)**2 +0* Gamma_FCA*n ) )*a \
            \
            +(P)**0.5
    
    # dn_dt = abs(a)**4 - n/tau

    # dT_dt = epsilon_T*abs(a)**2 *(eta_lin*eta_c + 2*alpha_TPA*abs(a)**2 + 2*Gamma_FCA*n) -T/tau_theta   

    # return [da_dt, dn_dt, dT_dt]
    return [da_dt]

a_init = 0+1j*0
n_init = 0
T_init = 0

dt = 1e-12
print("delta t in normalized time = ",dt*gamma0)
t = np.arange(0,int(tmax+buffer),dt*gamma0)
sol = solve_ivp(cmt, [0,int(tmax+buffer)], [a_init, n_init, T_init], method="RK45", t_eval=t,atol=1e-20,rtol=1e-15)

a = sol.y[0]
# n = sol.y[1]
# T = sol.y[2]
u = a/(sigma*1e-6*beta)**0.25
# N = n/sigma         #1/cm
# delta_T = n0*gamma0/(wL*kappa_theta)*T
# tpa_loss = np.max(beta2*(299792458)/(ng*Atpa*L*1e-6)*abs(u)**2  *1e-2) #1/cm
# fca_loss = np.max( sigma_FCA*1e4*N)       #1/cm
# print("linear loss = ", d*gamma0*ng/(c*1e-4)," 1/cm")
# print(beta2*(299792458)/(ng*Atpa*L*1e-6))
# print("max N = ",np.max(N)," 1/cm^3")
# print(np.max(abs(u)**2))
# print(np.max(abs(u*(sigma*1e-6*beta)**0.25)**2))
# spm_shift = np.max( wL*n2*(299792458)/(ng**2*Akerr*L*1e-6) *abs(u)**2 )*(2*np.pi*c*1e3/(wL)**2)

s_minus = (Pin)**0.5 - (gamma_c)**0.5*u

print(abs(s_minus)**2)
def f_res(t):
        """
        return the resonant frequency according to specified time t
        """
        if isinstance(t,np.ndarray):
            ans = -(f_end_bar-f_start_bar)/tmax*(t-buffer) + f0_bar+ (f_in_bar-f_start_bar)
            return np.where(t<=buffer, f0_bar+(f_in_bar-f_start_bar), ans)
# plt.plot(t,f_res(t)/1e12 + (n + sigma_FCD*n**0.8)/2/np.pi*gamma0/1e12)
# plt.show()

with open('file.txt', 'w') as f:
    x = dict(globals())
    for name, val in x.items():
        f.write( name+ " = " + str(val)+"\n" )
def mapping(data):
        """
        Convert the frequency scanning time domain signal to transfer function of ring
        """
        i = int( np.argwhere(abs(t-buffer)<dt*gamma0/10)[0] )
        f_res_bar = f_res(t)
        L = len(f_res_bar)
        wl = c*1e3/(f0_bar + f_in_bar -f_res_bar)/gamma0
        wl = wl[i:L-1]
        data = data[i:L-1]
        return wl, data
import os
os.chdir("./cmt_test/")
wl, Trans = mapping(10*np.log10(abs( s_minus)**2/Pin))
# ploting(t/gamma0,10*np.log10(abs( s_minus)**2/Pin) , x_label="time (normalized by gamma0)",title="Transmission scanning",filename="Transmission scanning")
ploting(wl,Trans , x_label="wavelength (nm)",title="Transmission scanning")
ploting(t/gamma0, abs(u)**2, x_label="time (normalized by gamma0)",title="Energy in Ring (J)",filename="Energy in Ring")
# ploting(t/gamma0, delta_T, x_label="time (normalized by gamma0)",title="T (Kelvin)",filename="T (Kelvin)")
# ploting(t/gamma0, N, x_label="time (normalized by gamma0)",title="Free Carrier conctration (1/cm^3)",filename="Free Carrier conctration)")
ploting(t/gamma0,beta2*(299792458)/(ng*Atpa*L*1e-6)*abs(u)**2*1e-2,x_label="t",title="tpa loss",filename="tpa loss")
ploting(t/gamma0,(sigma_FCA)*((299792458)**2*beta2*tau_car/(ng**2*2*h*(wL/2/np.pi)*(Afca*L*1e-6)**2))*abs(u)**4*1e-2,\
        x_label="t",title="FCA loss (1/cm)",filename="FCA loss")
print("TPA coeff = ",beta2*(299792458)/(ng*Atpa*L*1e-6)*1e-2," 1/cm/J")
print("FCA coeff = ",(sigma_FCA)*((299792458)**2*beta2*tau_car/(ng**2*2*h*(wL/2/np.pi)*(Afca*L*1e-6)**2))*1e-2, " 1/cm/J^2")
# ploting(t/gamma0,alpha_TPA*abs(a)**2 , Gamma_FCA*n,x_label="t",title="Normalized Nonlinear loss ",leg=["TPA","FCA"],filename="Normalized Nonlinear loss ")
