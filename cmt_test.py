# refer：Multibistability and self-pulsation in nonlinear high-Q silicon microring resonators considering
from utility import ploting, h
import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt

c=299792458e6 #um/s
R = 50  #um
Pin = 1e-3 #W
n0 = 2.588
# lambda0 = 2*np.pi*R*n0/125  #um
lambda0 = 2*np.pi*R*n0/600  #um
w0 = 2*np.pi*c/lambda0
alpha_ring = 0.16  #1/cm
alpha_c = alpha_ring #critical couple 1/cm
tau_ph = 1/(c*1e-4*(alpha_ring+alpha_c)/n0) # (energy life time) c*(alpha_ring+alpha_c)/n0
eta_lin = 0.4
Q = ((1/(tau_ph))/w0)**-1
print("Q = ",Q)
print("Linewidth = ",lambda0/Q*1000," nm")
print("lambda0 = ",lambda0)
print("f0 = ",w0/2/np.pi)
wl_in = lambda0*1000-lambda0/Q/4*1000 #nm
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

gamma_c = c*1e-4*alpha_c/n0 #energy
L = 2*np.pi*R
alpha = alpha_ring+alpha_c
eta_lin =  0.5 #alpha_abs/alpha_ring

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////
# Normalization

# gamma0 = c*1e-4*alpha/(2*n0)   #1/s amplitude
gamma0 = 1e10
sigma = sigma_r1*wL/(n0*gamma0) #cm^3 * rad
beta = (299792458)**2*beta2/(gamma0*2*h*(wL/2/np.pi)*n0**2*(Afca*L*1e-6)**2) #1/J^2
Kin = (sigma*beta)**0.5*gamma_c/(gamma0**2)
P = Kin*Pin
T_front = wL*kappa_theta/(n0*gamma0) # *dT dummy
tau = gamma0*tau_car
sigma_FCD = wL*sigma_r2/(gamma0*n0*sigma**0.8) # no unit
Gamma_FCA = sigma_FCA*(299792458)/(2*n0*gamma0*sigma*1e-6)
alpha_TPA = beta2*(299792458**2)/(2*gamma0*n0**2*Atpa*L*1e-6* (sigma*1e-6*beta)**0.5)
nkerr = wL*n2*(299792458)/(gamma0*n0**2*Akerr*L*1e-6*(sigma*1e-6*beta)**0.5)
epsilon_T = wL*kappa_theta/(n0*gamma0*rho_si*csi*Aeff*L*1e-6*(sigma*1e-6*beta)**0.5)
eta_c = 2*alpha_ring/alpha

P = Kin*Pin
tau_theta = tau_th*gamma0



# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////

buffer = 5000
wl_start = lambda0*1000 - lambda0/Q*1000*10
wl_end = lambda0*1000 + lambda0/Q*1000*10

f_in_bar = wL/gamma0/2/np.pi
f0_bar = w0/gamma0/2/np.pi
f_start_bar = c*1e3/wl_start/gamma0
f_end_bar = c*1e3/wl_end/gamma0
tmax = 30000
d = (c*1e-4*alpha/2/n0)/gamma0
def cmt(t, eqs):
    if t<(buffer):
        fres = f0_bar+(f_in_bar-f_start_bar)
    else:
        fres = -(f_end_bar-f_start_bar)/tmax*(t-buffer)+ f0_bar+ (f_in_bar-f_start_bar)
    delta = (fres - f_in_bar)
    a, n, T = eqs
    da_dt = ( 1j*delta - 1j*(0*nkerr*abs(a)**2 - 0*(n + sigma_FCD*n**0.8) + 0*T) \
            \
            - ( d + alpha_TPA*abs(a)**2 + 0*Gamma_FCA*n ) )*a \
            \
            +(P)**0.5
    
    dn_dt = abs(a)**4 - n/tau

    dT_dt = epsilon_T*abs(a)**2 *(eta_lin*eta_c + 2*alpha_TPA*abs(a)**2 + 2*Gamma_FCA*n) -T/tau_theta   

    return [da_dt, dn_dt, dT_dt]

a_init = 0+1j*0
n_init = 0
T_init = 0

dt = 1e-12
print("delta t in normalized time = ",dt*gamma0)
t = np.arange(0,int(tmax+buffer),dt*gamma0)
sol = solve_ivp(cmt, [0,int(tmax+buffer)], [a_init, n_init, T_init], method="RK45", t_eval=t,atol=1e-20,rtol=1e-15)

a = sol.y[0]
n = sol.y[1]
T = sol.y[2]
u = a/(sigma*1e-6*beta)**0.25
N = sigma*n
tpa_loss = np.max(beta2*(299792458)/(n0*Atpa*L*1e-6)*abs(u)**2  *1e-2) #1/cm
fca_loss = np.max( sigma_FCA*1e6*N)        #1/cm
print("max N = ",)
print(np.max(abs(a)**2))
print(np.max(abs(u*(sigma*1e-6*beta)**0.25)**2))
spm_shift = np.max( wL*n2*(299792458)/(n0**2*Akerr*L*1e-6) *abs(u)**2 )*(2*np.pi*c*1e3/(wL)**2)

s_minus = (Pin)**0.5 - (gamma_c)**0.5*sol.y[0] /(sigma*beta)**0.25


def w_res(t):
        """
        return the resonant frequency according to specified time t
        """
        if isinstance(t,np.ndarray):
            ans = -(f_end_bar-f_start_bar)/tmax*(t-buffer) + f0_bar+ (f_in_bar-f_start_bar)
            return np.where(t<buffer, f0_bar+(f_in_bar-f_start_bar), ans)
# plt.plot(t,w_res(t))
# plt.show()

with open('file.txt', 'w') as f:
    x = dict(globals())
    for name, val in x.items():
        f.write( name+ " = " + str(val)+"\n" )
def mapping(data):
        """
        Convert the frequency scanning time domain signal to transfer function of ring
        """
        i = int( np.argwhere(t==buffer)[0] )
        f_res_bar = w_res(t)
        L = len(f_res_bar)
        wl = c*1e3/(f0_bar + f_in_bar -f_res_bar)/gamma0
        wl = wl[i:L-1]
        data = data[i:L-1]
        return wl, data
wl, Trans = mapping(10*np.log10(abs( s_minus)**2/Pin))
ploting(wl, Trans, x_label="time (normalized by gamma0)",title="Transmission scanning")
ploting(t, abs(u)**2, x_label="time (normalized by gamma0)",title="Energy in Ring (J)")



