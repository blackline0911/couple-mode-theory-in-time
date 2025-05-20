# refer：Multibistability and self-pulsation in nonlinear high-Q silicon microring resonators considering
from utility import ploting
import numpy as np

c=299792458e6 #cm/s
R = 50  #um
Pin = 1 #mW
n0 = 2.588
lambda0 = 2*np.pi*R*n0/620  #um
w0 = 2*np.pi*c/lambda0
alpha_ring = 0.16  #1/cm
alpha_c = alpha_ring #critical couple
tau_ph = 0.27 #ns (energy life time) c*(alpha_ring+alpha_c)/n0
eta_lin = 0.4
Q = ((1/(tau_ph*1e-9))/w0)**-1
print("Q = ",Q)
print("Linewidth = ",lambda0/Q*1000," nm")
print("lambda0 = ",lambda0)
print("f0 = ",w0/2/np.pi)
wl_in = 
wL = 



# //////////////////////////////////////////////////////
n2 = 4.5e-18
beta2 = 0.75e-11
sigma_FCA = 14.5e-22
Aeff=0.204e-12
Atpa = 0.1289e-12
Afca = 0.116e-12
rho_si = 2.329e6
csi = 0.713
sigma_r1 = 8.8e-22 #cm^3 ? 
sigma_r2 = 8.5e-18 #cm^3 ? 
kappa_theta = 1.86e-4



# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////


alpha = alpha_ring+alpha_c
gamma0 = c*1e-4*alpha/(2*n0)
sigma = sigma_r1
