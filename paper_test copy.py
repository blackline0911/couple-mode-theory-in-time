from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *

# refering:
# Raman gain and nonlinear optical absorption measurements in a low-loss silicon waveguide
wl_in = 1.32105
wl_res = 1.321045
Pin = 2.5 #mW
kd = sqrt(0.091)
alpha_linear = 10**(0.22/10)
tau_eff=14
wl_min =  1.32078
wl_max =  1.32085
# experiment_condition ={"mode":"voltage_drive",
#                         "lambda_incident":wl_in,
#                         "Pin":Pin} 
experiment_condition ={"mode":"scan_frequency",
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)
ring_mod = ring(2*np.pi*50,  
            np.exp(-alpha_linear/2*2*np.pi*50*1e-4),
            0,
            0.2,
            wl_in,
            sqrt(1-kd**2),
            sqrt(1-kd**2),
            ng = 3.93,
            lambda0=wl_res,
            tau_eff=tau_eff,
            sigma_FCA = 1.45,
            )

v = driver(f_drive=100,
           v_bias=0,
           vpp=0,
           R=53.9)
if sim.mode == "scan_frequency":
    t = time(mode = sim.mode)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=3000,resolution=2,buffer=550)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
else:
    t = time(mode = sim.mode)
    t.main(ring_mod,v,N=50,resolution=2)

v.create_voltage(time=t)
v.varying_Cj()
sim.save_data(ring_mod,t,v)

b,Q,s_minus = solving(sim,ring_mod,v,t)
b0 = np.real(sqrt(t0)*sqrt(Pin))

experiment_condition ={"mode":"scan_frequency",
                        "lambda_incident":wl_in,
                        "Pin":1}
sim1 = simulation()
sim1.main(experiment_condition=experiment_condition)
b1,Q1,s1 = solving(sim1,ring_mod,v,t)
os.chdir("./data/")
if sim.mode=="voltage_drive":
    ploting(t.t_total,v.v,x_label='time',title='voltage',filename='voltage')
else:
    T = Transfer_function(ring_mod,t)
    print(sim.Pin)
    wl,data = T.mapping(abs(s_minus)**2/sim.Pin)     
    ploting(wl*1000,data,x_label='wavelength (nm)',title='Transfer function',filename='T_NL')
    # ploting(wl*1000,data_v1,x_label='wavelength (nm)',title='Transfer function',filename='T')
    ploting(t.t_total,abs(b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')
    ploting(t.t_total,wl_scan,x_label='time',title='lambda res (um)',filename='lambda_res')
    
# ploting(t.t_total,Q,x_label='time (ps)',title='Q',filename='Q_test')
# ploting(t.t_total,abs(s_minus)**2,x_label='time (ps)',title=r"$|s_-^2|$",filename='s_minus_power')
# ploting(t.t_total,180/np.pi*np.angle(s_minus),x_label='time (ps)',title='s_minus phase',filename='s_minus phase')
ploting(t.t_total,(ring_mod.TPA_coeff*t0*sim.Pin*abs(b/b0)**2+ring_mod.FCA_coeff*t0**2*sim.Pin**2*abs(b/b0)**4),x_label='time (ps)',title="NonLinear Loss (1/cm)",filename='NonLinear Loss')

np.savetxt('b.txt',b/b0,fmt="%.8f", delimiter="\n")
np.savetxt('T_NL.txt',data,fmt="%.8f", delimiter="\n")
np.savetxt('wl.txt',wl,fmt="%.8f", delimiter="\n")


