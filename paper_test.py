from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *

# refering 
# "Measurement_of_the_Nonlinear_Loss_and_Effective_Free_Carrier_Lifetime_in_Silicon_Microring_Resonators"
wl_in = 1.32105
wl_res = 1.321045
Pin = 10 #mW
kd = sqrt(0.091)
alpha_linear = 0.63
tau_eff=20
radius = 50
# mode = "scan_frequency"
mode = "voltage_drive"

experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)

ring_mod = ring(2*np.pi*radius,  
            np.exp(-alpha_linear/2*2*np.pi*radius*1e-4),
            0,
            0.2,
            wl_in,
            sqrt(1-kd**2),
            sqrt(1-kd**2),
            ng = 3.93,
            # FSR=0.00141,
            lambda0=wl_res,
            tau_eff=tau_eff,
            TPA_fit_factor=1,
            FCA_fit_factor=1)
wl_min =  ring_mod.lambda0 - ring_mod.lambda0/ring_mod.Q*2
wl_max =  ring_mod.lambda0 + ring_mod.lambda0/ring_mod.Q*2 
v = driver(f_drive=100,
           v_bias=0,
           vpp=0,
           R=53.9)
if sim.mode == "scan_frequency":
    t = time(mode = sim.mode)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=5000,resolution=2,buffer=1000)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
else:
    t = time(mode = sim.mode)
    t.main(ring_mod,v,N=10000,resolution=1)

v.create_voltage(time=t)
v.varying_Cj()
sim.save_data(ring_mod,t,v)
b,Q,s_minus,N = solving(sim,ring_mod,v,t)
b0 = np.real(sqrt(t0)*sqrt(sim.Pin))

experiment_condition["Pin"] = 1
sim1 = simulation()
sim1.main(experiment_condition=experiment_condition)
b1,Q1,s1,N1= solving(sim1,ring_mod,v,t)

os.chdir("./paper_test/")
sim.save_data(ring_mod,t,v)


if sim.mode=="voltage_drive":
    ploting(t.t_total,v.v,x_label='time',title='voltage',filename='voltage')
    ploting(t.t_total,N,x_label='time (ps)',title='Free carrier density (cm^-3)',filename='Free_carrier_density')
    ploting(t.t_total,abs(b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')
    
else:
    T = Transfer_function(ring_mod,t)
    wl,data_NL = T.mapping(abs(s_minus)**2/sim.Pin)     
    wl,data = T.mapping(abs(s1)**2/sim1.Pin)     
    ploting(wl*1000,10*np.log10(data_NL),10*np.log10(data),x_label='wavelength (nm)',title='Transfer function',filename='T_v2',leg=['Pin = '+str(sim.Pin),'Pin = '+str(sim1.Pin)])
    ploting(wl*1000,(data_NL),(data),x_label='wavelength (nm)',title='Transfer function',filename='T_v2',leg=['Pin = '+str(sim.Pin),'Pin = '+str(sim1.Pin)])
    ploting(t.t_total,abs(b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')
    ploting(t.t_total,N,x_label='time (ps)',title='Free carrier density (cm^-3)',filename='Free_carrier_density')

    
# ploting(t.t_total,Q,x_label='time (ps)',title='Q',filename='Q_test')
ploting(t.t_total,abs(s_minus)**2,x_label='time (ps)',title=r"$|s_-^2|$",filename='s_minus_power')
# ploting(t.t_total,180/np.pi*np.angle(s_minus),x_label='time (ps)',title='s_minus phase',filename='s_minus phase')

# ploting(t.t_total,(ring_mod.TPA_coeff*t0*sim.Pin*abs(b/b0)**2+ring_mod.FCA_coeff*t0**2*sim.Pin**2*abs(b/b0)**4),x_label='time (ps)',title="NonLinear Loss (1/cm)",filename='NonLinear Loss')
ploting(t.t_total,(ring_mod.TPA_coeff*t0*sim.Pin*abs(b/sim.b0)**2 + N*ring_mod.sigma_FCA*1e-17 ) ,x_label='time (ps)',title="NonLinear Loss (1/cm)",filename='NonLinear Loss')




