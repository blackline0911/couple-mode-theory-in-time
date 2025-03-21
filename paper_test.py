from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *

wl_in = 1.32105
# wl_in = 1.5488
Pin = 10 #mW
# experiment_condition ={"mode":"voltage_drive",
#                         "lambda_incident":wl_in,
#                         "Pin":Pin} 
experiment_condition ={"mode":"scan_frequency",
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)
wl_min = 1.32095
wl_max = 1.3211
# wl_min = 1.55
# wl_max = 1.548
ring_mod = ring(L=2*np.pi*50, 
            ng=3.93, 
            gamma = sqrt(1-0.091), 
            # alpha = np.exp(-0.63*2*np.pi*50*1e-4),
            alpha=sqrt(1-0.091),
            me=0,
            cross_section=0.2,
            lambda_incident=wl_in,
            lambda0=1.321045
            )

v = driver(f_drive=100,
           v_bias=0,
           vpp=0,
           R=53.9)
if sim.mode == "scan_frequency":
    t = time(mode = sim.mode)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=2000,resolution=2,buffer=550)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
else:
    t = time(mode = sim.mode)
    t.main(ring_mod,v,N=50,resolution=2)

v.create_voltage(time=t)
v.varying_Cj()

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
    wl,data_NL = T.mapping(abs(s_minus)**2/sim.Pin)     
    wl,data = T.mapping(abs(s1)**2/sim1.Pin)     
    ploting(wl*1000,10*np.log10(data_NL),10*np.log10(data),x_label='wavelength (nm)',title='Transfer function',filename='T')
    # ploting(wl*1000,data_v1,x_label='wavelength (nm)',title='Transfer function',filename='T')
    ploting(t.t_total,abs(b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')
    ploting(t.t_total,wl_scan,x_label='time',title='lambda res (um)',filename='lambda_res')
    
# ploting(t.t_total,Q,x_label='time (ps)',title='Q',filename='Q_test')
# ploting(t.t_total,abs(s_minus)**2,x_label='time (ps)',title=r"$|s_-^2|$",filename='s_minus_power')
# ploting(t.t_total,180/np.pi*np.angle(s_minus),x_label='time (ps)',title='s_minus phase',filename='s_minus phase')
ploting(t.t_total,(ring_mod.TPA_coeff*t0*sim.Pin*abs(b/b0)**2+ring_mod.FCA_coeff*t0**2*sim.Pin**2*abs(b/b0)**4),x_label='time (ps)',title="NonLinear Loss (1/cm)",filename='NonLinear Loss')

np.savetxt('b.txt',b/b0,fmt="%.8f", delimiter="\n")
np.savetxt('T.txt',data_NL,fmt="%.8f", delimiter="\n")
np.savetxt('wl.txt',wl,fmt="%.8f", delimiter="\n")


sim.save_data(ring_mod,t,v)