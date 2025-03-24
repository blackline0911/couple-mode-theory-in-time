from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *

wl_in = 1.5488
Pin = 1 #mW
experiment_condition ={"mode":"scan_frequency",
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
# experiment_condition ={"mode":"voltage_drive",
#                         "lambda_incident":wl_in,
#                         "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)
wl_start = 1.55
wl_end = 1.548
ring_mod = ring(L=2*np.pi*5, 
            ng=4.3, 
            gamma = 0.95012, 
            alpha = 0.95246,
            me=37.2,
            cross_section=0.5*0.22,
            lambda_incident=wl_in,
            neff=2.51464,
            )

v = driver(f_drive=100,
           v_bias=0,
           vpp=0,
           R=53.9)
if sim.mode == "scan_frequency":
    t = time(mode = "scan_frequency")
    ring_mod.scan_frequency(wl_start ,wl_end,t)
    t.main(ring_mod,t_max=300,resolution=2)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    v.create_voltage(time=t)
if sim.mode == "voltage_drive":
    t = time(mode = "voltage_drive")
    t.main(ring_mod,N=10,driver=v)
    v.create_voltage(time=t)
    v.varying_Cj()


b,Q,s_minus = solving(sim,ring_mod,v,t)
T = Transfer_function(ring_mod,t)
wl,data = T.mapping(abs(s_minus)**2/Pin)

os.chdir("./data/")
# ploting(t.t_total,ring_mod.w_res(t.t_total),x_label='time',title='f res',filename='f_res')
ploting(t.t_total,abs(b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')
# ploting(t.t_total,wl_scan,x_label='time',title='lambda res (um)',filename='lambda_res')
# ploting(t.t_total,Q,x_label='time (ps)',title='Q',filename='Q_test')

ploting(t.t_total,abs(s_minus)**2,x_label='time (ps)',title=r"$|s_-^2|$ (mW)",filename='s_minus_power')
ploting(t.t_total,abs(s_minus)**2,x_label='time (ps)',title=r"$|s_-^2|$ (mW)",filename='s_minus_power_2')
# ploting(t.t_total,180/np.pi*np.angle(s_minus),x_label='time (ps)',title='s_minus phase',filename='s_minus phase')
# np.savetxt('b.txt',b/sim.b0,fmt="%.8f", delimiter="\n")
ploting(t.t_total,(ring_mod.TPA_coeff*t0*sim.Pin*abs(b/sim.b0)**2+ring_mod.FCA_coeff*t0**2*sim.Pin**2*abs(b/sim.b0)**4),x_label='time (ps)',title="NonLinear Loss (1/cm)",filename='NonLinear Loss')

Pin = 1
experiment_condition ={"mode":"scan_frequency",
                        "lambda_incident":wl_in,
                        "Pin":Pin}
sim1 = simulation() 
sim1.main(experiment_condition=experiment_condition)
b1,Q1,s1 = solving(sim1,ring_mod,v,t)
T1 = Transfer_function(ring_mod,t)
wl,data_Pin_1 = T1.mapping(abs(s1)**2/sim1.Pin)
ploting(wl*1000,(data),(data_Pin_1),x_label='wavelength (nm)',title='Transfer function',filename='T')
# np.savetxt('T.txt',data_v1,fmt="%.8f", delimiter="\n")
# np.savetxt('T_NL_0.txt',data,fmt="%.8f", delimiter="\n")
sim.save_data(ring_mod,t,v)