from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *

wl_in = 1.5468
Pin = 1 #mW
FSR = 0.0195
# mode = "scan_frequency"
mode = "voltage_drive"
bit_num = 20
v_bias = -1
vpp = 1
experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)
ring_mod = ring(2*np.pi*5, 
            (0.95124),
            39.5,
            0.5*0.22,
            wl_in,
            (0.95105),
            FSR = FSR,
            neff=2.5111,
            FSR_shift=0,)

wl_min =  ring_mod.lambda0 - ring_mod.lambda0/ring_mod.Q/2
wl_max =  ring_mod.lambda0 + ring_mod.lambda0/ring_mod.Q/2 

v = driver(f_drive=50,
           v_bias=v_bias,
           vpp=vpp,
           R=53.9,
           raise_cosine=1,
           sine_wave=0,
           PRBS=1)

os.chdir("./eye_diagram_test/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=5000,resolution=1,buffer=100)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    v.create_voltage(time=t)
    sim.save_data(ring_mod,t,v)
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    T = Transfer_function(ring_mod,t)
    wl,data = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    wl,data_phase = T.mapping(180/np.pi*np.angle(s_minus))
    v.v_bias=-1.5
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    T = Transfer_function(ring_mod,t)
    wl,data_v1 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    wl,data_phase_v1 = T.mapping(180/np.pi*np.angle(s_minus))
    v.v_bias=0
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    T = Transfer_function(ring_mod,t)
    wl,data_v2 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    wl,data_phase_v2 = T.mapping(180/np.pi*np.angle(s_minus))
    v.v_bias=-0.5
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    T = Transfer_function(ring_mod,t)
    wl,data_v3 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    wl,data_phase_v3 = T.mapping(180/np.pi*np.angle(s_minus))
    ploting(wl,data_v2,data_v3,data,data_v1 ,x_label='time (ps)',title='Transfer function',filename='Transfer function_vs_voltage',leg=['V=0','V=-0.5','V=-1','V=-1.5'])
    
if sim.mode == "voltage_drive":
    sim.save_data(ring_mod,t,v)
    t.main(ring_mod,N=bit_num,driver=v)
    v.create_voltage(time=t)
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    ploting(t.t_total,v.v,x_label='time (ps)',title='voltage (V)',filename='voltage')
    ploting(t.t_total,v.Cj,x_label='time (ps)',title='capacitance (C)',filename='capacitance')
    v.varying_Cj()
    b1,Q1,s_minus1,N1 = solving(sim,ring_mod,v,t)
    ploting(t.t_total,Q,Q1,x_label='time (ps)',title='Q',filename='Q',leg=['fixed Cj','varying Cj'])
    ploting(t.t_total,v.Cj,x_label='time (ps)',title='capacitance (C)',filename='capacitance_varying')
    voltage = np.linspace(-2,0,1000)
    ploting(voltage,v.Cj_V(voltage),x_label='time (ps)',title='capacitance vs voltage',filename='capacitance _vs_voltage')
    ploting(t.t_total,abs(s_minus)**2,abs(s_minus1)**2,x_label='time (ps)',title='output Power',filename='output Power')
    ploting(t.t_total,10*np.log10(abs(s_minus)**2),10*np.log10(abs(s_minus1)**2),x_label='time (ps)',title='output Power (dB)',filename='output Power dB')
    ploting(t.t_total,180/np.pi*np.angle(s_minus),180/np.pi*np.angle(s_minus1),x_label='time (ps)',title='s_minus phase',filename='s_minus phase',leg=['varying Cj','Cj=20fF'])
    sim.eye_diagram(t,v,abs(s_minus)**2,filename='eye')
    sim.eye_diagram(t,v,abs(s_minus1)**2,filename='eye2')

ploting(t.t_total,(ring_mod.TPA_coeff*t0*sim.Pin*abs(b/sim.b0)**2 + N*ring_mod.sigma_FCA*1e-17 ) ,x_label='time (ps)',title="NonLinear Loss (1/cm)",filename='NonLinear Loss')
