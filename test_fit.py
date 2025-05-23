from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *
import read_excel

wl_in = 1.558
Pin = 10 #mW
FSR = 0.019431
mode = "scan_frequency"
experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)
ring_mod = ring(2*np.pi*5, 
            (0.953553),
            39.3,
            0.5*0.22,
            wl_in,
            (0.947582),
            FSR = FSR,
            neff=2.680503,
            tau_eff=20,
            sigma_FCA=1.04,
            beta_TPA=5,
            FSR_shift=0,)

wl_min =  1.562
wl_max =  1.557 
# wl_min =  ring_mod.lambda0 - ring_mod.lambda0/ring_mod.Q/2
# wl_max =  ring_mod.lambda0 + ring_mod.lambda0/ring_mod.Q/2

v = driver(f_drive=50,
           v_bias=0,
           vpp=0,
           R=53.9)

os.chdir("./test_fit/")
#exp_data = read_excel.load_excel_data("D:/Homework/Master_degree/ring/CMT/nonlinear/ring spectrum.xlsx", '10dbm')
# exp_data = read_excel.load_excel_data("D:/Master_degree/paper/微分方程/ring spectrum.xlsx", '10dbm')
# wl = exp_data[:,0]
# idx_1 = np.argmin(np.abs(wl-1.56*1e-6))
# idx_2 = np.argmin(np.abs(wl-1.57*1e-6))
# ploting(wl[:],
#         exp_data[:,1]-exp_data[:,3],
#         x_label='wl',title='spectrum',filename='exp_data_all')

# ploting(wl[idx_1:idx_2],
#         exp_data[idx_1:idx_2,1]-exp_data[idx_1:idx_2,3],
#         x_label='wl',title='spectrum',filename='exp_data',leg=['0V'])


if sim.mode == "scan_frequency":
    t = time(mode = "scan_frequency")
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=5000,resolution=2,buffer=100,driver=v)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    v.create_voltage(time=t)
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    T = Transfer_function(ring_mod,t)
    wl,data = T.mapping(10*np.log10(abs(s_minus)**2/Pin))
    wl,data_phase = T.mapping(180/np.pi*np.angle(s_minus))
if sim.mode == "voltage_drive":
    t = time(mode = "voltage_drive")
    t.main(ring_mod,N=100,driver=v)
    v.create_voltage(time=t)
    v.varying_Cj()
    b,Q,s_minus = solving(sim,ring_mod,v,t)
    ploting(t.t_total,Q,x_label='time (ps)',title='Q',filename='Q_test')

sim.save_data(ring_mod,t,v)

ploting(t.t_total,abs(b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')
ploting(t.t_total,N,x_label='time (ps)',title='free_carrier_density (1/cm^3)',filename='free_carrier_density')
ploting(t.t_total,(ring_mod.TPA_coeff*abs(b)**2 + N*ring_mod.sigma_FCA*1e-17 ) ,x_label='time (ps)',title="NonLinear Loss (1/cm)",filename='NonLinear Loss')


""" Estimating Self phase modulation"""
n2 = 5.6e-9 # um^2/mW
dn_SPM = n2*abs(b)**2/(ring_mod.round_trip_time*(1e-12)*ring_mod.cross_section)
# dw = -dn_SPM*2*np.pi*ring_mod.f_res_bar/(ring_mod.neff+dn_SPM)
df = -dn_SPM*ring_mod.f_res_bar/ring_mod.ng
ploting(t.t_total, df,x_label='time (ps)',title="w0 changed by Self Phase modulation (THz)",filename='Self_Phase_modulation')
ploting(t.t_total, df/(ring_mod.f_res_bar),x_label='time (ps)',title="Self Phase modulation (compared to original w0)",filename='Self_Phase_modulation_compared')
ploting(t.t_total, df/(c/(ring_mod.ng*ring_mod.L)*t0),x_label='time (ps)',title="Self Phase modulation (compared to FSR)",filename='Self_Phase_modulation_compared_FSR')

Pin = 1
experiment_condition["Pin"]=Pin
sim1 = simulation() 
sim1.main(experiment_condition=experiment_condition)
b1,Q1,s1,N1= solving(sim1,ring_mod,v,t)
if sim1.mode =="scan_frequency":
    T1 = Transfer_function(ring_mod,t)
    wl,data_Pin_1 = T1.mapping(10*np.log10(abs(s1)**2/sim1.Pin))
    wl,data_Pin_1_phase = T1.mapping(180/np.pi*np.angle(s1))
ploting(wl*1000,(data),data_Pin_1,x_label='wavelength (nm)',title='Transfer function',filename='T2',leg=['Pin = '+str(sim.Pin),'Pin = '+str(sim1.Pin)])
# ploting(wl*1000,10**(data/10),10**(data_Pin_1/10),x_label='wavelength (nm)',title='Transfer function',filename='T_dB',leg=['Pin = '+str(sim.Pin),'Pin = '+str(sim1.Pin)])

