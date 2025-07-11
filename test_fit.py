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
radius = 5.0 #um
Amp_RoundTripLoss_pdk = [0.953553, 0.953921, 0.954268, 0.954521, 0.954725]
a = -2/(2*np.pi*radius*1e-4)*np.log(Amp_RoundTripLoss_pdk)+4.5
Amp_RoundTripLoss_pdk = np.exp(-2*np.pi*radius*1e-4*a/2)
neff_pdk = [2.680466, 2.680503, 2.680534, 2.680565, 2.680591, 2.680612]
mode_area = 0.22*0.5
gamma = 0.947582
mode = "scan_frequency"


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)
a_fit = alpha_fit(RoundTripLoss=Amp_RoundTripLoss_pdk,L = 2*np.pi*radius,La = 2*np.pi*radius,input2 = "amp")
n_fit = neff_fit(neff_data=neff_pdk)
ring_mod = ring(L=2*np.pi*radius, 
            L_active=[2*np.pi*radius],
            alpha=a_fit.alpha_V,
            neff=n_fit.neff_V,
            cross_section=mode_area,
            lambda_incident=wl_in,
            gamma=[gamma],
            FSR = FSR,
            tau_eff=0.5,
            sigma_FCA=1.04,
            beta_TPA=5,
            n2=4.5e-9,
            FSR_shift=0,)

H = Heater(300,1.3325,0.5*150/0.6,tau_th=17.5)

# wl_min =  1.562
# wl_max =  1.557 
wl_min =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 - ring_mod.lambda0/ring_mod.Q*5
wl_max =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 + ring_mod.lambda0/ring_mod.Q*5

# dummy Cox, Rsi, Cp, level
Cjs = [23.6e-15, 20e-15]
v = driver(f_drive=50,
           v_bias=0,
           vpp=0,
           Rs=[53.9],
           cj=Cjs,
           Cox = [15e-15],
           Rsi = [1000],
           Cp = [6e-15],
           level=2)

os.chdir("./test_fit/")
# exp_data = read_excel.load_excel_data("D:/Homework/Master_degree/ring/CMT/nonlinear/ring spectrum.xlsx", '10dbm')
exp_data = read_excel.load_excel_data("D:/Master_degree/paper/微分方程/ring spectrum.xlsx", '10dbm')
wl = exp_data[:,0]
idx_1 = np.argmin(np.abs(wl-1.56*1e-6))
idx_2 = np.argmin(np.abs(wl-1.57*1e-6))


if sim.mode == "scan_frequency":
    t = time(mode = "scan_frequency")
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=10000,resolution=1,buffer=100,driver=v)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    v.create_voltage(time=t)
    sim.b,s_minus = solving(sim,ring_mod,v,t,H)
    T = Transfer_function(ring_mod,t)
    wl_sim,data = T.mapping(10*np.log10(abs(s_minus)**2/Pin))
sim.save_data(ring_mod,t,v,H)

# ploting(t.t_total,abs(sim.b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')


# """ Estimating Self phase modulation"""
# n2 = 5.6e-9 # um^2/mW
# dn_SPM = n2*abs(b)**2/(ring_mod.round_trip_time*(1e-12)*ring_mod.cross_section)
# # dw = -dn_SPM*2*np.pi*ring_mod.f_res_bar/(ring_mod.neff+dn_SPM)
## df = -dn_SPM*ring_mod.f_res_bar/ring_mod.ng
# ploting(t.t_total, df,x_label='time (ps)',title="w0 changed by Self Phase modulation (THz)",filename='Self_Phase_modulation')
# ploting(t.t_total, df/(ring_mod.f_res_bar),x_label='time (ps)',title="Self Phase modulation (compared to original w0)",filename='Self_Phase_modulation_compared')
# ploting(t.t_total, df/(c/(ring_mod.ng*ring_mod.L)*t0),x_label='time (ps)',title="Self Phase modulation (compared to FSR)",filename='Self_Phase_modulation_compared_FSR')

# Pin = 1
# experiment_condition["Pin"]=Pin
# sim1 = simulation() 
# sim1.main(experiment_condition=experiment_condition)
# b1,s1= solving(sim,ring_mod,v,t,H)
# if sim1.mode =="scan_frequency":
#     T1 = Transfer_function(ring_mod,t)
#     wl,data_Pin_1 = T1.mapping(10*np.log10(abs(s1)**2/sim1.Pin))
#     wl,data_Pin_1_phase = T1.mapping(180/np.pi*np.angle(s1))
# ploting(wl*1000,(data),data_Pin_1,x_label='wavelength (nm)',title='Transfer function',filename='T2',leg=['Pin = '+str(sim.Pin),'Pin = '+str(sim1.Pin)])
# ploting(wl*1000,(data),data_Pin_1,x_label='wavelength (nm)',title='Transfer function',filename='T2',leg=['Pin = '+str(sim.Pin),'Pin = '+str(sim1.Pin)])
# ploting(wl*1000,10**(data/10),10**(data_Pin_1/10),x_label='wavelength (nm)',title='Transfer function',filename='T_dB',leg=['Pin = '+str(sim.Pin),'Pin = '+str(sim1.Pin)])
# ploting(wl[idx_1:idx_2],
#         exp_data[idx_1:idx_2,1]-exp_data[idx_1:idx_2,3],
#         x_label='wl',title='spectrum',filename='exp_data',leg=['0V'])

# plt.figure()
# plt.plot( wl_sim*1000,data,label='simulation')
# # plt.plot(wl[idx_1:idx_2]*1e9,exp_data[idx_1:idx_2,1]-exp_data[idx_1:idx_2,3],label='experiment')
# plt.xlabel('wavelength (nm)')
# plt.ylabel('Transmission (dB)')
# plt.title('Transfer function')
# plt.legend()
# plt.savefig('Transfer_function_comparison.png')
# plt.show()


wl_min =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 + ring_mod.lambda0/ring_mod.Q*5
wl_max =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 - ring_mod.lambda0/ring_mod.Q*5
t = time(mode = "scan_frequency")
ring_mod.scan_frequency(wl_min ,wl_max,t)
t.main(ring_mod,t_max=10000,resolution=1,buffer=100,driver=v)
wl_scan =  c/ring_mod.w_res(t.t_total)*t0
v.create_voltage(time=t)
sim.b,s_minus = solving(sim,ring_mod,v,t,H)
T = Transfer_function(ring_mod,t)
wl_sim_reverse,data_reverse = T.mapping(10*np.log10(abs(s_minus)**2/Pin))

plt.figure()
plt.plot( wl_sim*1000,data,label='simulation')
plt.plot( wl_sim_reverse*1000,data_reverse,label='scan reverse')
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.xlabel('wavelength (nm)')
plt.ylabel('Transmission (dB)')
plt.title('Transfer function')
plt.legend()
plt.savefig('Transfer_function_comparison.png')
plt.show()
