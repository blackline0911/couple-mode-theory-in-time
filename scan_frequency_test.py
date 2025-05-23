from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *
import time as timer
from Heater import Heater


# 這個模擬檔是用來與不同比對cmt的架構、參數算出來的Transmission的
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
mode = "scan_frequency"
# mode = "voltage_drive"
wl_in =  1309.0917632973328/1000
Pin = 1 #mW
FSR = 0.002258
radius = 50
energy_RoundTripLoss_data = np.ones(5)*np.exp(-30*2*np.pi*radius*1e-4)
neff_pdk =  np.ones(6)*2.588
mode_area = 0.22*0.5
gamma = np.exp(-30/2*2*np.pi*radius*1e-4)

bit_num = 50
v_bias = -0.5
vpp = 1
Rs = 53.9
Cjs = [23.6e-15, 20e-15]
f_drive=50
level = "NRZ"


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////

a_fit = alpha_fit(RoundTripLoss=energy_RoundTripLoss_data,L = 2*np.pi*radius,input = "energy")
n_fit = neff_fit(neff_data=neff_pdk)
ring_mod = ring(L=2*np.pi*radius, 
            L_active=2*np.pi*radius, 
            alpha=a_fit.alpha_V,
            neff=n_fit.neff_V,
            cross_section=mode_area,
            lambda_incident=wl_in,
            gamma=[gamma],
            FSR = FSR,
            band="O",
            beta_TPA=7.5,
            sigma_FCA=1.45,
            n2 = 4.5e-9,
            tau_eff=10,
            FSR_shift=0,
            FCA_fit_factor=1,
            TPA_fit_factor=1,
            SPM_fit_factor=1,
            )
print("w0 = ",ring_mod.f_res_bar/t0*2*np.pi)
# wl_min =  1.5464
# wl_max =  1.55
# H = Heater(300,2.42,0.5*150/0.6)
H = Heater(300,0,0.5*150/0.6)
wl_min =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 - ring_mod.lambda0/ring_mod.Q/1.5
wl_max =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 + ring_mod.lambda0/ring_mod.Q/2
print("wl_min = ",wl_min ," um")
print("wl_max = ",wl_max ," um")

v = driver(f_drive=f_drive,
           v_bias=v_bias,
           vpp=vpp,
           Rs=Rs,
           raise_cosine=1,
           cj = Cjs,
           PRBS=1,
           level = level)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
os.chdir("./scan_frequency_test/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    vbias = np.array([v_bias+vpp/2,v_bias,v_bias-vpp/2])
    # vbias = np.array([0,-0.5,-1,-1.5,-2])
    vbias = np.arange(0,-0.5,-0.5)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=10000,resolution=1,buffer=100,driver=v)
    print("LineWidth = ",ring_mod.lambda0/ring_mod.Q*1000," nm")
    print(H.P*ring_mod.HE/1e6)
    print("detuning wabvelength = ",(ring_mod.lambda0+H.P*ring_mod.HE/1e6-wl_in)*1000," nm")
    print("dt  = ",t.dt)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    
    T_record = np.zeros( (int(len(t.t_total)-t.buffer*t0/t.dt-1),len(vbias)))
    plt.figure()
    for vb in vbias:
        v.v_bias = vb
        sim.b,s_minus = solving(sim,ring_mod,v,t,H)
        T = Transfer_function(ring_mod,t)
        wl,data = T.mapping(dB(abs(s_minus)**2/sim.Pin))
        wl,data_phase = T.mapping(180/np.pi*np.angle(s_minus))
        plt.plot(wl,data,
                 label="V = "+str(vb))
        T_record[:,int(np.argwhere(vbias==vb))] = data
    plt.grid(color='g',linestyle='--', alpha=0.5)
    plt.xlabel('wavelength(um)')
   
    plt.legend()
    if ring_mod.FCA_fit_factor==1:
        plt.title('Transfer function (with NL absorb) (func fit alpha)')
        plt.savefig("Transmission_vs_voltage (with NL) (func fit alpha)")
    else:
        plt.title('Transfer function (no NL absorb) (func fit alpha)')
        plt.savefig("Transmission_vs_voltage (no NL) (func fit alpha)")
    plt.show()
    ploting(t.t_total,abs(sim.b)**2,x_label="time (ps)",title="Energy in Ring (mJ)")
    ploting(t.t_total,ring_mod.TPA_coeff*abs(sim.b)**2,x_label="time (ps)",title="tpa loss (1/cm)")
    ploting(t.t_total,ring_mod.FCA_coeff*abs(sim.b)**4,x_label="time (ps)",title="fca loss (1/cm)")
    print("ring_mod.TPA_coeff = ",ring_mod.TPA_coeff)
    print("ring_mod.FCA_coeff = ",ring_mod.FCA_coeff)
    sim.save_data(ring_mod,t,v)