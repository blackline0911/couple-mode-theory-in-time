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

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# mode = "scan_frequency"
mode = "voltage_drive"
wl_in = 1.54685
Pin = 1 #mW
FSR = 0.0195
radius = 5
Amp_RoundTripLoss_pdk = [0.95124, 0.95248, 0.95273, 0.95292, 0.95305]
neff_pdk = [2.51105, 2.5111, 2.51113, 2.51116, 2.51118, 2.5112]
mode_area = 0.22*0.5

bit_num = 5
v_bias = -1
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

a_fit = alpha_fit(RoundTripLoss=Amp_RoundTripLoss_pdk,L = 2*np.pi*radius)
n_fit = neff_fit(neff_data=neff_pdk)
ring_mod = ring(L=2*np.pi*radius, 
            alpha=a_fit.alpha_V,
            neff=n_fit.neff_V,
            cross_section=mode_area,
            lambda_incident=wl_in,
            gamma=[0.95105],
            FSR = FSR,
            FCA_fit_factor=1
            )

wl_min =  1.5464
wl_max =  1.55
# wl_min =  ring_mod.lambda0 - ring_mod.lambda0/ring_mod.Q/2
# wl_max =  ring_mod.lambda0 + ring_mod.lambda0/ring_mod.Q/2

v = driver(f_drive=f_drive,
           v_bias=v_bias,
           vpp=vpp,
           Rs=Rs,
           raise_cosine=1,
           cj = Cjs,
           PRBS=1,
           level = level)
H = Heater(300,0,0.5*150/0.6)
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
os.chdir("./eye_diagram_test/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    vbias = np.array([v_bias+vpp/2,v_bias,v_bias-vpp/2])
    vbias = np.array([0,-0.5,-1,-1.5,-2])
    # vbias = np.arange(0,-0.5,-0.5)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=10000,resolution=1,buffer=100,driver=v)
    print("dt = ",t.dt)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    sim.save_data(ring_mod,t,v)
    T_record = np.zeros( (int(len(t.t_total)-t.buffer*t0/t.dt-1),len(vbias)))
    plt.figure()
    for vb in vbias:
        v.v_bias = vb
        b,s_minus = solving(sim,ring_mod,v,t,H)
        T = Transfer_function(ring_mod,t)
        wl,data = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
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
    ploting(t.t_total,abs(b)**2,x_label="time (ps)",title="Energy in Ring (mJ)")

    
    highlevel_arg = int(np.argwhere(vbias==v_bias+vpp/2))
    lowlevel_arg = int(np.argwhere(vbias==v_bias-vpp/2))
    Extinction_ratio = ER(t,dB_inv(T_record[:,lowlevel_arg]),dB_inv(T_record[:,highlevel_arg]))
    Tranmission_penalty = TP(dB_inv(T_record[:,highlevel_arg]),dB_inv(T_record[:,lowlevel_arg]),1)
    Insertion_Loss = IL(t,dB_inv(T_record[:,highlevel_arg]),dB_inv(T_record[:,lowlevel_arg]))
    ploting(wl*1000,Extinction_ratio,Tranmission_penalty,Insertion_Loss,x_label="wavelength(nm)",\
            title="vswing = "+str(vbias[highlevel_arg])+"~"+str(vbias[lowlevel_arg]),filename="Transmission ER TP",leg=["ER","TP","IL"])
    

    # v.v_bias=-1.5
    # b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    # T = Transfer_function(ring_mod,t)
    # wl,data_v1 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    # wl,data_phase_v1 = T.mapping(180/np.pi*np.angle(s_minus))
    # v.v_bias=0
    # b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    # T = Transfer_function(ring_mod,t)
    # wl,data_v2 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    # wl,data_phase_v2 = T.mapping(180/np.pi*np.angle(s_minus))
    # v.v_bias=-0.5
    # b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    # T = Transfer_function(ring_mod,t)
    # wl,data_v3 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    # wl,data_phase_v3 = T.mapping(180/np.pi*np.angle(s_minus))

    V = np.linspace(-5,0,1000)
    ploting(V,ring_mod.alpha(V),x_label="voltage (V)",title="Amplitude absorption coefficient (1/cm)",filename="alpha_V")
    ploting(V,ring_mod.neff(V),x_label="voltage (V)",title="neff vs Voltage (1/cm)",filename="neff_V")
if sim.mode == "voltage_drive":

    t.main(ring_mod,N=bit_num,driver=v)
    print("\n\nSimulation at ",str(v.f_drive/1e9)," GHz, ",str(v.vpp),"V vpp, ",str(v.v_bias),"V vbias\n\n")
    filename = ('sim_'+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV"+"_vbias_"+str(v_bias)+str(v.level))
    sim.save_data(ring_mod,t,v,file_name=filename)
    v.method = "large_signal"
    sim.eye_diagram(t,v,v.v,
                    filename="voltage_eye_"+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp))+"_vbias_"+str(v_bias),
                    title="vpos" ,
                    plot_bit_num=2)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    print("\n\n.................Start Large Signal Simulation.................\n\n")
    b,Q,s_minus,N ,vneg, vj, i2= solving(sim,ring_mod,v,t,H)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    name = 'eye_LargeSignal_'+str(bit_num)+'bits_with_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp*1000))+"mV_vbias_"+str(v.v_bias)
    sim.eye_diagram(t,v,abs(s_minus)**2,filename=name+str(v.level),plot_bit_num=2)
    sim.eye_diagram(t,v,vneg,filename="vneg LargeSignal"+str(v.level),plot_bit_num=2,title="vneg")
    sim.eye_diagram(t,v,vj,filename="vj LargeSignal"+str(v.level),plot_bit_num=2,title="vj")

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # print("\n\n.................Start Small Signal Simulation.................\n\n")
    # v.method = "small_signal"
    # b1,Q1,s_minus1,N1,vneg, vj, i2 = solving(sim,ring_mod,v,t,H)

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # ploting(t.t_total,v.v,v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='voltage (V)',filename='voltage_and_Large_signal_vj',leg=['External Voltage','Large Signal Vj'])
    # ploting(t.t_total,v.v,Q1/v.cj_normalizing,x_label='time (ps)',title='voltage (V)',filename='voltage_and_small_signal_vj',leg=['External Voltage','Small Signal Vj'])
    # ploting(t.t_total,v.v-v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='Resistor Voltage (V)',filename='Rv')
    # # print("len of t_total = ",len(t.t_total))
    # # print("len of Q = ",len(Q))
    # # print("len of Q1 = ",len(Q1))
    # ploting(t.t_total,Q,Q1,x_label='time (ps)',title='Q',filename='Q',leg=['Large Signal','Small Signal'])
    # ploting(t.t_total,v.V_Q(Q/v.cj_normalizing),Q1/v.cj_normalizing,x_label='time (ps)',title='junction voltage',filename='junction voltage',leg=['large signal','small signal'])
    # ploting(t.t_total,abs(s_minus)**2,abs(s_minus1)**2,x_label='time (ps)',title='output Power',filename='output Power',leg=['large signal','small signal'])

    # # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # name = 'eye_SmallSignal_'+str(bit_num)+'bits_with_TPA_FCA_f'+   \
    #     str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp))+"_vbias_"+str(v.v_bias)
    # sim.eye_diagram(t,v,abs(s_minus1)**2,filename=name+str(v.level),plot_bit_num=2)
    # sim.eye_diagram(t,v,vneg,filename="vneg SmallSignal"+str(v.level),plot_bit_num=2,title="vneg")
    # sim.eye_diagram(t,v,vj,filename="vj SmallSignal"+str(v.level),plot_bit_num=2,title="vj")