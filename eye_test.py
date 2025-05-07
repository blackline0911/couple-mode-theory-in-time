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
mode = "scan_frequency"
# mode = "voltage_drive"
wl_in = 1.5467
Pin = 1 #mW
FSR = 0.0195
radius = 5
RoundTripLoss_pdk = [0.95124, 0.95248, 0.95273, 0.95292, 0.95305]
me_data = [39.5, 37.9]

bit_num = 300
v_bias = -1
vpp = 1
Rs = 53.9
Cjs = [23.6e-15, 20e-15]
f_drive=50

L = 2*np.pi*radius*1e-4


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////

a_fit = alpha_fit(RoundTripLoss=RoundTripLoss_pdk,L = 2*np.pi*radius)
m_fit = me_fit(me_data=me_data)
ring_mod = ring(L=2*np.pi*radius, 
            alpha=a_fit.alpha_V,
            me=m_fit.me_V,
            cross_section=0.5*0.22,
            lambda_incident=wl_in,
            gamma=[0.95105],
            FSR = FSR,
            neff=2.5111,
            )

wl_min =  ring_mod.lambda0 - ring_mod.lambda0/ring_mod.Q/8
wl_max =  ring_mod.lambda0 + ring_mod.lambda0/ring_mod.Q/8 

v = driver(f_drive=f_drive,
           v_bias=v_bias,
           vpp=vpp,
           Rs=Rs,
           raise_cosine=1,
           cj = Cjs,
           PRBS=1)
H = Heater(300,0,0.5*150/0.6)
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
os.chdir("./eye_diagram_test/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    vbias = np.array([0,-1,-2,-3])
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=5000,resolution=1,buffer=100,driver=v)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    sim.save_data(ring_mod,t,v)
    plt.figure()
    for vb in vbias:
        v.v_bias = vb
        b,Q,s_minus,N = solving(sim,ring_mod,v,t,H)
        T = Transfer_function(ring_mod,t)
        wl,data = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
        wl,data_phase = T.mapping(180/np.pi*np.angle(s_minus))
        plt.plot(wl,data,
                 label="V = "+str(vb))
    plt.grid(color='g',linestyle='--', alpha=0.5)
    plt.xlabel('wavelength(um)')
    plt.title('Transfer function')
    plt.legend()
    plt.savefig("Transmission_vs_voltage")
    plt.show()
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
V = np.linspace(-2,0,1000)
ploting(V,ring_mod.alpha(V),x_label="voltage (V)",title="absorption coefficient (1/cm)",filename="alpha_V")
if sim.mode == "voltage_drive":
    
    # F = np.arange(50,51)
    # VPP = np.arange(1,2)
    # VBIAS = np.array([-1])
    # for f in F:
    #     for vpp in VPP:
    #         for v_bias in VBIAS:
    f = 50
    vpp=6
    v_bias=-2
    v.v_bias = v_bias
    v.vpp = vpp
    v.f_drive = f*1e9
    v.renew()
    t.main(ring_mod,N=bit_num,driver=v)
    print("\n\nSimulation at ",str(v.f_drive/1e9)," GHz, ",str(v.vpp),"V vpp, ",str(v.v_bias),"V vbias\n\n")
    filename = ('sim_'+str(int(f))+"GHz_vpp_"+str(int(vpp))+"_vbias_"+str(v_bias))
    sim.save_data(ring_mod,t,v,file_name=filename)
    v.method = "large_signal"
    sim.eye_diagram(t,v,v.v,
                    filename="voltage_eye_"+str(int(f))+"GHz_vpp_"+str(int(vpp))+"_vbias_"+str(v_bias), 
                    plot_bit_num=2)
    
    print("\n\n.................Start Large Signal Simulation.................\n\n")
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    

    name = 'eye_LargeSignal_'+str(bit_num)+'bits_without_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp))+"_vbias_"+str(v.v_bias)
    sim.eye_diagram(t,v,abs(s_minus)**2,filename=name,plot_bit_num=2)

    print("\n\n.................Start Small Signal Simulation.................\n\n")
    v.method = "small_signal"
    b1,Q1,s_minus1,N1 = solving(sim,ring_mod,v,t)
    ploting(t.t_total,v.v,v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='voltage (V)',filename='voltage_and_Large_signal_vj',leg=['External Voltage','Large Signal Vj'])
    ploting(t.t_total,v.v,Q1/v.cj_normalizing,x_label='time (ps)',title='voltage (V)',filename='voltage_and_small_signal_vj',leg=['External Voltage','Small Signal Vj'])
    ploting(t.t_total,v.v-v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='Resistor Voltage (V)',filename='Rv')
    # print("len of t_total = ",len(t.t_total))
    # print("len of Q = ",len(Q))
    # print("len of Q1 = ",len(Q1))
    ploting(t.t_total,Q,Q1,x_label='time (ps)',title='Q',filename='Q',leg=['Large Signal','Small Signal'])
    ploting(t.t_total,v.V_Q(Q/v.cj_normalizing),Q1/v.cj_normalizing,x_label='time (ps)',title='junction voltage',filename='junction voltage',leg=['large signal','small signal'])
    ploting(t.t_total,abs(s_minus)**2,abs(s_minus1)**2,x_label='time (ps)',title='output Power',filename='output Power',leg=['large signal','small signal'])

    name = 'eye_SmallSignal_'+str(bit_num)+'bits_without_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp))+"_vbias_"+str(v.v_bias)
    sim.eye_diagram(t,v,abs(s_minus1)**2,filename=name,plot_bit_num=2)