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


# Refer：Ultra-Wide Free-Spectral-Range Silicon Microring Modulator for High Capacity WDM
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# mode = "scan_frequency"
mode = "voltage_drive"
wl_in = 1.54969
Pin = 1 #mW
FSR = 0.0324
radius = 3
La_Lc_ratio = 1
dneff_dV = 8.248835006293938e-05

G_plus_alpha = 151.87144816774455
ratio = 165/(192*2)
alpha_energy0 =  G_plus_alpha*(1-ratio) #1/cm
G_energy = G_plus_alpha*ratio #1/cm

dalpha_dv = 0.5
def func(V,alpha0,b,c):
    return b*V/(abs(V)+c)**0.5+alpha0
V = np.array([0,-0.5,-1,-1.5,-2])
alpha_energy_data = func(V,alpha_energy0,1.4,2)
Amp_RoundTripLoss_data = np.exp(-alpha_energy_data/2*2*np.pi*radius*1e-4)

gamma = np.exp(-G_energy/2*2*np.pi*radius*1e-4)
print("gamma = ",gamma)
neff0 = 2.631361725786003
neff_calculated = [neff0+dneff_dV*(-0.5),neff0, neff0+dneff_dV*0.5, neff0+dneff_dV*1, neff0+dneff_dV*1.5, neff0+dneff_dV*2]
print("neff_calculated = ",neff_calculated)
mode_area = 0.22*0.5


bit_num = 500
v_bias = -1.5
vpp = 2
Rs =250
a_cj = 20e-15
b_cj = a_cj**2/(6e-15)**2 - 3
Cjs = [a_cj/(b_cj)**0.5, a_cj/(b_cj + 1)**0.5]
f_drive= 100
level = "PAM4"
Cox = 19.2e-15
Rsi = 760.0
Cpad = 20.1e-15


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////

a_fit = alpha_fit(RoundTripLoss=Amp_RoundTripLoss_data,L = 2*np.pi*radius,input2 = "amp",fit_mode="func")
n_fit = neff_fit(neff_data=neff_calculated)
ring_mod = ring(L=2*np.pi*radius, 
            L_active = [2*np.pi*radius*La_Lc_ratio],
            alpha=a_fit.alpha_V,
            neff=n_fit.neff_V,
            cross_section=mode_area,
            lambda_incident=wl_in,
            gamma=[gamma],
            FSR = FSR,
            FSR_shift=0,
            FCA_fit_factor=0,
            TPA_fit_factor=0,
            SPM_fit_factor=0,
            self_heating_factor=0,
            band="C",
            HE = 73
            )



H = Heater(300,0,0.5*150/0.6)
wl_min =  1.548
wl_max =  1.552
# wl_min =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 - ring_mod.lambda0/ring_mod.Q/1.5
# wl_max =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 + ring_mod.lambda0/ring_mod.Q/2

v = driver(f_drive=f_drive,
           v_bias=v_bias,
           vpp=vpp,
           Rs=[Rs],
           raise_cosine=1,
           cj = [Cjs],
           PRBS=1,
           level = level,
           Cox = [Cox],
           Rsi=[Rsi],
           Cp=[Cpad],
           )
V = np.linspace(-5,0,1000)
ploting(V,v.Cj_V(V,v.a,v.b),x_label="voltage (V)",title="junction capacitance (F)",filename="Cj_V")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
os.chdir("./eye_diagram_test_v3/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    print("wl_min = ",wl_min ," um")
    print("wl_max = ",wl_max ," um")
    # vbias = np.array([v_bias+vpp/2,v_bias,v_bia…s-vpp/2])
    vbias = np.array([0,-1,-2,-3,-4])
    # vbias = np.arange(-3,-3.5,-0.5)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=10000,resolution=0,buffer=50,driver=v)
    print("LineWidth = ",ring_mod.lambda0/ring_mod.Q*1000," nm")
    print(H.P*ring_mod.HE/1e6)
    print("detuning wabvelength = ",(ring_mod.lambda0+H.P*ring_mod.HE/1e6-wl_in)*1000," nm")
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    
    T_record = np.zeros( (int(len(t.t_total)-t.buffer*t0/t.dt-1),len(vbias)))
    plt.figure()
    for vb in vbias:
        v.v_bias = vb
        sim.b,s_minus = solving(sim,ring_mod,v,t,H)
        T = Transfer_function(ring_mod,t)
        wl,data = T.mapping(dB(abs(s_minus)**2/sim.Pin))
        wl,data_phase = T.mapping(180/np.pi*np.angle(s_minus))
        plt.plot(wl*1000,data,
                 label="V = "+str(vb))
    plt.grid(color='g',linestyle='--', alpha=0.5)
    plt.xlabel('wavelength(nm)')
   
    plt.legend()
    if ring_mod.FCA_fit_factor==1:
        plt.title('Transfer function (with NL absorb) (func fit alpha)')
        plt.savefig("Transmission_vs_voltage (with NL) (func fit alpha)")
    else:
        plt.title('Transfer function (no NL absorb) (func fit alpha)')
        plt.savefig("Transmission_vs_voltage (no NL) (func fit alpha)")
    plt.show()
    ploting(t.t_total,abs(sim.b)**2,x_label="time (ps)",title="Energy in Ring (mJ)")

    sim.save_data(ring_mod,t,v)

    V = np.linspace(-5,0,1000)
    
    ploting(V,ring_mod.alpha(V),x_label="voltage (V)",title="Energy absorption coefficient (1/cm)",filename="alpha_V")
    ploting(V,ring_mod.neff(V),x_label="voltage (V)",title="neff vs Voltage",filename="neff_V")
    ploting(V,1e6*c*1e-12/sim.f_pround_bar/ring_mod.ng*( ring_mod.neff(V) - ring_mod.neff(0))*(ring_mod.L_active[0]/ring_mod.L),\
            x_label="voltage (V)",title="resonant wavelength vs Voltage (pm/V)",filename="lambda_V")
if sim.mode == "voltage_drive":
    t.main(ring_mod,N=bit_num,driver=v)
    print("\n\nSimulation at ",str(v.f_drive/1e9)," GHz, ",str(v.vpp),"V vpp, ",str(v.v_bias),"V vbias\n\n")
    filename = ('sim_'+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV"+"_vbias_"+str(v_bias)+str(v.level))
    
    # v.method = "large_signal"
    sim.eye_diagram(t,v,v.v,
                    filename="voltage_eye_"+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV_vbias_"+str(int(1000*v_bias))+"mV",
                    title="vpos" ,
                    plot_bit_num=2)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    print("\n\n.................Start Small Signal Simulation.................\n\n")
    sim.b,s_minus,N ,delta_T,vneg, vj, i2= solving(sim,ring_mod,v,t,H)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    sim.save_data(ring_mod,t,v,file_name=filename)
    name = 'eye_SmallSignal_'+str(bit_num)+'bits_with_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(-int(v.vpp*1000))+"mV_vbias_"+str(-int(v.v_bias*1000))+"mV"
    ploting(t.t_total,delta_T,x_label="time",title="delta T")
    sim.eye_diagram(t,v,abs(s_minus)**2,filename=name+str(v.level),plot_bit_num=4,title="Output power (mW)")
    sim.eye_diagram(t,v,vneg,filename="vneg_SmallSignal"+str(v.level),plot_bit_num=4,title="vneg")
    sim.eye_diagram(t,v,vj,filename="vj_SmallSignal"+str(v.level),plot_bit_num=4,title="vj")
    sim.eye_diagram(t,v,v.v+vneg,filename="vpos + vneg SmallSignal"+str(v.level),plot_bit_num=4,title="vpos_plus_vneg")
