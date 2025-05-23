import numpy as np
import matplotlib.pyplot as plt
from utility import *
# Refer to imec paper : 
# Electro-Opto-Thermal Dynamic Compact Model for Ring Modulator Arrays
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


# Refer：O-Band Silicon Ring Modulators With Highly Efficient Electro- and Thermo-Optic Modulation
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
mode = "scan_frequency"
# mode = "voltage_drive"
wl_in = 1.3054
Pin = 1 #mW
FSR = 0.0135
radius = 5
La_Lc_ratio = 1
dneff_dV = 0.00010048812689651478

G_plus_alpha = 38.680029629629686
ratio = 1-16/48
alpha_energy0 =  G_plus_alpha*(1-ratio) #1/cm
G_energy = G_plus_alpha*ratio #1/cm
print("alpha_energy0 = ",alpha_energy0, "1/cm")
print("G_energy = ",G_energy, "1/cm")

dalpha_dv = 0.5
def func(V,alpha0,b,c):
    return b*V/(abs(V)+c)**0.5+alpha0
V = np.array([0,-0.5,-1,-1.5,-2])
alpha_energy_data = func(V,alpha_energy0,0.5,2e-06)
Amp_RoundTripLoss_data = np.exp(-alpha_energy_data/2*2*np.pi*radius*1e-4)
print(alpha_energy_data)


# gamma = 0.947582
gamma = np.exp(-G_energy/2*2*np.pi*radius*1e-4)
print("gamma(amplitude) = ",gamma)
neff0 = 2.700996734985352
neff_calculated = [neff0+dneff_dV*(-0.5),neff0, neff0+dneff_dV*0.5, neff0+dneff_dV*1, neff0+dneff_dV*1.5, neff0+dneff_dV*2]
print("neff_calculated = ",neff_calculated)
mode_area = 0.22*0.5


bit_num = 50
v_bias = -0.3
vpp = 1.6
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

a_fit = alpha_fit(RoundTripLoss=Amp_RoundTripLoss_data,L = 2*np.pi*radius,input = "amp",fit_mode="func")
n_fit = neff_fit(neff_data=neff_calculated)
ring_mod = ring(L=2*np.pi*radius, 
            L_active = 2*np.pi*radius*La_Lc_ratio,
            alpha=a_fit.alpha_V,
            neff=n_fit.neff_V,
            cross_section=mode_area,
            lambda_incident=wl_in,
            gamma=[gamma],
            FSR = FSR,
            FSR_shift=0,
            FCA_fit_factor=1,
            TPA_fit_factor=1,
            SPM_fit_factor=1,
            band="O",
            HE = 222
            )
print("G_energy + alpha_energy = ",(2*np.pi*229646656979082.34*ring_mod.ng/(c*1e-4))/5000," 1/cm")
# print("delta neff/delta V = ",11e-6*ring_mod.D_bar*( -ring_mod.ng/(2*np.pi*ring_mod.f_res_bar) )*(1/La_Lc_ratio))
print("ng = ",ring_mod.ng)
# print("lambda0 = ",ring_mod.lambda0)
# print("f0 = ",ring_mod.f_res_bar," THz")

# f0 =  229.08878656875603  THz
# ng =  3.984719655112499
# lambda0 =  1.3086299966498962 um

# wl_min =  1.5464
# wl_max =  1.55
# H = Heater(300,2.42,0.5*150/0.6)
H = Heater(300,0,0.5*2*np.pi*5/0.6)
wl_min =  1.3052
wl_max =  1.3058
# wl_min =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 - ring_mod.lambda0/ring_mod.Q/1.5
# wl_max =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 + ring_mod.lambda0/ring_mod.Q/2
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
os.chdir("./thermal_model_test_v2/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    vbias = np.array([v_bias+vpp/2,v_bias,v_bias-vpp/2])
    # vbias = np.array([0,-1,-2,-3,-4])
    vbias = np.arange(0,-0.5,-0.5)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=10000,resolution=2,buffer=100,driver=v)
    print("LineWidth = ",ring_mod.lambda0/ring_mod.Q*1000," nm")
    print("Thermal tuning wavelength = ",H.P*ring_mod.HE/1e6," um")
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
        plt.plot(wl*1000,data,
                 label="V = "+str(vb))
        T_record[:,int(np.argwhere(vbias==vb))] = data
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

    
    # highlevel_arg = int(np.argwhere(vbias==v_bias+vpp/2))
    # lowlevel_arg = int(np.argwhere(vbias==v_bias-vpp/2))
    # Extinction_ratio = ER(t,dB_inv(T_record[:,lowlevel_arg]),dB_inv(T_record[:,highlevel_arg]))
    # Tranmission_penalty = TP(dB_inv(T_record[:,highlevel_arg]),dB_inv(T_record[:,lowlevel_arg]),1)
    # Insertion_Loss = IL(t,dB_inv(T_record[:,highlevel_arg]),dB_inv(T_record[:,lowlevel_arg]))
    # ploting(wl*1000,Extinction_ratio,Tranmission_penalty,Insertion_Loss,x_label="wavelength(nm)",\
    #         title="vswing = "+str(vbias[highlevel_arg])+"~"+str(vbias[lowlevel_arg]),filename="Transmission ER TP",leg=["ER","TP","IL"])
    
    sim.save_data(ring_mod,t,v)

    V = np.linspace(-5,0,1000)
    ploting(V,ring_mod.alpha(V),x_label="voltage (V)",title="Energy absorption coefficient (1/cm)",filename="alpha_V")
    ploting(V,ring_mod.neff(V),x_label="voltage (V)",title="neff vs Voltage",filename="neff_V")
    ploting(V,1e6*c*1e-12/sim.f_pround_bar/ring_mod.ng*( ring_mod.neff(V) - ring_mod.neff(0))*(ring_mod.L_active/ring_mod.L),\
            x_label="voltage (V)",title="resonant wavelength vs Voltage (pm/V)",filename="lambda_V")
if sim.mode == "voltage_drive":
    t.main(ring_mod,N=bit_num,driver=v)
    print("\n\nSimulation at ",str(v.f_drive/1e9)," GHz, ",str(v.vpp),"V vpp, ",str(v.v_bias),"V vbias\n\n")
    filename = ('sim_'+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV"+"_vbias_"+str(v_bias)+str(v.level))
    sim.save_data(ring_mod,t,v,file_name=filename)
    sim.eye_diagram(t,v,v.v,
                    filename="voltage_eye_"+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp))+"_vbias_"+str(v_bias),
                    title="vpos" ,
                    plot_bit_num=2)
    def CMT_thermal(t_bar,eqs,driver:driver,ring:ring,Heater:Heater,SPM=None,TPA=None,FCA=None,sim=None,T_args=None):
        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        b_bar , Q_pround ,N_bar,delta_T, v_neg, vj, i2= eqs
        voltage = np.real(driver.refering_v(t_bar))
        dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
        alpha_linear = ring.alpha(vj)
        cj = driver.Cj_V(vj)
        dlambda = ring.lambda0/ring.ng*( ring.neff(vj) - ring.neff(0))
        tau_heater = 5e-6

        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,delta_T,ring.f_res_bar,alpha_linear,TPA,SPM,T_args,dlambda,Heater)
        
        dQ_dt = (voltage/(driver.Rs * driver.cj_normalizing )*t0) \
            - (1/( driver.Rs ) )*driver.V_Q(Q_pround)*t0/driver.cj_normalizing
        
        dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)

        dT_dt = Heater.T_rate_equation(b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim)
        return [ da_dt,dQ_dt ,dN_dt, dT_dt, dvneg_dt, dvj_dt, di2_dt]

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    b,Q,s_minus,N ,delta_T,vneg, vj, i2= solving(sim,ring_mod,v,t,H)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    name = 'eye_LargeSignal_'+str(bit_num)+'bits_with_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp*1000))+"mV_vbias_"+str(v.v_bias)
    ploting(t.t_total,delta_T,x_label="time",title="delta T")
    sim.eye_diagram(t,v,abs(s_minus)**2,filename=name+str(v.level),plot_bit_num=2)
    sim.eye_diagram(t,v,vneg,filename="vneg LargeSignal"+str(v.level),plot_bit_num=2,title="vneg")
    sim.eye_diagram(t,v,vj,filename="vj LargeSignal"+str(v.level),plot_bit_num=2,title="vj")

    