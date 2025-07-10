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


# Refer：Ring modulator in IMEC PDK，以PDK 3.6.1為例
# 此模擬檔案為自行設計ring modulator，以PDK 3.6.1的ring modulator為參考
# ring內的round trip loss、調變效率訂為跟PDK 3.6.1的ring modulator相同，然後調整ring radius和couple coefficient，為了有更高Q和更深deep
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
mode = "scan_frequency"
# mode = "voltage_drive"
delta_f = 8 # GHz
wl_in = 1.559453431986187
# wl_in = 1.5593-0.0008/100*delta_f
Pin = 1.0 #mW
ng = 3.9697242191825954
radius = 10 #um
cavity_length = 2*np.pi*radius 
La_Lc_ratio = 2*np.pi*5/(2*np.pi*radius)
mode_area = 0.22*0.5
lambda_res =  1.559453431986187
me = 38 # pm/V

x = dB(0.994898) # energy round trip loss in dB
y = x*4*((radius-5)/radius)
z = 10**(y/10) 
print("x = ",x)
print("z = ",z)
print("total round trip loss = ",)
Amp_RoundTripLoss_data = np.array([0.95223,0.95248,0.95273,0.95292,0.95305])*(z)**0.5
print("Amp_RoundTripLoss_data = ",Amp_RoundTripLoss_data)
a_fit = alpha_fit(RoundTripLoss=Amp_RoundTripLoss_data,La=cavity_length*La_Lc_ratio,L = cavity_length,input2 = "amp")

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the coupling coefficient
# gamma = 0.95105
# gamma = 0.94898839
gamma = 0.94737172

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the index (varing with voltage)
w_res = 2*np.pi*c/lambda_res
D_bar = -2*np.pi*c/lambda_res**2
print("delta neff/delta V = ",me*1e-6*D_bar*( -ng/(w_res) )*(1/La_Lc_ratio))
neff0 = 2.680503
dneff_dV = me*1e-6*D_bar*( -ng/(w_res) )*(1/La_Lc_ratio)
neff_calculated = [neff0+dneff_dV*(-0.5),neff0, neff0+dneff_dV*0.5, neff0+dneff_dV*1, neff0+dneff_dV*1.5, neff0+dneff_dV*2]

n_fit = neff_fit(neff_data=neff_calculated)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the equivalent circuit parameters
# These parameters are based on the paper's data
# The values are chosen to match the simulation conditions in the paper
# The values can be adjusted based on the specific device characteristics
bit_num = 1000
v_bias = -1/2
vpp = 1/2
Rs =53.9
a_cj = 37.5e-15
b_cj = a_cj**2/(20e-15)**2 - 1
Cjs = np.real([a_cj/(b_cj)**0.5, a_cj/(b_cj + 1)**0.5])
f_drive= 70
level = "PAM4"
Cox = 34.7e-15
Rsi = 1439.0
Cpad = 6.6e-15


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////

ring_mod = ring(L=cavity_length, 
            L_active = [cavity_length*La_Lc_ratio],
            alpha=a_fit.alpha_V,
            neff=n_fit.neff_V,
            cross_section=mode_area,
            lambda_incident=wl_in,
            gamma=[gamma],
            ng=ng,
            FSR_shift=0,
            FCA_fit_factor=0,
            TPA_fit_factor=0,
            SPM_fit_factor=0,
            self_heating_factor=0,
            band="C",
            HE = 73,
            )
print("FSR = ",(ring_mod.lambda0*1000)**2/(ng*ring_mod.L*1000)," nm")
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the heater parameters
H = Heater(300,0,0.5*150/0.6)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the wavelength range for the scan
# The range is set based on the resonant wavelength and the heater power
# The values can be adjusted based on the specific device characteristics
wl_min =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 - ring_mod.lambda0/ring_mod.Q
wl_max =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 + ring_mod.lambda0/ring_mod.Q

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
           beta=0.6,
           )
V = np.linspace(-2,0,1000)
ploting(V,v.Cj_V(V,v.a,v.b),x_label="voltage (V)",title="junction capacitance (F)",filename="Cj_V")

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

V = np.linspace(-2,0,1000)
ploting(V,ring_mod.alpha(V),x_label="voltage (V)",title="Energy absorption coefficient (1/cm)",filename="alpha_V")
ploting(V,ring_mod.neff(V),x_label="voltage (V)",title="neff vs Voltage",filename="neff_V")
plt.scatter(np.array([0.5,0,-0.5,-1,-1.5,-2]), neff_calculated)
plt.show()
ploting(V,1e6*c*1e-12/sim.f_pround_bar/ring_mod.ng*( ring_mod.neff(V) - ring_mod.neff(0))*(ring_mod.L_active[0]/ring_mod.L),\
            x_label="voltage (V)",title="resonant wavelength vs Voltage (pm/V)",filename="lambda_V")

os.chdir("./eye_test_design/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    # vbias = np.array([v_bias+vpp/2,v_bias,v_bias-vpp/2])
    vbias = np.array([0,-0.5,-1,-1.5,-2])
    # vbias = np.arange(-0,-0.5,-0.5)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=10000,resolution=1,buffer=50,driver=v)

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
    # ploting(t.t_total,abs(sim.b)**2,x_label="time (ps)",title="Energy in Ring (mJ)")

    sim.save_data(ring_mod,t,v,H)
    
if sim.mode == "voltage_drive":
    t.main(ring_mod,N=bit_num,driver=v,resolution=1)
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
    sim.eye_diagram(t,v,abs(s_minus)**2,filename=name+str(v.level),plot_bit_num=2,title="Output power (mW)")
    sim.eye_diagram(t,v,vneg,filename="vneg_SmallSignal"+str(v.level),plot_bit_num=2,title="vneg")
    sim.eye_diagram(t,v,vj,filename="vj_SmallSignal"+str(v.level),plot_bit_num=2,title="vj")
    sim.eye_diagram(t,v,v.v+vneg,filename="vpos + vneg SmallSignal"+str(v.level),plot_bit_num=2,title="vpos_plus_vneg")

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
