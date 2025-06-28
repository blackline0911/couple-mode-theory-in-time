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


# Refer：Efficient 330-Gb/s PAM-8 modulation using silicon microring modulators
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

mode = "scan_frequency"
# mode = "voltage_drive"
wl_in = 1.555
# Pin = 3.68 #mW
Pin = 1 #mW
FSR = 0.012667
radius = 7.5
La_LSB_ratio = 1/3
La_MSB_ratio = 2/3
segment = 1
SB = "LSB"
mode_area = 0.22*0.5
lambda_res = 1.555
Q = 5200
me_LSB = 30 # pm/V

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Cox_LSB = 57.6e-15
Rs_LSB = 109
Rsi_LSB = 416
Cpad_LSB = 17.6e-15
bit_num = 1000
v_bias = -2.5/2
vpp = 1.8/2

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Calculate group index
cavity_length = 2*np.pi*radius
ng = lambda_res**2/(FSR*cavity_length)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Calculate total loss (intrinsic loss + coupling loss)
w_res = 2*np.pi*c/(lambda_res) # THz
print("G_energy + alpha_energy = ",( w_res*ng/(c*1e-4))/Q," 1/cm")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the energy absorption coefficient (varing with voltage)

G_plus_alpha = ( w_res*ng/(c*1e-4))/Q
ratio = 1-181/(192*2)
alpha_energy0 =  G_plus_alpha*(1-ratio) #1/cm
G_energy = G_plus_alpha*ratio #1/cm

def func(V,alpha0,b,c):
    return b*V/(abs(V)+c)**0.5+alpha0
V = np.array([0,-0.5,-1,-1.5,-2])
# d = 0.23355001135484965
d = 0.28
# alpha_energy_data = func(V,alpha_energy0,0.6,0.6)
alpha_energy_data = V*d + alpha_energy0
Amp_RoundTripLoss_data = np.exp(-alpha_energy_data/2*2*np.pi*radius*1e-4)
a_fit = alpha_fit(RoundTripLoss=Amp_RoundTripLoss_data,L = 2*np.pi*radius,input2 = "amp",fit_mode="linear")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the coupling coefficient
gamma = np.exp(-G_energy/2*2*np.pi*radius*1e-4)


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Define the index (varing with voltage)
neff0 = 2.639849989417571
D_bar = -2*np.pi*c/lambda_res**2
dneff_dV = me_LSB*1e-6*D_bar*( -ng/(w_res) )*(1/La_LSB_ratio)
print("delta neff/delta V = ",me_LSB*1e-6*D_bar*( -ng/(w_res) )*(1/La_LSB_ratio))

neff_calculated = [neff0+dneff_dV*(-0.5),neff0, neff0+dneff_dV*0.5, neff0+dneff_dV*1, neff0+dneff_dV*1.5, neff0+dneff_dV*2]
# print("neff_calculated = ",neff_calculated)
n_fit = neff_fit(neff_data=neff_calculated)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



a_cj = 60e-15
b_cj = a_cj**2/(27e-15)**2 - 3
Cjs_LSB = [a_cj/(b_cj)**0.5, a_cj/(b_cj + 1)**0.5]
a_cj = 100e-15
b_cj = a_cj**2/(8.5e-15)**2 - 2.5
Cjs_MSB = [a_cj/(b_cj)**0.5, a_cj/(b_cj + 1)**0.5]
f_drive=100
level = "PAM4"
print("Cjs_LSB = ",Cjs_LSB)
print("Cjs_MSB = ",Cjs_MSB)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if segment==2:
    Cox = [Cox_LSB,Cox_MSB]
    Rsi = [Rsi_LSB,Rsi_MSB]
    Cpad = [Cpad_LSB,Cpad_MSB]
    Rs = [Rs_LSB,Rs_MSB]
    Cj = [Cjs_LSB,Cjs_MSB]
    La = [2*np.pi*radius*La_LSB_ratio,2*np.pi*radius*La_MSB_ratio]
if segment==1:
    if SB=="LSB":
        Cox = [Cox_LSB]
        Rsi = [Rsi_LSB]
        Cpad = [Cpad_LSB]
        Rs = [Rs_LSB]
        Cj = [Cjs_LSB]
        La = [2*np.pi*radius*La_LSB_ratio]
    if SB=="MSB":
        Cox = [Cox_MSB]
        Rsi = [Rsi_MSB]
        Cpad = [Cpad_MSB]
        Rs = [Rs_MSB]
        Cj = [Cjs_MSB]
        La = [2*np.pi*radius*La_MSB_ratio]



# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":Pin} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////


ring_mod = ring(L=2*np.pi*radius, 
            L_active = La,
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
            band="C",
            HE = 73,
            segment=segment,
            self_heating_factor=0
            )
ploting(V,1e6*c*1e-12/sim.f_pround_bar/ring_mod.ng*( ring_mod.neff(V) - ring_mod.neff(0))*(ring_mod.L_active[0]/ring_mod.L),\
            x_label="voltage (V)",title="resonant wavelength vs Voltage (pm/V)",filename="lambda_V")

H = Heater(300,0,0.5*150/0.6)
wl_min =  1.5545
wl_max =  1.5555
# wl_min =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 - ring_mod.lambda0/ring_mod.Q
# wl_max =  ring_mod.lambda0+ring_mod.HE*H.P*1e-6 + ring_mod.lambda0/ring_mod.Q
print("wl_min = ",wl_min ," um")
print("wl_max = ",wl_max ," um")

v = driver(f_drive=f_drive,
           v_bias=v_bias,
           vpp=vpp,
           Rs=Rs,
           raise_cosine=1,
           cj = Cj,
           level = level,
           Cox = Cox,
           Rsi=Rsi,
           Cp=Cpad,
           segment=segment)

V = np.linspace(-5,0,1000)
if segment==2:
    ploting(V,v.Cj_V(V,v.a,v.b),v.Cj_V(V,v.c,v.d),x_label="voltage (V)",title="junction capacitance LSB (F)",filename="Cj_V_LSB",leg=["LSB","MSB"])
if segment==1:
    if SB=="LSB":
        ploting(V,v.Cj_V(V,v.a,v.b),x_label="voltage (V)",title="junction capacitance LSB (F)",filename="Cj_V_LSB")
    if SB=="MSB":
        ploting(V,v.Cj_V(V,v.a,v.b),x_label="voltage (V)",title="junction capacitance MSB (F)",filename="Cj_V_MSB")
        
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
os.chdir("./eye_diagram_test_v6/")
t = time(mode = sim.mode)

if sim.mode == "scan_frequency":
    # vbias = np.array([v_bias+vpp/2,v_bias,v_bias-vpp/2])
    vbias = np.array([0,-1,-2,-3,-4,-5])
    # vbias = np.arange(0,-0.5,-0.5)
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=15000,resolution=0,buffer=200,driver=v)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    
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
    # ploting(V,1e6*c*1e-12/sim.f_pround_bar/ring_mod.ng*( ring_mod.neff(V) - ring_mod.neff(0))*(ring_mod.L_active/ring_mod.L),\
    #         x_label="voltage (V)",title="resonant wavelength vs Voltage (pm/V)",filename="lambda_V")
if sim.mode == "voltage_drive":
    t.main(ring_mod,N=bit_num,driver=v,resolution=1)
    print("\n\nSimulation at ",str(v.f_drive/1e9)," GHz, ",str(v.vpp),"V vpp, ",str(v.v_bias),"V vbias\n\n")
    filename = ('sim_'+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV"+"_vbias_"+str(v_bias)+str(v.level))
    
    # v.method = "large_signal"
    if segment==2:
        sim.eye_diagram(t,v,v.v_LSB,
                        filename="voltage_LSB_eye_"+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV_vbias_"+str(int(1000*v_bias))+"mV",
                        title="vpos LSB" ,
                        plot_bit_num=2)
        sim.eye_diagram(t,v,v.v_MSB,
                        filename="voltage_LSB_eye_"+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV_vbias_"+str(int(1000*v_bias))+"mV",
                        title="vpos MSB" ,
                        plot_bit_num=2)
    if segment==1:
        sim.eye_diagram(t,v,v.v,
                        filename="voltage_eye_"+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp*1000))+"mV_vbias_"+str(int(1000*v_bias))+"mV",
                        title="vpos" ,
                        plot_bit_num=2)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    print("\n\n.................Start Small Signal Simulation.................\n\n")
    if segment==2:
        sim.b,s_minus,N ,delta_T,vneg_LSB, vj_LSB, i2_LSB, \
            vneg_MSB, vj_MSB, i2_MSB,= solving(sim,ring_mod,v,t,H)
    else:
        if SB=="LSB":
            sim.b,s_minus,N ,delta_T,vneg_LSB, vj_LSB, i2_LSB = solving(sim,ring_mod,v,t,H)
        if SB=="MSB":
            sim.b,s_minus,N ,delta_T,vneg_MSB, vj_MSB, i2_MSB = solving(sim,ring_mod,v,t,H)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    sim.save_data(ring_mod,t,v,file_name=filename)
    
    if segment==2:
        name = 'eye_SmallSignal_'+str(bit_num)+'bit_two_segment_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(-int(v.vpp*1000))+"mV_vbias_"+str(-int(v.v_bias*1000))+"mV"
        sim.eye_diagram(t,v,abs(s_minus)**2,filename=name+str(v.level),plot_bit_num=2,title="Output power (mW)")
        sim.eye_diagram(t,v,vneg_LSB,filename="vneg_MSB_SmallSignal"+str(v.level),plot_bit_num=2,title="vneg_LSB")
        sim.eye_diagram(t,v,vneg_MSB,filename="vneg_MSB_SmallSignal"+str(v.level),plot_bit_num=2,title="vneg_MSB")
        sim.eye_diagram(t,v,vj_LSB,filename="vj_MSB_SmallSignal"+str(v.level),plot_bit_num=2,title="vj_LSB")
        sim.eye_diagram(t,v,vj_MSB,filename="vj_MSB_SmallSignal"+str(v.level),plot_bit_num=2,title="vj_MSB")
    else:
        if SB=="LSB":
            name = 'eye_SmallSignal_'+str(bit_num)+'bit_LSB_f'+   \
            str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(-int(v.vpp*1000))+"mV_vbias_"+str(-int(v.v_bias*1000))+"mV"
            sim.eye_diagram(t,v,abs(s_minus)**2,filename=name+str(v.level),plot_bit_num=2,title="Output power (mW)")
            sim.eye_diagram(t,v,vneg_LSB,filename="vneg_SmallSignal"+str(v.level),plot_bit_num=2,title="vneg")
            sim.eye_diagram(t,v,vj_LSB,filename="vj_SmallSignal"+str(v.level),plot_bit_num=2,title="vj")
        if SB=="MSB":
            name = 'eye_SmallSignal_'+str(bit_num)+'bits_MSB_f'+   \
            str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(-int(v.vpp*1000))+"mV_vbias_"+str(-int(v.v_bias*1000))+"mV"
            sim.eye_diagram(t,v,abs(s_minus)**2,filename=name+str(v.level),plot_bit_num=2,title="Output power (mW)")
            sim.eye_diagram(t,v,vneg_MSB,filename="vneg_SmallSignal"+str(v.level),plot_bit_num=2,title="vneg")
            sim.eye_diagram(t,v,vj_MSB,filename="vj_SmallSignal"+str(v.level),plot_bit_num=2,title="vj")

    # sim.eye_diagram(t,v,v.v+vneg,filename="vpos + vneg SmallSignal"+str(v.level),plot_bit_num=2,title="vpos_plus_vneg")

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    