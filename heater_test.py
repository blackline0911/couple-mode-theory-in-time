from ring import * 
from driver import *
from time_class import *
from cmt_solver import *
from utility import *
from Heater import Heater


wl_in = 1.5595
Pin = 1
mode = "scan_frequency"
radius = 5
FSR = 0.019431
ng = None
wl_res = None
neff = 2.680503
me = 39.3
mode_cross_section = 0.22*0.5

a = 0.953553     # Round trip amplitude loss
gamma = 0.947582  # couple through amplitude coefficient

# ///////////////////////////////////////
# ///////////////////////////////////////
# Equivalent circuit
Rs=68.4
cj = [23.8e-15, 20.1e-15]

# ///////////////////////////////////////
# ///////////////////////////////////////
# Driver
v_bias = 0
vpp = 0
f_drive_GHz = 50
bit_num=30

# ///////////////////////////////////////
# ///////////////////////////////////////
# Heater
T_surround = 300
V_heater = 0
R_heater = 0.5*150/0.6

# ///////////////////////////////////////
# ///////////////////////////////////////
# Simulation Parameter
wl_min =  1.5591
wl_max =  1.565
result_folder = 'heater_test'
time_resolution = 1
t_max = 5000


# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

exp_condition = {
    'mode':mode,
    'Pin' :Pin,
    'lambda_incident':wl_in
}
sim = simulation()
sim.main(experiment_condition=exp_condition)

ring_mod = ring(2*np.pi*radius,
                a,
                me,
                mode_cross_section,
                wl_in,
                gamma,
                FSR=FSR,
                lambda0=wl_res,
                neff = neff,
                ng = ng
                )

v = driver(f_drive=f_drive_GHz,
           v_bias=v_bias,
           vpp=vpp,
           R=Rs,
           raise_cosine=1,
           sine_wave=0,
           cj = np.array(cj),
           PRBS=1)

H = Heater(T_surround=T_surround,
           V = V_heater,
           R=R_heater)


if not os.path.isdir(result_folder):
    os.mkdir(result_folder)
os.chdir(result_folder)

t = time(mode = sim.mode)
if sim.mode == "scan_frequency":
    
    V_Heater = np.array([0,-0.5,1])
    
    plt.figure()
    i = 0
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=t_max,resolution=time_resolution,buffer=500,driver=v)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    sim.save_data(ring_mod,t,v)
    for V in V_Heater:
        H.V = V
        H.renew()
        print(H.V)
        print(H.R)
        print("Power = ",H.P," mW")
        # ring_mod.lambda_incident = wl_in+ring_mod.HE*1e-6*H.P
        
        print(ring_mod.lambda_incident)
        b,Q,s_minus,N = solving(sim,ring_mod,v,t,H)
        # ploting(t.t_total,abs(s_minus)**2,x_label="time (ps)",title="|S-|^2")
        
        T = Transfer_function(ring_mod,t)
        wl,data = T.mapping((abs(s_minus)**2/sim.Pin))
        wl,data_phase = T.mapping(180/np.pi*np.angle(s_minus))
        plt.plot(wl*1000,dB(data),
                    label=str((V)))
    i+=1
    
    plt.grid(color='g',linestyle='--', alpha=0.5)
    plt.xlabel('wavelength(nm)')
    plt.title('Transfer function')
    plt.legend()
    plt.savefig("Tranfer_function")
    plt.show()
if sim.mode == "voltage_drive":
    t.main(ring_mod,N=bit_num,driver=v)
    print("\n\nSimulation at ",str(v.f_drive/1e9)," GHz, ",str(v.vpp),"V vpp, ",str(v.v_bias),"V vbias\n\n")
    filename = ('sim_'+str(int(f_drive_GHz))+"GHz_vpp_"+str(int(vpp))+"_"+str(int(10*vpp))+"_vbias_"+str(v_bias))
    sim.save_data(ring_mod,t,v,file_name=filename)
    v.method = "large_signal"
    sim.eye_diagram(t,v,v.v,
                    filename="voltage_eye_"+str(int(f_drive_GHz))+"GHz_vpp_"+str(int(vpp))+"_"+str(int(10*vpp))+"_vbias_"+str(v_bias), 
                    plot_bit_num=2)
    
    print("\n\n.................Start Large Signal Simulation.................\n\n")
    b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    

    name = 'eye_LargeSignal_'+str(bit_num)+'bits_without_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp))+"_"+str(int(10*vpp))+"_vbias_"+str(v.v_bias)
    sim.eye_diagram(t,v,abs(s_minus)**2,filename=name,plot_bit_num=2)

    print("\n\n.................Start Small Signal Simulation.................\n\n")
    v.method = "small_signal"
    b1,Q1,s_minus1,N1 = solving(sim,ring_mod,v,t)
    ploting(t.t_total,v.v,v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='voltage (V)',filename='voltage_and_Large_signal_vj',leg=['External Voltage','Large Signal Vj'])
    ploting(t.t_total,v.v,Q1/v.cj_normalizing,x_label='time (ps)',title='voltage (V)',filename='voltage_and_small_signal_vj',leg=['External Voltage','Small Signal Vj'])
    ploting(t.t_total,v.v-v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='Resistor Voltage (V)',filename='Rv')
    ploting(t.t_total,Q,Q1,x_label='time (ps)',title='Q',filename='Q',leg=['Large Signal','Small Signal'])
    ploting(t.t_total,v.V_Q(Q/v.cj_normalizing),Q1/v.cj_normalizing,x_label='time (ps)',title='junction voltage',filename='junction voltage',leg=['large signal','small signal'])
    ploting(t.t_total,abs(s_minus)**2,abs(s_minus1)**2,x_label='time (ps)',title='output Power',filename='output Power',leg=['large signal','small signal'])

    name = 'eye_SmallSignal_'+str(bit_num)+'bits_without_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp))+"_vbias_"+str(v.v_bias)
    sim.eye_diagram(t,v,abs(s_minus1)**2,filename=name,plot_bit_num=2)
