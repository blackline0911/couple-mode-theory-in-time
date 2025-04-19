from driver import driver
from time_class import *
from utility import *
from scipy.integrate import *
import os
import time as timer


# This simulation file only simulate driver performance

mode = "voltage_drive"
bit_num = 100
v_bias = -3
vpp = 5
wl_in = 1.5467
experiment_condition ={"mode":mode,
                        "lambda_incident":wl_in,
                        "Pin":1} 


sim = simulation()
sim.main(experiment_condition=experiment_condition)

# Simulating driver equation (Large Signal and Small Signal)
def Large_signal_Q(t_bar, Q_pround):
    voltage = np.real(v.refering_v(t_bar))
    # print(v.V_Q(Q_pround))
    f2 = (voltage/(v.R * v.cj_normalizing )*t0) \
        - (1/( v.R ) )*v.V_Q(Q_pround)*t0/v.cj_normalizing
    return f2 

def Small_signal_Q(t_bar, Q_pround):
    voltage = np.real(v.refering_v(t_bar))
    cj = v.Cj
    # print("Vj = ",Q_pround," V")
    f2 = (voltage/(v.R * v.cj_normalizing )*t0) \
        - (1/( v.R *cj) )*Q_pround*t0
    return f2 

# DeFine Solving process My Different simulation modes
def CMT_voltage_driving(driver,
                        time,func):
    Q_record = np.array([])
    Q_init=0
    # notes: time_range argument should be slightly exclude the t_eval
    sol =  solve_ivp(func,[0,time.t_all_segment[0][-1]],
                     [Q_init], 
                     t_eval=time.t_all_segment[0] ,
                     atol=atol,rtol=rtol)
    q = sol.y[0]
    Q_record = np.append(Q_record,sol.y[0])
    Q_init = sol.y[0][-1]
    
    for i in range(1,time.N):
        sol =  solve_ivp(func,[time.t_all_segment[i-1][-1] ,\
                        time.t_all_segment[i][-1]] ,
                         [Q_init], 
                         t_eval=np.append(  np.array([time.t_all_segment[i-1][-1]]), \
                                            np.array(time.t_all_segment[i])  ),
                         atol=atol,rtol=rtol)
        q = sol.y[0]
        Q_record = np.append(Q_record,q[1::])
        Q_init = sol.y[0][-1]
    Q_bar = Q_record

    return  Q_bar*driver.cj_normalizing   

v = driver(f_drive=50,
           v_bias=v_bias,
           vpp=vpp,
           R=53.9,
           raise_cosine=1,
           sine_wave=0,
           PRBS=1)

# Dummy ring component
ring_mod = ring(2*np.pi*5, 
            (0.95124),
            39.5,
            0.5*0.22,
            wl_in,
            (0.95105),
            FSR = 0.0195,
            neff=2.5111,
            FSR_shift=0,
            FCA_fit_factor=1,
            TPA_fit_factor=1,
            )

t = time(mode = sim.mode)
t.main(ring_mod,N=bit_num,driver=v)

sim.eye_diagram(t,v,v.v,
                filename="voltage_eye_"+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(vpp))+"_vbias_"+str(v_bias), 
                plot_bit_num=2)
Q_init = 0
print("\n\n.................Start Large Signal Simulation.................\n\n")
t1 = timer.time()
Q = CMT_voltage_driving(v,
                        t,Large_signal_Q)
t2 = timer.time()
print("Time for Large signal simulation is ",t2-t1," second")
print("\n\n.................Start Small Signal Simulation.................\n\n")
t1 = timer.time()
Q1 = CMT_voltage_driving(v,
                        t,Small_signal_Q)
t2 = timer.time()
print("Time for Small signal simulation is ",t2-t1," second")

os.chdir("./driver_test/")
filename = ('sim_'+str(int(v.f_drive/1e9))+"GHz_vpp_"+str(int(v.vpp))+"_vbias_"+str(v.v_bias))
sim.save_data(t,v,file_name=filename)

ploting(t.t_total,v.v,v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='voltage (V)',filename='Exteral_voltage_and_Large_signal_vj',leg=['External Voltage','Large Signal Vj'])
ploting(t.t_total,v.v,Q1/v.cj_normalizing,x_label='time (ps)',title='voltage (V)',filename='Exteral_voltage_and_small_signal_vj',leg=['External Voltage','Small Signal Vj'])
ploting(t.t_total,v.v-v.V_Q(Q/v.cj_normalizing),x_label='time (ps)',title='Resistor Voltage (V)',filename='Rv')
ploting(t.t_total,v.V_Q(Q/v.cj_normalizing),Q1/v.cj_normalizing,x_label='time (ps)',title='junction voltage',filename='junction voltage',leg=['large signal','small signal'])
ploting(t.t_total,Q,Q1,x_label='time (ps)',title='Q',filename='Q_vbias_'+str(int(v.v_bias))+"V",leg=['large signal','small signal'])
    
name = 'eye_LargeSignal_'+str(bit_num)+'bits_without_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp))+"_vbias_"+str(v.v_bias)
sim.eye_diagram(t,v,v.V_Q(Q/v.cj_normalizing),filename=name,plot_bit_num=2)

name = 'eye_SmallSignal_'+str(bit_num)+'bits_without_TPA_FCA_f'+   \
        str(int(v.f_drive*1e-9))+"GHz_vpp_"+str(int(v.vpp))+"_vbias_"+str(v.v_bias)
sim.eye_diagram(t,v,Q1/v.cj_normalizing,filename=name,plot_bit_num=2)
