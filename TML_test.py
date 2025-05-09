import time as timer
from scipy.integrate import *
import numpy as np
from utility import *
from driver import *
from time_class import *
from ring  import *


v_bias = -1
vpp = 1
f_drive = 50
Cjs = [23.6*1e-15,20*1e-15]
Rs = 53.9
Rsi = 1439.8
Cox = 34.7e-15
Cp = 6.6e-15
f0 = 1e12
bit_num = 1000
level='NRZ'

# Some Dummy components
# //////////////////////////////////////////////////////////////////////////////////////////////////////////
experiment_condition ={"mode":"voltage_drive",
                        "lambda_incident":1.55,
                        "Pin":1} 
wl_in = 1.5467
Pin = 1 #mW
FSR = 0.0195
radius = 5
Amp_RoundTripLoss_pdk = [0.95124, 0.95248, 0.95273, 0.95292, 0.95305]
neff_pdk = [2.51105, 2.5111, 2.51113, 2.51116, 2.51118, 2.5112]
mode_area = 0.22*0.5

sim = simulation()
sim.main(experiment_condition=experiment_condition)

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
# ///////////////////////////////////////////////////////////////////////////////////////////////

v = driver(v_bias=v_bias,
           vpp=vpp,
           f_drive=f_drive,
           Rs = Rs,
           cj = Cjs,
           raise_cosine=1,
        #    sine_wave=1,
           PRBS=1,
           beta = 0.6,
           level=level)

t_class = time(mode = sim.mode)
t_class.main(ring_mod,v,N=bit_num,resolution=2)
sim.eye_diagram(t_class,v,title="vpos",signal=v.v,plot_bit_num=2,filename="sine_tml_test")

Z0 = 50
Cp_bar = Cp/t0
Cox_bar = Cox/t0
# Simple RC circuit test
# def TML(t,vj,driver):
#     voltage = driver.refering_v(t)
#     v_neg = ( 1/driver.Rs*vj - (1/driver.Rs-1/Z0)*voltage ) /(1/driver.Rs+1/Z0)
    
#     f1 = t0/driver.Cj_V(driver.v_bias)/Z0*(voltage - v_neg)

#     return f1

def TML_2(t,para,driver):
    v_neg , vj, i2 = para
    voltage = driver.refering_v(t)
    dvpos_dt = driver.refering_dv_dt(voltage,t_class,t)
    Rs = driver.Rs
    Cj_bar = driver.Cj_V(driver.v_bias)/t0
    dvneg_dt = ( 1/Cp_bar*1/Z0*(voltage - v_neg) \
                -dvpos_dt\
                -1/Cp_bar/Rs*(voltage + v_neg - vj)\
                -1/Cp_bar*i2)
    
    dvj_dt = 1/(Rs*Cj_bar)*(voltage + v_neg - vj)

    di2_dt = (1/Rsi*dvpos_dt + \
              1/Rsi*dvneg_dt\
               - 1/Rsi/Cox_bar*i2)
    
    return [dvneg_dt, dvj_dt, di2_dt]

# sol = solve_ivp(TML,t_span=[0,t_class.t_total[-1]], y0=[0+1j*0], t_eval=t_class.t_total,args= (v,),atol=1e-20,rtol=1e-15)
t1 = timer.time()
sol = solve_ivp(TML_2,t_span=[0,t_class.t_total[-1]], y0=[0+1j*0,0+1j*0,0+1j*0], t_eval=t_class.t_total,args= (v,),atol=1e-20,rtol=1e-15)
t2 = timer.time()
print("spent time = ",t2-t1," second")
v_neg = sol.y[0]
# vj = sol.y[0]
# v_neg = ( 1/v.Rs*vj - (1/v.Rs-1/Z0)*v.v ) /(1/v.Rs+1/Z0)
i1 = Cp*FDM(t_class.dt/t0,v.v+v_neg)
vj = sol.y[1]
i2 = sol.y[2]
i3 = 1/v.Rs*(v.v+v_neg-vj)
dis = 0
voltage = v.v[dis::]
v_neg = v_neg[dis::]
vj = vj[dis::]
i2 = i2[dis::]
i3 = i3[dis::]
t_array = t_class.t_total[dis::]
import os
os.chdir("TML")
ploting(t_array,i1,x_label="time (ps)",title="i1",filename='i1')
ploting(t_array,i2,x_label="time (ps)",title="i2",filename='i2')
ploting(t_array,i3,x_label="time (ps)",title="i3",filename='i3')
ploting(t_array,voltage+v_neg,x_label="time (ps)",title="v_on_modulator",filename='v_on_modulator')

sim.eye_diagram(t_class,v,title="Vj",signal=vj,plot_bit_num=2,filename="vj_test "+v.level)
sim.eye_diagram(t_class,v,title="Vneg",signal=v_neg,plot_bit_num=2,filename="vneg_test"+v.level)
sim.eye_diagram(t_class,v,title="V total",signal=v.v+v_neg,plot_bit_num=2,filename="v total_test"+v.level)