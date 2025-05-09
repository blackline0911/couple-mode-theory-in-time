# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
from scipy.integrate import *
import numpy as np
from cmath import *
from utility import *
from driver import driver
from ring import ring
from time_class import time
from functools import partial

# relative solver tolerance 
rtol = 1e-14
# absolute solver tolerance
atol = 1e-20
# accuracy = atol + abs(y)*rtol

def CMT_large_signal(t_bar,eqs,SPM=None,FCA=None,vg_in_cm=None,ring=None,sim=None,driver=None,Heater=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    b_bar , Q_pround ,N_bar, v_neg, vj, i2 = eqs
    voltage = np.real(driver.refering_v(t_bar))
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    alpha_linear = ring.alpha(voltage)
    cj = driver.Cj_V(voltage)
    cj_bar = cj/t0
    dlambda = ring.lambda0/ring.ng*( ring.neff(vj) - ring.neff(0))
    Rs = driver.Rs
    Cox_bar = driver.Cox/t0
    Rsi = driver.Rsi
    Cp_bar = driver.Cp/t0

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # Coupled Transmission Line Equations
    dvneg_dt, dvj_dt, di2_dt = driver.TML(v_neg, vj, voltage, dvpos_dt, i2, \
                                          driver.Z0, Rs, cj_bar, Cox_bar, Rsi, Cp_bar)
    
    f1 = 1j*2*np.pi*(ring.f_res_bar-sim.f_pround_bar + \
                     SPM*abs(b_bar)**2)*b_bar- \
    \
        (ring.tu_e_bar_total_inv + \
            ring.vg_in_cm*alpha_linear*(1 + ring.TPA_ratio/2*(sim.b0)**2*abs(b_bar)**2\
    \
            + N_bar*1e-5/2/(ring.vg_in_cm*alpha_linear) ) )*b_bar + \
                ring.input_kappa *1 + \
    \
        1j*ring.D_bar*(dlambda)*b_bar + \
    \
        1j*ring.D_bar*( ring.dlambda_dT*(Heater.T_surround-300) + (ring.HE*1e-6)*Heater.P)
        
    f2 = (voltage/(driver.Rs * driver.cj_normalizing )*t0) \
        - (1/( driver.Rs ) )*driver.V_Q(Q_pround)*t0/driver.cj_normalizing

    f3 = -N_bar/ring.tau_eff*(t0/1e-9) + FCA*abs(b_bar)**4

    return [ f1,f2 ,f3, dvneg_dt, dvj_dt, di2_dt]

def CMT_small_signal(t_bar,eqs,SPM=None,FCA=None,ring=None,sim=None,driver=None,Heater=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    b_bar , Q_pround ,N_bar, v_neg, vj, i2= eqs
    voltage = np.real(driver.refering_v(t_bar))
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    cj = driver.Cj
    cj_bar = cj/t0
    alpha_linear = ring.alpha(voltage)
    dlambda = ring.lambda0/ring.ng*( ring.neff(vj) - ring.neff(0))
    Rs = driver.Rs
    Cox_bar = driver.Cox/t0
    Rsi = driver.Rsi
    Cp_bar = driver.Cp/t0

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # Coupled Transmission Line Equations
    dvneg_dt, dvj_dt, di2_dt = driver.TML(v_neg, vj, voltage, dvpos_dt, i2, \
                                          driver.Z0,Rs, cj_bar,Cox_bar,Rsi,Cp_bar)
    
    # Couple mode Equations in time domain
    f1 = 1j*2*np.pi*(ring.f_res_bar - sim.f_pround_bar + \
                      SPM*abs(b_bar)**2)*b_bar- \
    \
        (ring.tu_e_bar_total_inv + \
            ring.vg_in_cm*alpha_linear*(1 + ring.TPA_ratio/2*(sim.b0)**2*abs(b_bar)**2\
    \
            + N_bar*1e-5/2/(ring.vg_in_cm*alpha_linear) ) )*b_bar + \
                ring.input_kappa *1 + \
    \
        1j*ring.D_bar*(dlambda)*b_bar + \
    \
        1j*ring.D_bar*( ring.dlambda_dT*(Heater.T_surround-300) + (ring.HE*1e-6)*Heater.P)

    f2 = (voltage/(driver.Rs * driver.cj_normalizing )*t0) \
        - (1/( driver.Rs*cj)*t0 )*Q_pround

    f3 = -N_bar/ring.tau_eff*(t0/1e-9) + FCA*abs(b_bar)**4

    return [ f1,f2 ,f3, dvneg_dt, dvj_dt, di2_dt]

def CMT_scan_frequency(t_bar,eqs,SPM,FCA,ring,sim,driver,Heater=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    f_res_bar = ring.w_res(t_bar)
    b_bar , N_bar = eqs
    voltage = driver.v_bias
    alpha_linear = ring.alpha(voltage)
    dlambda = ring.lambda0/ring.ng*( ring.neff(voltage) - ring.neff(0))
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # 不要忘記這裡的tau_o吸收是amplitude的吸收，不是energy的
    f1 = 1j*2*np.pi*(f_res_bar-sim.f_pround_bar + SPM*abs(b_bar)**2)*b_bar \
        - (ring.tu_e_bar_total_inv + \
        \
        ring.vg_in_cm*alpha_linear*(1 + ring.TPA_ratio/2*(sim.b0)**2*abs(b_bar)**2\
        + N_bar*1e-5/2/(ring.vg_in_cm*alpha_linear) ) )*b_bar + \
        \
        ring.input_kappa *1 + \
        \
        1j*ring.D_bar*dlambda*b_bar +\
        \
        1j*ring.D_bar*( ring.dlambda_dT*(Heater.T_surround-300) + (ring.HE*1e-6)*Heater.P)*b_bar

    # 這裡不加Q的方程式並不會影響結果

    f3 = -N_bar/ring.tau_eff*(t0/1e-9) + FCA*abs(b_bar)**4
    return [ f1,f3]

# DeFine Solving process My Different simulation modes
def CMT_voltage_driving(sim,ring, 
                        driver,
                        time,
                        Heater):
    method_dict = {
        'large_signal' : CMT_large_signal,
        'small_signal' : CMT_small_signal,
    }
    b_record = np.array([])
    Q_record = np.array([])
    N_record = np.array([])
    vneg_record = np.array([])
    vj_record = np.array([])
    i2_record = np.array([])
    b_init=0+1j*0
    Q_init=0
    N_init=0
    vneg_init=0
    vj_init=0
    i2_init=0
    # notes: time_range argument should be slightly exclude the t_eval
    FCA = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*(sim.b0)**4
    SPM = ring.dw_SPM_coeff * (sim.b0)**2
    ode_func = partial(method_dict[driver.method],
                       SPM=SPM,FCA=FCA,
                       ring = ring,sim=sim,
                       driver=driver,Heater=Heater)
    sol =  solve_ivp(ode_func ,[0 ,time.t_all_segment[0][-1]], 
                            [b_init, Q_init,N_init, vneg_init, vj_init, i2_init],
                            method=sim.algorithm,
                            t_eval = time.t_all_segment[0] ,
                            atol = atol,rtol = rtol,)
    b = sol.y[0]
    b_record = np.append(b_record,sol.y[0])
    q = sol.y[1]
    Q_record = np.append(Q_record,sol.y[1])
    n = sol.y[2]
    N_record = np.append(N_record,sol.y[2])
    vneg = sol.y[3]
    vneg_record = np.append(vneg_record,vneg)
    vj = sol.y[4]
    vj_record = np.append(vj_record,vj)
    i2 = sol.y[5]
    i2_record = np.append(i2_record,i2)

    b_init = sol.y[0][-1]
    Q_init = sol.y[1][-1]
    N_init = sol.y[2][-1]
    vneg_init = sol.y[3][-1]
    vj_init = sol.y[4][-1]
    i2_init = sol.y[5][-1]
    
    for i in range(1,time.N):
        sol =  solve_ivp(ode_func ,
                        [time.t_all_segment[i-1][-1] ,\
                        time.t_all_segment[i][-1]], 
                        [b_init, Q_init,N_init, vneg_init, vj_init, i2_init],
                        method=sim.algorithm,
                        t_eval = np.append(  np.array([time.t_all_segment[i-1][-1]]), \
                                            np.array(time.t_all_segment[i])  ),
                        atol = atol,
                        rtol = rtol,)
        b = sol.y[0]
        b_record = np.append(b_record,b[1::])
        q = sol.y[1]
        Q_record = np.append(Q_record,q[1::])
        n = sol.y[2]
        N_record = np.append(N_record,n[1::])
        vneg = sol.y[3]
        vneg_record = np.append(vneg_record,vneg[1::])
        vj = sol.y[4]
        vj_record = np.append(vj_record,vj[1::])
        i2 = sol.y[5]
        i2_record = np.append(i2_record,i2[1::])
        b_init = sol.y[0][-1]
        Q_init = sol.y[1][-1]
        N_init = sol.y[2][-1]
        vneg_init = sol.y[3][-1]
        vj_init = sol.y[4][-1]
        i2_init = sol.y[5][-1]
    b_bar = b_record
    Q_bar = Q_record
    N_bar = N_record
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*sim.b0, Q_bar*driver.cj_normalizing , s_minus_bar*sim.S0    ,N_bar/(ring.sigma_FCA*1e-17)*1e-5  ,\
            vneg_record, vj_record, i2_record

def solve_scan_frequency(sim,ring, 
                        driver,
                        time,
                        Heater):
    SPM = ring.dw_SPM_coeff *(sim.b0)**2
    FCA = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*(sim.b0)**4
    b_init = 0+1j*0
    N_init = 0
    sol = solve_ivp(CMT_scan_frequency ,
                    [0,time.t_max+time.buffer], 
                    [b_init,N_init],
                    method=sim.algorithm,
                    t_eval = time.t_total,atol = atol,rtol = rtol,
                    args=(SPM,FCA,ring,sim,driver,Heater))
    b_bar = sol.y[0]
    N_bar = sol.y[1]
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*sim.b0, s_minus_bar*sim.S0    ,N_bar/(ring.sigma_FCA*1e-17)*1e-5

def solving(sim,
            ring, 
            driver,
            time,
            Heater
            ):
    """
    ring : ring object
    driver: driver object
    time: time object
    lambda_incident: incident laser wavelength (micron)
    Pin : Power of input laser (mWatt)
    """
#     """Since the length of solution function array may not be the same as t_eval argument we specified when the time in single solve_ivp is long.
#     Hence, we divide the time according to the Baud Rate, and solve coupled differential equation by each time segments. 
#     """
    mode_dict = {
        "voltage_drive" : CMT_voltage_driving,
        "scan_frequency" : solve_scan_frequency
    }

    solver = partial(mode_dict[sim.mode],sim,ring,driver,time,Heater)
    
    if sim.mode == "scan_frequency":
        b,  s_minus, N= solver()
        return b,  s_minus   ,N
    else:
        b, Q, s_minus, N, vneg, vj, i2 = solver()
        return b, Q, s_minus   ,N, vneg, vj, i2
