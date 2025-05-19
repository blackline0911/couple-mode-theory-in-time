# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
from scipy.integrate import *
import numpy as np
from cmath import *
from utility import *
from driver import driver
from ring import ring
from time_class import time
from functools import partial
from sim import simulation
from Heater import Heater

# relative solver tolerance 
rtol = 1e-14
# absolute solver tolerance
atol = 1e-20
# accuracy = atol + abs(y)*rtol

def CMT_large_signal(t_bar,eqs,driver:driver,ring:ring,Heater:Heater,SPM=None,TPA=None,FCA=None,sim=None,T_args=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    b_bar , Q_pround ,N_bar,delta_T, v_neg, vj, i2= eqs
    voltage = np.real(driver.refering_v(t_bar))
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    alpha_linear = ring.alpha(vj)
    cj = driver.Cj_V(vj)
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
    
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,delta_T,ring.f_res_bar,alpha_linear,TPA,SPM,T_args,dlambda,Heater)
    
    dQ_dt = (voltage/(driver.Rs * driver.cj_normalizing )*t0) \
        - (1/( driver.Rs ) )*driver.V_Q(Q_pround)*t0/driver.cj_normalizing
    
    dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)

    dT_dt = Heater.T_rate_equation(b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim)
    return [ da_dt,dQ_dt ,dN_dt, dT_dt, dvneg_dt, dvj_dt, di2_dt]

def CMT_small_signal(t_bar,eqs,driver:driver,SPM=None,TPA=None,FCA=None,ring=None,sim=None,Heater=None,T_args=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    b_bar , Q_pround ,N_bar,delta_T, v_neg, vj, i2= eqs
    voltage = np.real(driver.refering_v(t_bar))
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    cj = driver.Cj
    cj_bar = cj/t0
    alpha_linear = ring.alpha(vj)
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
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,delta_T,ring.f_res_bar,alpha_linear,TPA,SPM,T_args,dlambda,Heater)

    dQ_dt = (voltage/(driver.Rs * driver.cj_normalizing )*t0) \
        - (1/( driver.Rs*cj)*t0 )*Q_pround

    dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)

    dT_dt = Heater.T_rate_equation(b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim)

    return [ da_dt,dQ_dt  ,dN_dt ,dT_dt , dvneg_dt, dvj_dt, di2_dt]

def CMT_scan_frequency(t_bar,eqs,SPM,TPA,FCA,T_args,ring:ring,sim:simulation,driver:driver,Heater:Heater):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    f_res_bar = ring.w_res(t_bar)
    b_bar  = eqs
    voltage = driver.v_bias
    alpha_linear = ring.alpha(voltage)
    df = -sim.f_pround_bar/ring.ng*( ring.neff(voltage) - ring.neff(0))*(ring.L_active/ring.L)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    #         + tu_o_bar*FCA/2 * abs(b_bar)**4 ) )*b_bar + \
    # 不要忘記這裡的tau_o吸收是energy的吸收，不是amplitude的

    N_bar = FCA* abs(b_bar)**4

    delta_T = (sim.b0*abs(b_bar))**2 / (T_args[1]*T_args[2]*T_args[3]) * (\
        \
        +ring.vg_in_cm*alpha_linear \
        \
        + TPA*abs(b_bar)**2 \
        \
        + ring.vg_in_cm * N_bar*1e-5 \
        ) * Heater.tau_th
    
    # print("FCA absorption = ",N_bar*1e-5*abs(b_bar)**4,' 1/cm')
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,delta_T,f_res_bar,alpha_linear,TPA,SPM,T_args,\
                     (-c*1e-12/sim.f_pround_bar**2)*df , Heater)
    # print((-c*1e-12/sim.f_pround_bar**2)*df)

    # N_bar = sigma_FCA*N*1e5

    # dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)
    # print("N =",N_bar*1/(ring.sigma_FCA*1e-17)*1e-5)
    # print("FCA absorb =",N_bar*1e-5," 1/cm")
    # print("FCA absorb 1 =",ring.FCA_coeff*abs(sim.b0*b_bar)**4," 1/cm")
    return [ da_dt ]

# DeFine Solving process My Different simulation modes
def CMT_voltage_driving(sim,ring:ring, 
                        driver,
                        time,
                        Heater):
    method_dict = {
        'large_signal' : CMT_large_signal,
        'small_signal' : CMT_small_signal,
    }
    Aeff = 0.204 #um^2
    Veff = Aeff*ring.L #um^3
    ro_si = 2.329e-12 # g/um^3
    cSi = 0.713*1000 # mJ/(g*K)
    kappa_theta = 1.86e-4 # 1/K
    T_args = [kappa_theta, ro_si, cSi, Veff]
    b_record = np.array([])
    Q_record = np.array([])
    N_record = np.array([])
    T_record = np.array([])
    vneg_record = np.array([])
    vj_record = np.array([])
    i2_record = np.array([])
    b_init=0+1j*0
    Q_init=0
    N_init=0
    delta_T_init = 0
    vneg_init=0
    vj_init=0
    i2_init=0
    # notes: time_range argument should be slightly exclude the t_eval
    TPA = ring.TPA_coeff*(sim.b0)**2
    FCA = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*(sim.b0)**4
    SPM = ring.df_SPM_coeff *(sim.b0)**2
    ode_func = partial(method_dict[driver.method],
                       SPM=SPM,FCA=FCA,TPA=TPA,
                       T_args=T_args,
                       ring = ring,sim=sim,
                       driver=driver,Heater=Heater)
    sol =  solve_ivp(ode_func ,[0 ,time.t_all_segment[0][-1]], 
                            [b_init, Q_init,N_init, delta_T_init, vneg_init, vj_init, i2_init],
                            method=sim.algorithm,
                            t_eval = time.t_all_segment[0] ,
                            atol = atol,rtol = rtol,)
    b = sol.y[0]
    b_record = np.append(b_record,sol.y[0])
    q = sol.y[1]
    Q_record = np.append(Q_record,sol.y[1])
    n = sol.y[2]
    N_record = np.append(N_record,sol.y[2])
    T = sol.y[3]
    T_record = np.append(T_record,sol.y[3])
    vneg = sol.y[4]
    vneg_record = np.append(vneg_record,vneg)
    vj = sol.y[5]
    vj_record = np.append(vj_record,vj)
    i2 = sol.y[6]
    i2_record = np.append(i2_record,i2)

    b_init = sol.y[0][-1]
    Q_init = sol.y[1][-1]
    N_init = sol.y[2][-1]
    delta_T_init = sol.y[3][-1]
    vneg_init = sol.y[4][-1]
    vj_init = sol.y[5][-1]
    i2_init = sol.y[6][-1]
    
    for i in range(1,time.N):
        sol =  solve_ivp(ode_func ,
                        [time.t_all_segment[i-1][-1] ,\
                        time.t_all_segment[i][-1]], 
                        [b_init, Q_init,N_init, delta_T_init, vneg_init, vj_init, i2_init],
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
        T = sol.y[3]
        T_record = np.append(T_record,T[1::])
        vneg = sol.y[4]
        vneg_record = np.append(vneg_record,vneg[1::])
        vj = sol.y[5]
        vj_record = np.append(vj_record,vj[1::])
        i2 = sol.y[6]
        i2_record = np.append(i2_record,i2[1::])
        b_init = sol.y[0][-1]
        Q_init = sol.y[1][-1]
        N_init = sol.y[2][-1]
        delta_T_init = sol.y[3][-1]
        vneg_init = sol.y[4][-1]
        vj_init = sol.y[5][-1]
        i2_init = sol.y[6][-1]
    b_bar = b_record
    Q_bar = Q_record
    N_bar = N_record
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*sim.b0, Q_bar*driver.cj_normalizing , s_minus_bar*sim.S0    ,N_bar/(ring.sigma_FCA*1e-17)*1e-5  , T_record,\
            vneg_record, vj_record, i2_record

def solve_scan_frequency(sim,ring:ring, 
                        driver:driver,
                        time:time,
                        Heater:Heater):
    
    # 非線性吸收項記得乘以vg
    SPM = ring.df_SPM_coeff *(sim.b0)**2
    # FCA項因為還要算dN_dt，所以先不用乘以vg
    FCA = ring.FCA_coeff*1e5*(sim.b0)**4
    # FCA = ring.FCA_coeff*t0*(sim.b0)**4
    TPA = ring.TPA_coeff*(sim.b0)**2 
    Aeff = 0.204 #um^2
    Veff = Aeff*ring.L #um^3
    ro_si = 2.329e-12 # g/um^3
    cSi = 0.713*1000 # mJ/(g*K)
    kappa_theta = 1.86e-4 # 1/K
    T_args = [kappa_theta, ro_si, cSi, Veff]
    b_init = 0+1j*0
    N_init = 0
    delta_T_init = 0
    sol = solve_ivp(CMT_scan_frequency ,
                    [0,time.t_max+time.buffer], 
                    [b_init],
                    method=sim.algorithm,
                    t_eval = time.t_total,atol = atol,rtol = rtol,
                    args=(SPM,TPA,FCA,T_args,ring,sim,driver,Heater))
    b_bar = sol.y[0]
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*sim.b0, s_minus_bar*sim.S0    
    # return b_bar*sim.b0, s_minus_bar*sim.S0

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
        b,  s_minus  = solver()
        return b,  s_minus
    else:
        b, Q, s_minus, N, T,vneg, vj, i2 = solver()
        return b, Q, s_minus   ,N, T,vneg, vj, i2
