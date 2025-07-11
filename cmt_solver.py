# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
from scipy.integrate import *
from numpy import append, array
from cmath import *
from utility import t0
from driver import driver
from ring import ring
from time_class import time
from functools import partial
from Heater import Heater
from CMTs import CMT_large_signal_single_segment, CMT_large_signal_dual_segment, \
    CMT_small_signal_single_segment, CMT_small_signal_dual_segment, CMT_scan_frequency

# relative solver tolerance 
rtol = 1e-14
# absolute solver tolerance
atol = 1e-20
# accuracy = atol + abs(y)*rtol


# DeFine Solving process by Different simulation modes
def CMT_voltage_driving(sim,ring:ring, 
                        driver,
                        time,
                        Heater):
    method_dict = {
        'large_signal_single_segment' : CMT_large_signal_single_segment,
        'large_signal_dual_segment' : CMT_large_signal_dual_segment,
        'small_signal_single_segment' : CMT_small_signal_single_segment,
        'small_signal_dual_segment' : CMT_small_signal_dual_segment,
    }
    Veff = Heater.Aeff*ring.L #um^3
    T_coeff = sim.b0**2/( ring.ro_si*ring.cSi*Veff)
    T_args = [ring.kappa_thermal,T_coeff]
    
    b_init=0+1j*0
    N_init=0
    delta_T_init = 0
    # notes: time_range argument should be slightly exclude the t_eval
    TPA = ring.TPA_coeff*(sim.b0)**2
    FCA = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*(sim.b0)**4
    SPM = ring.df_SPM_coeff *(sim.b0)**2
    if ring.segment == 1:
        b_record    =  array([])
        N_record    =  array([])
        T_record    =  array([])
        vneg_record =  array([])
        vj_record   =  array([])
        i2_record   =  array([])
        vneg_init   =   0
        vj_init     =   0
        i2_init     =   0
        func_name = driver.method+ '_single_segment'
        init = [b_init, N_init, delta_T_init, vneg_init, vj_init, i2_init]
        ode_func = partial(method_dict[func_name],
                        SPM=SPM,FCA=FCA,TPA=TPA,
                        T_args=T_args,
                        ring = ring,sim=sim,
                        driver=driver,Heater=Heater)
        sol =  solve_ivp(ode_func ,[0 ,time.t_all_segment[0][-1]], 
                                init,
                                method=sim.algorithm,
                                t_eval = time.t_all_segment[0] ,
                                atol = atol,rtol = rtol,)
        b    =   sol.y[0]
        n    =   sol.y[1]
        T    =   sol.y[2]
        vneg =   sol.y[3]
        vj   =   sol.y[4]
        i2   =   sol.y[5]
        b_record    =   append(b_record,b)
        N_record    =   append(N_record,n)
        T_record    =   append(T_record,T)
        vneg_record =   append(vneg_record,vneg)
        vj_record   =   append(vj_record,vj)
        i2_record   =   append(i2_record,i2)

        b_init       =  sol.y[0][-1]
        N_init       =  sol.y[1][-1]
        delta_T_init =  sol.y[2][-1]
        vneg_init    =  sol.y[3][-1]
        vj_init      =  sol.y[4][-1]
        i2_init      =  sol.y[5][-1]
        
        for i in range(1,time.N):
            sol =  solve_ivp(ode_func ,
                            [time.t_all_segment[i-1][-1] ,\
                            time.t_all_segment[i][-1]], 
                            [b_init, N_init, delta_T_init, vneg_init, vj_init, i2_init],
                            method=sim.algorithm,
                            t_eval = append(array([time.t_all_segment[i-1][-1]]), \
                                            array(time.t_all_segment[i])),
                            atol   = atol,
                            rtol   = rtol,)
            b =     sol.y[0]
            n =     sol.y[1]
            T =     sol.y[2]
            vneg =  sol.y[3]
            vj =    sol.y[4]
            i2 =    sol.y[5]
            b_record     =   append(b_record,b[1::])
            N_record     =   append(N_record,n[1::])
            T_record     =   append(T_record,T[1::])
            vneg_record  =   append(vneg_record,vneg[1::])
            vj_record    =   append(vj_record,vj[1::])
            i2_record    =   append(i2_record,i2[1::])
            b_init       =   sol.y[0][-1]
            N_init       =   sol.y[1][-1]
            delta_T_init =   sol.y[2][-1]
            vneg_init    =   sol.y[3][-1]
            vj_init      =   sol.y[4][-1]
            i2_init      =   sol.y[5][-1]

        s_minus_bar = (1-ring.input_kappa*b_record)

        return b_record*sim.b0, s_minus_bar*sim.S0,             \
               N_record/(ring.sigma_FCA*1e-17)*1e-5, T_record,  \
               vneg_record, vj_record, i2_record
    if ring.segment == 2:
        vneg_LSB_init   =   0
        vj_LSB_init     =   0
        i2_LSB_init     =   0
        vneg_MSB_init   =   0
        vj_MSB_init     =   0
        i2_MSB_init     =   0
        b_record =          array([])
        N_record =          array([])
        T_record =          array([])
        vneg_LSB_record =   array([])
        vj_LSB_record   =   array([])
        i2_LSB_record   =   array([])
        vneg_MSB_record =   array([])
        vj_MSB_record   =   array([])
        i2_MSB_record   =   array([])
        func_name = driver.method+ '_dual_segment'
        init = [b_init, N_init, delta_T_init, vneg_LSB_init, vj_LSB_init, i2_LSB_init,
                vneg_MSB_init, vj_MSB_init, i2_MSB_init]
        ode_func = partial(method_dict[func_name],
                        SPM=SPM,FCA=FCA,TPA=TPA,
                        T_args=T_args,
                        ring = ring,sim=sim,
                        driver=driver,Heater=Heater,)
        sol =  solve_ivp(ode_func ,[0 ,time.t_all_segment[0][-1]], 
                                init,
                                method=sim.algorithm,
                                t_eval = time.t_all_segment[0] ,
                                atol = atol,rtol = rtol,)
        b               =  sol.y[0]
        n               =  sol.y[1]
        T               =  sol.y[2]
        vneg_LSB        =  sol.y[3]
        vj_LSB          =  sol.y[4]
        i2_LSB          =  sol.y[5]
        vneg_MSB        =  sol.y[6]
        vj_MSB          =  sol.y[7]
        i2_MSB          =  sol.y[8]
        b_record    =   append(b_record,b)
        N_record    =   append(N_record,n)
        T_record    =   append(T_record,T)
        vneg_LSB_record =   append(vneg_LSB_record,vneg_LSB)
        vj_LSB_record   =   append(vj_LSB_record,vj_LSB)
        i2_LSB_record   =   append(i2_LSB_record,i2_LSB)
        vneg_MSB_record =   append(vneg_MSB_record,vneg_MSB)
        vj_MSB_record   =   append(vj_MSB_record,vj_MSB)
        i2_MSB_record   =   append(i2_MSB_record,i2_MSB)

        b_init       =  sol.y[0][-1]
        N_init       =  sol.y[1][-1]
        delta_T_init =  sol.y[2][-1]
        vneg_LSB_init = sol.y[3][-1]
        vj_LSB_init =   sol.y[4][-1]
        i2_LSB_init =   sol.y[5][-1]
        vneg_MSB_init = sol.y[6][-1]
        vj_MSB_init =   sol.y[7][-1]
        i2_MSB_init =   sol.y[8][-1]
        for i in range(1,time.N-driver.shift_bit):
            sol =  solve_ivp(ode_func ,
                            [time.t_all_segment[i-1][-1] ,\
                            time.t_all_segment[i][-1]], 
                            [b_init, N_init, delta_T_init,\
                              vneg_LSB_init, vj_LSB_init, i2_LSB_init,\
                              vneg_MSB_init, vj_MSB_init, i2_MSB_init,],
                            method=sim.algorithm,
                            t_eval =  append(array([time.t_all_segment[i-1][-1]]), \
                                            array(time.t_all_segment[i])),
                            atol = atol,
                            rtol = rtol,)
            b               =  sol.y[0]
            n               =  sol.y[1]
            T               =  sol.y[2]
            vneg_LSB        =  sol.y[3]
            vj_LSB          =  sol.y[4]
            i2_LSB          =  sol.y[5]
            vneg_MSB        =  sol.y[6]
            vj_MSB          =  sol.y[7]
            i2_MSB          =  sol.y[8]
            b_record        =   append(b_record,b[1::])
            N_record        =   append(N_record,n[1::])
            T_record        =   append(T_record,T[1::])
            vneg_LSB_record =   append(vneg_LSB_record,vneg_LSB[1::])
            vj_LSB_record   =   append(vj_LSB_record,vj_LSB[1::])
            i2_LSB_record   =   append(i2_LSB_record,i2_LSB[1::])
            vneg_MSB_record =   append(vneg_MSB_record,vneg_MSB[1::])
            vj_MSB_record   =   append(vj_MSB_record,vj_MSB[1::])
            i2_MSB_record   =   append(i2_MSB_record,i2_MSB[1::])
            b_init =        sol.y[0][-1]
            N_init =        sol.y[1][-1]
            delta_T_init =  sol.y[2][-1]
            vneg_LSB_init = sol.y[3][-1]
            vj_LSB_init =   sol.y[4][-1]
            i2_LSB_init =   sol.y[5][-1]
            vneg_MSB_init = sol.y[6][-1]
            vj_MSB_init =   sol.y[7][-1]
            i2_MSB_init =   sol.y[8][-1]
        s_minus_bar = (1-ring.input_kappa*b_record)
        return  b_record*sim.b0, s_minus_bar*sim.S0,               \
                N_record/(ring.sigma_FCA*1e-17)*1e-5,  T_record,   \
                vneg_LSB_record, vj_LSB_record, i2_LSB_record,     \
                vneg_MSB_record, vj_MSB_record, i2_MSB_record

def solve_scan_frequency(sim,ring:ring, 
                        driver:driver,
                        time:time,
                        Heater:Heater):
    # method_dict = {
    #     'single_segment' : CMT_scan_frequency_single_segment,
    #     'dual_segment' : CMT_scan_frequency,
    # }
    SPM = ring.df_SPM_coeff *(sim.b0)**2
    # 計算FCA穩態解
    FCA     = ring.FCA_coeff*1e5*(sim.b0)**4
    TPA     = ring.TPA_coeff*(sim.b0)**2 
    Veff    = Heater.Aeff*ring.L
    T_coeff = sim.b0**2/( ring.ro_si*ring.cSi*Veff)
    T_args  = [ring.kappa_thermal,T_coeff]
    b_init  = 0+1j*0
    
    sol = solve_ivp( CMT_scan_frequency ,
                    [0,time.t_max+time.buffer], 
                    [b_init],
                    method=sim.algorithm,
                    t_eval = time.t_total,atol = atol,rtol = rtol,
                    args=(SPM,TPA,FCA,T_args,ring,sim,driver,Heater))
    b_bar = sol.y[0]
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*sim.b0, s_minus_bar*sim.S0

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
        if ring.segment == 1:
            b, s_minus, N, T, vneg, vj, i2 = solver()
            return b, s_minus, N, T, vneg, vj, i2
        if ring.segment == 2:
            b, s_minus, N, T,vneg_LSB, vj_LSB, i2_LSB,\
            vneg_MSB, vj_MSB, i2_MSB  = solver()
        return  b, s_minus, N, T,vneg_LSB, vj_LSB, i2_LSB,\
                vneg_MSB, vj_MSB, i2_MSB
