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




def CMT_large_signal(t_bar,eqs,SPM=None,FCA=None,vg_in_cm=None,ring=None,sim=None,driver=None):
    b_bar , Q_pround ,N_bar= eqs
    voltage = np.real(driver.refering_v(t_bar))
    f1 = 1j*2*np.pi*(ring.f_res_bar-sim.f_pround_bar + \
                     SPM*abs(b_bar)**2)*b_bar- \
    \
        (ring.tu_e_bar_total_inv + \
            ring.vg_in_cm*ring.alpha_linear*(1 + ring.TPA_ratio*(sim.b0)**2*abs(b_bar)**2\
    \
            + N_bar*1e-5/(ring.vg_in_cm*ring.alpha_linear) ) )*b_bar + \
                ring.input_kappa *1 + \
    \
        1j*ring.D_bar*(-ring.me*1e-12/1e-6)*driver.V_Q(Q_pround)*b_bar    
        
    f2 = (voltage/(driver.R * driver.cj_normalizing )*t0) \
        - (1/( driver.R ) )*driver.V_Q(Q_pround)*t0/driver.cj_normalizing

    f3 = -N_bar/ring.tau_eff*(t0/1e-9) + FCA*abs(b_bar)**4

    return [ f1,f2 ,f3]

def CMT_small_signal(t_bar,eqs,SPM,FCA,ring,sim,driver):
    b_bar , Q_pround ,N_bar= eqs
    voltage = np.real(driver.refering_v(t_bar))
    cj = driver.Cj
    f1 = 1j*2*np.pi*(ring.f_res_bar - sim.f_pround_bar + \
                      SPM*abs(b_bar)**2)*b_bar- \
    \
        (ring.tu_e_bar_total_inv + \
            ring.vg_in_cm*ring.alpha_linear*(1 + ring.TPA_ratio*(sim.b0)**2*abs(b_bar)**2\
    \
            + N_bar*1e-5/(ring.vg_in_cm*ring.alpha_linear) ) )*b_bar + \
                ring.input_kappa *1 + \
    \
        1j*ring.D_bar*(-ring.me*1e-12/1e-6)*driver.cj_normalizing* (Q_pround/cj)*b_bar

    f2 = (voltage/(driver.R * driver.cj_normalizing )*t0) \
        - (1/( driver.R*cj)*t0 )*Q_pround

    f3 = -N_bar/ring.tau_eff*(t0/1e-9) + FCA*abs(b_bar)**4

    return [ f1,f2 ,f3]

def CMT_scan_frequency(t_bar,eqs,SPM,FCA,ring,sim,driver):
    f_res_bar = ring.w_res(t_bar)
    b_bar , Q_pround, N_bar = eqs
    voltage = driver.v_bias
    cj = driver.Cj_V(voltage)
    f1 = 1j*2*np.pi*(f_res_bar-sim.f_pround_bar)*b_bar \
        - (ring.tu_e_bar_total_inv + \
            ring.vg_in_cm*ring.alpha_linear*(1 + ring.TPA_ratio*(sim.b0)**2*abs(b_bar)**2\
            + N_bar*1e-5/(ring.vg_in_cm*ring.alpha_linear) ) )*b_bar + \
        ring.input_kappa *1 + \
        1j*ring.D_bar*(-ring.me*1e-12/1e-6)*Q_pround*b_bar +\
        SPM*abs(b_bar)**2
    
    f2 = (voltage/(driver.R * cj )*t0) - (1/( driver.R*cj ) )*Q_pround*t0

    f3 = -N_bar/ring.tau_eff*(t0/1e-9) + FCA*abs(b_bar)**4
    return [ f1,f2,f3]

# DeFine Solving process My Different simulation modes
def CMT_voltage_driving(sim,ring, 
                        driver,
                        time):
    method_dict = {
        'large_signal' : CMT_large_signal,
        'small_signal' : CMT_small_signal,
    }
    b_record = np.array([])
    Q_record = np.array([])
    N_record = np.array([])
    b_init=0+1j*0
    Q_init=0
    N_init=0
    # notes: time_range argument should be slightly exclude the t_eval
    FCA = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*(sim.b0)**4
    SPM = ring.dw_SPM_coeff * (sim.b0)**2
    ode_func = partial(method_dict[driver.method],
                       SPM=SPM,FCA=FCA,
                       ring = ring,sim=sim,
                       driver=driver)
    sol =  solve_ivp(ode_func ,[0 ,time.t_all_segment[0][-1]], 
                            [b_init, Q_init,N_init],
                            method=sim.algorithm,
                            t_eval = time.t_all_segment[0] ,
                            atol = atol,rtol = rtol,)
    b = sol.y[0]
    b_record = np.append(b_record,sol.y[0])
    q = sol.y[1]
    Q_record = np.append(Q_record,sol.y[1])
    n = sol.y[2]
    N_record = np.append(N_record,sol.y[2])
    b_init = sol.y[0][-1]
    Q_init = sol.y[1][-1]
    N_init = sol.y[2][-1]
    
    for i in range(1,time.N):
        sol =  solve_ivp(ode_func ,
                        [time.t_all_segment[i-1][-1] ,\
                        time.t_all_segment[i][-1]], 
                        [b_init, Q_init, N_init],
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
        b_init = sol.y[0][-1]
        Q_init = sol.y[1][-1]
        N_init = sol.y[2][-1]
    b_bar = b_record
    Q_bar = Q_record
    N_bar = N_record
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*sim.b0, Q_bar*driver.cj_normalizing , s_minus_bar*sim.S0    ,N_bar/(ring.sigma_FCA*1e-17)*1e-5   

def solve_scan_frequency(sim,ring, 
                        driver,
                        time):
    SPM = ring.dw_SPM_coeff *(sim.b0)**2
    FCA = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*(sim.b0)**4
    b_init = 0+1j*0
    Q_init = 0
    N_init = 0
    sol = solve_ivp(CMT_scan_frequency ,
                    [0,time.t_max+time.buffer], 
                    [b_init, Q_init,N_init],
                    method=sim.algorithm,
                    t_eval = time.t_total,atol = atol,rtol = rtol,
                    args=(SPM,FCA,ring,sim,driver))
    b_bar = sol.y[0]
    Q_bar = sol.y[1]
    N_bar = sol.y[2]
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*sim.b0, Q_bar*driver.cj_normalizing , s_minus_bar*sim.S0    ,N_bar/(ring.sigma_FCA*1e-17)*1e-5

def solving(sim,
            ring, 
            driver,
            time,
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

    solver = partial(mode_dict[sim.mode],sim,ring,driver)
    b, Q, s_minus, N = solver(time=time)
    return b, Q, s_minus   ,N
