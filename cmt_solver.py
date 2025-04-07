# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
from scipy.integrate import *
import numpy as np
from cmath import *
from utility import *
from driver import driver
from ring import ring
from time_class import time

# relative solver tolerance 
rtol = 1e-14
# absolute solver tolerance
atol = 1e-20
# accuracy = atol + abs(y)*rtol

# method of solving algorithm
method = 'RK45'

     

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

    # normalize factors
    S0 = np.real(sqrt(sim.Pin))

    b0 = sqrt(t0)*S0

    
    # normalized dw/dlambda
    D_bar = -2*np.pi*c/ring.lambda0**2 * t0


#     """Since the length of solution function array may not be the same as t_eval argument we specified when the time in single solve_ivp is long.
#     Hence, we divide the time according to the Baud Rate, and solve coupled differential equation by each time segments. 
#     """
    b_record = np.array([])
    Q_record = np.array([])
    N_record = np.array([])
    b_init=0+1j*0
    Q_init=0
    N_init=0
    vg_in_cm = ring.vg*1e-4
    # notes: time_range argument should be slightly exclude the t_eval
    if time.mode == "voltage_drive":
        A = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*t0**2*S0**4
        def CMT(t_bar,eqs):
            b_bar , Q_pround ,N= eqs
            # b_bar , Q_pround = eqs
            voltage = driver.refering_v(t_bar)
            cj = driver.refering_Cj(voltage)

            f1 = 1j*2*np.pi*(ring.f_res_bar-sim.f_pround_bar)*b_bar- \
                (ring.tu_e_bar_total_inv + vg_in_cm*ring.alpha_linear*(1 + ring.TPA_ratio*t0*S0**2*abs(b_bar)**2\
                    + N*1e-5/(vg_in_cm*ring.alpha_linear) ) )*b_bar + \
                      ring.input_kappa *1 + \
                1j*D_bar*(-ring.me*1e-12/1e-6)*Q_pround*b_bar
            # f1 = 1j*2*np.pi*(ring.f_res_bar-sim.f_pround_bar)*b_bar- \
            #     (1/ring.tu_e_bar + vg_in_cm*ring.alpha_linear*(1 + ring.TPA_ratio*t0*S0**2*abs(b_bar)**2\
            #         + ring.FCA_ratio*t0**2*S0**4*abs(b_bar)**4) )*b_bar +\
            #           sqrt(2/ring.tu_e_bar) *1 + \
            #     1j*D_bar*(-ring.me*1e-12/1e-6)*Q_pround*b_bar

            f2 = (voltage/(driver.R * cj )*t0) - (1/( driver.R*cj ) )*Q_pround*t0

            f3 = -N/ring.tau_eff*(t0/1e-9) + A*abs(b_bar)**4
            return [ f1,f2 ,f3]
        for i in range(time.N):
            if i==0:
                sol =  solve_ivp(CMT ,[0 ,time.t_all_segment[0][-1]], [b_init, Q_init,N_init],
                    method=method,t_eval = time.t_all_segment[i] ,atol = atol,rtol = rtol)
                b = sol.y[0]
                b_record = np.append(b_record,sol.y[0])
                q = sol.y[1]
                Q_record = np.append(Q_record,sol.y[1])
                n = sol.y[2]
                N_record = np.append(N_record,sol.y[2])
            else:
                sol =  solve_ivp(CMT ,[time.t_all_segment[i-1][-1] ,time.t_all_segment[i][-1]], 
                                [b_init, Q_init,N_init],
                                method=method,
                                t_eval = np.append(  np.array([time.t_all_segment[i-1][-1]]), np.array(time.t_all_segment[i])  ),
                                atol = atol,rtol = rtol)
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

    if time.mode == "scan_frequency":
        A = ring.FCA_coeff/(ring.tau_eff*1e-9)*1e5*t0*t0**2*S0**4
        # TPA = ring.TPA_coeff*t0*S0**2
        # FCA = ring.FCA_coeff*t0**2*S0**4
        def CMT(t_bar,eqs):
            f_res_bar = ring.w_res(t_bar)
            b_bar , Q_pround, N = eqs
            voltage = driver.refering_v(t_bar)
            if driver.varying_cj:
                cj = driver.refering_Cj(voltage)
            else:
                cj = driver.Cj
            f1 = 1j*2*np.pi*(f_res_bar-sim.f_pround_bar)*b_bar \
                - (ring.tu_e_bar_total_inv + \
                   vg_in_cm*ring.alpha_linear*(1 + ring.TPA_ratio*t0*S0**2*abs(b_bar)**2\
                    + N*1e-5/(vg_in_cm*ring.alpha_linear) ) )*b_bar + \
                ring.input_kappa *1 + \
                1j*D_bar*(-ring.me*1e-12/1e-6)*Q_pround*b_bar
            # f1 = 1j*2*np.pi*(f_res_bar-sim.f_pround_bar)*b_bar \
            #     - (ring.tu_e_bar_total_inv + \
            #        vg_in_cm*ring.alpha_linear*(1 + ring.TPA_ratio*t0*S0**2*abs(b_bar)**2\
            #         + ring.FCA_ratio*t0**2*S0**4*abs(b_bar)**4) )*b_bar + \
            #     ring.input_kappa *1 + \
            #     1j*D_bar*(-ring.me*1e-12/1e-6)*Q_pround*b_bar
            # f1 = 1j*2*np.pi*(f_res_bar-sim.f_pround_bar)*b_bar \
            #     - (ring.tu_e_bar_total_inv + \
            #        vg_in_cm*( ring.alpha_linear + TPA*abs(b_bar)**2\
            #         + FCA*abs(b_bar)**4) )*b_bar + \
            #     ring.input_kappa *1 + \
            #     1j*D_bar*(-ring.me*1e-12/1e-6)*Q_pround*b_bar

            f2 = (voltage/(driver.R * cj )*t0) - (1/( driver.R*cj ) )*Q_pround*t0

            f3 = -N/ring.tau_eff*(t0/1e-9) + A*abs(b_bar)**4
            return [ f1,f2,f3]
        sol = solve_ivp(CMT ,[0,time.t_max+time.buffer], [b_init, Q_init,N_init],method=method,t_eval = time.t_total,atol = atol,rtol = rtol)
        b_bar = sol.y[0]
        Q_bar = sol.y[1]
        N_bar = sol.y[2]
    s_minus_bar = (1-ring.input_kappa*b_bar)

    return b_bar*b0, Q_bar*driver.Cj , s_minus_bar*S0    ,N_bar/(ring.sigma_FCA*1e-17)*1e-5   
