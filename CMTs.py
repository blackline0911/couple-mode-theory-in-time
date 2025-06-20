from utility import t0, c
from driver import driver
from ring import ring
from sim import simulation
from Heater import Heater
from numpy import real
from time_class import time

def CMT_large_signal_single_segment(t_bar,eqs,driver:driver,ring:ring,Heater:Heater,SPM=None,TPA=None,FCA=None,sim=None,T_args=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    b_bar  ,N_bar,delta_T, v_neg, vj, i2= eqs
    voltage = real(driver.refering_v(t_bar))
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    alpha_linear = ring.alpha(vj)
    cj = driver.Cj_V(vj,driver.a,driver.b)
    cj_bar = cj/t0
    dlambda = ring.lambda0/ring.ng*( ring.neff(vj) - ring.neff(0))*(ring.L_active[0]/ring.L)
    Rs = driver.Rs
    Cox_bar = driver.Cox/t0
    Rsi = driver.Rsi
    Cp_bar = driver.Cp/t0

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    # Coupled Transmission Line Equations
    dvneg_dt, dvj_dt, di2_dt = driver.TML(v_neg, vj, voltage, dvpos_dt, i2, \
                                          driver.Z0, Rs, cj_bar, Cox_bar, Rsi, Cp_bar)
    
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,delta_T*ring.self_heating_factor,ring.f_res_bar,alpha_linear,TPA,SPM,T_args,dlambda,Heater)
    
    dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)

    dT_dt = Heater.T_rate_equation(b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim)

    return [ da_dt ,dN_dt, dT_dt, dvneg_dt, dvj_dt, di2_dt]

def CMT_large_signal_dual_segment(t_bar,eqs,driver:driver,SPM=None,TPA=None,FCA=None,ring:ring=None,sim:simulation=None,Heater=None,T_args=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    b_bar  ,N_bar,delta_T, v_neg_LSB, v_neg_MSB, vj_LSB, vj_MSB, i2_LSB, i2_MSB= eqs
    
    voltage = real(driver.refering_v(t_bar))
    voltage_shifted = real(driver.refering_v(t_bar-time.T_normalized*driver.shift_bit))/2
    v_combined = voltage + voltage_shifted
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    dvpos_dt_shifted = driver.refering_dv_dt(voltage,t_bar-time.T_normalized*driver.shift_bit)/2
    cj_LSB_bar = driver.Cj_V(vj_LSB,driver.a,driver.b)/t0
    cj_MSB_bar = driver.Cj_V(vj_MSB,driver.c,driver.d)/t0
    alpha_linear = ring.alpha(v_combined)
    dlambda_LSB = ring.lambda0/ring.ng*( ring.neff(vj_LSB) - ring.neff(0))*(ring.L_active[0]/ring.L)
    dlambda_MSB = ring.lambda0/ring.ng*( ring.neff(vj_MSB) - ring.neff(0))*(ring.L_active[1]/ring.L)
    Rs_LSB = driver.Rs[0]
    Rs_MSB = driver.Rs[1]
    Cox_LSB_bar = driver.Cox[0]/t0
    Cox_MSB_bar = driver.Cox[1]/t0
    Rsi_LSB = driver.Rsi[0]
    Rsi_MSB = driver.Rsi[1]
    Cp_LSB_bar = driver.Cp[0]/t0
    Cp_MSB_bar = driver.Cp[1]/t0

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # Coupled Transmission Line Equations
    dvneg_dt_LSB, dvj_dt_LSB, di2_dt_LSB = driver.TML(v_neg_LSB, vj_LSB, voltage, dvpos_dt, i2_LSB, \
                                        driver.Z0,Rs_LSB, cj_LSB_bar,Cox_LSB_bar,Rsi_LSB,Cp_LSB_bar)
    dvneg_dt_MSB, dvj_dt_MSB, di2_dt_MSB = driver.TML(v_neg_MSB, vj_MSB, voltage, dvpos_dt_shifted, i2_MSB, \
                                        driver.Z0,Rs_MSB, cj_MSB_bar,Cox_MSB_bar,Rsi_MSB,Cp_MSB_bar)
        
    # Couple mode Equations in time domain
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,\
                     delta_T*ring.self_heating_factor,\
                    ring.f_res_bar,alpha_linear,\
                    TPA,SPM,T_args,\
                    dlambda_LSB + dlambda_MSB,\
                    Heater)

    dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)

    dT_dt = Heater.T_rate_equation(b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim)
    
    return [ da_dt ,dN_dt ,dT_dt , dvneg_dt_LSB, dvj_dt_LSB, di2_dt_LSB, dvneg_dt_MSB, dvj_dt_MSB, di2_dt_MSB]

def CMT_small_signal_single_segment(t_bar,eqs,driver:driver,SPM=None,TPA=None,FCA=None,ring:ring=None,sim:simulation=None,Heater=None,T_args=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    b_bar  ,N_bar,delta_T, v_neg, vj, i2= eqs

    voltage = real(driver.refering_v(t_bar))
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    cj = driver.Cj_V(driver.v_bias,driver.a,driver.b)
    cj_bar = cj/t0
    alpha_linear = ring.alpha(vj)
    dlambda = ring.lambda0/ring.ng*( ring.neff(vj) - ring.neff(0))*(ring.L_active[0]/ring.L)
    Rs = driver.Rs[0]
    Cox_bar = driver.Cox[0]/t0
    Rsi = driver.Rsi[0]
    Cp_bar = driver.Cp[0]/t0

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # Coupled Transmission Line Equations
    dvneg_dt, dvj_dt, di2_dt = driver.TML(v_neg, vj, voltage, dvpos_dt, i2, \
                                        driver.Z0,Rs, cj_bar,Cox_bar,Rsi,Cp_bar)
        
    # Couple mode Equations in time domain
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,delta_T*ring.self_heating_factor,ring.f_res_bar,alpha_linear,TPA,SPM,T_args,dlambda,Heater)

    dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)

    dT_dt = Heater.T_rate_equation(b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim)
    
    return [ da_dt ,dN_dt ,dT_dt , dvneg_dt, dvj_dt, di2_dt]
    
def CMT_small_signal_dual_segment(t_bar,eqs,driver:driver,SPM=None,TPA=None,FCA=None,ring:ring=None,sim:simulation=None,Heater=None,T_args=None):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    b_bar  ,N_bar,delta_T, v_neg_LSB, v_neg_MSB, vj_LSB, vj_MSB, i2_LSB, i2_MSB= eqs
    
    voltage = real(driver.refering_v(t_bar))
    voltage_shifted = real(driver.refering_v(t_bar-time.T_normalized*driver.shift_bit))/2
    v_combined = voltage + voltage_shifted
    dvpos_dt = driver.refering_dv_dt(voltage,t_bar)
    dvpos_dt_shifted = driver.refering_dv_dt(voltage,t_bar-time.T_normalized*driver.shift_bit)/2
    cj_LSB_bar = driver.Cj_V(driver.v_bias,driver.a,driver.b)/t0
    cj_MSB_bar = driver.Cj_V(driver.v_bias,driver.c,driver.d)/t0
    alpha_linear = ring.alpha(v_combined)
    dlambda_LSB = ring.lambda0/ring.ng*( ring.neff(vj_LSB) - ring.neff(0))*(ring.L_active[0]/ring.L)
    dlambda_MSB = ring.lambda0/ring.ng*( ring.neff(vj_MSB) - ring.neff(0))*(ring.L_active[1]/ring.L)
    Rs_LSB = driver.Rs[0]
    Rs_MSB = driver.Rs[1]
    Cox_LSB_bar = driver.Cox[0]/t0
    Cox_MSB_bar = driver.Cox[1]/t0
    Rsi_LSB = driver.Rsi[0]
    Rsi_MSB = driver.Rsi[1]
    Cp_LSB_bar = driver.Cp[0]/t0
    Cp_MSB_bar = driver.Cp[1]/t0

    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    # Coupled Transmission Line Equations
    dvneg_dt_LSB, dvj_dt_LSB, di2_dt_LSB = driver.TML(v_neg_LSB, vj_LSB, voltage, dvpos_dt, i2_LSB, \
                                        driver.Z0,Rs_LSB, cj_LSB_bar,Cox_LSB_bar,Rsi_LSB,Cp_LSB_bar)
    dvneg_dt_MSB, dvj_dt_MSB, di2_dt_MSB = driver.TML(v_neg_MSB, vj_MSB, voltage, dvpos_dt_shifted, i2_MSB, \
                                        driver.Z0,Rs_MSB, cj_MSB_bar,Cox_MSB_bar,Rsi_MSB,Cp_MSB_bar)
        
    # Couple mode Equations in time domain
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,\
                     delta_T*ring.self_heating_factor,\
                    ring.f_res_bar,alpha_linear,\
                    TPA,SPM,T_args,\
                    dlambda_LSB + dlambda_MSB,\
                    Heater)

    dN_dt = ring.FC_rate_equation(b_bar,N_bar,FCA,ring.tau_eff)

    dT_dt = Heater.T_rate_equation(b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim)
    
    return [ da_dt ,dN_dt ,dT_dt , dvneg_dt_LSB, dvj_dt_LSB, di2_dt_LSB, dvneg_dt_MSB, dvj_dt_MSB, di2_dt_MSB]
    
def CMT_scan_frequency(t_bar,eqs,SPM,TPA,FCA,T_args,ring:ring,sim:simulation,driver:driver,Heater:Heater):
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    f_res_bar = ring.w_res(t_bar)
    b_bar  = eqs
    voltage = driver.v_bias
    alpha_linear = ring.alpha(voltage)
    df = -sim.f_pround_bar/ring.ng*( ring.neff(voltage) - ring.neff(0))*(ring.L_active[0]/ring.L)
    
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    # 不要忘記這裡的tau_o吸收是energy的吸收，不是amplitude的

    N_bar = FCA* abs(b_bar)**4

    delta_T = abs(b_bar)**2 * T_args[1] * (\
        \
        +ring.vg_in_cm*alpha_linear \
        \
        + ring.vg_in_cm *TPA*abs(b_bar)**2 \
        \
        + ring.vg_in_cm * N_bar*1e-5 \
        ) * Heater.tau_th/t0
    
    da_dt = ring.CMT(sim.f_pround_bar,b_bar,N_bar,delta_T*ring.self_heating_factor,f_res_bar,alpha_linear,TPA,SPM,T_args,\
                     (-c*1e-12/sim.f_pround_bar**2)*df , Heater)
    
    return [ da_dt]