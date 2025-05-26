from utility import *
import sim 

class Heater:
    id = "Heater"
    def __init__(self,T_surround,
                 V,
                 R,
                 tau_th = 100, #ns 
                 Aeff = 0.22*0.5,
                 ):
        self.V = V
        self.R = R
        self.T_surround = T_surround
        self.tau_th = tau_th*1e-9
        self.renew()

    def renew(self):
        self.P = self.V**2/self.R*1000 #mW

    def T_rate_equation(self,b_bar,N_bar,delta_T,T_args,alpha_linear,TPA,ring,sim):
        dT_dt = abs(b_bar)**2*T_args[1] * (\
        \
        +ring.vg_in_cm*alpha_linear \
        \
        + ring.vg_in_cm *TPA*abs(b_bar)**2 \
        \
        + ring.vg_in_cm * N_bar*1e-5 \
        ) \
        \
        - t0/self.tau_th*delta_T
        return dT_dt
        