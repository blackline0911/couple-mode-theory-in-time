import utility
import sim 

class Heater:
    id = "Heater"
    def __init__(self,T_surround,
                 V,
                 R,
                 tau_th = 100, #ns 
                 ):
        self.V = V
        self.R = R
        self.T_surround = T_surround
        self.tau_th = tau_th*1e-9
        self.renew()

    def renew(self):
        self.P = self.V**2/self.R*1000 #mW
        