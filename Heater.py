import utility
import sim 

class Heater:
    id = "Heater"
    def __init__(self,T_surround,
                 V,
                 R):
        self.V = V
        self.R = R
        self.T_surround = T_surround
        self.renew()

    def renew(self):
        self.P = self.V**2/self.R*1000 #mW
        