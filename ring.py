from utility import *
import numpy as np
from cmath import *
from time_class import *

class ring():
    def __init__(self,radius:float, 
                 neff:float,
                 ng:float,
                 gamma:float,
                 alpha:float,
                 me:float,
                 FCA_coeff_ratio = 1e-3):
        """
        input argments:
        radius : radius of ring (micron)
        neff : effective phase index of waveguide in ring
        ng: group index of waveguide in ring
        gamma:couple through coefficient of the bus waveguide of ring
        alpha:round trip amplitude attenuation coefficient (one for no loss)
        me : modulation efficiency (pm/reverse voltage)

        """
        self.L=2*np.pi*radius
        self.neff=neff
        self.ng=ng
        self.gamma=gamma
        self.alpha = alpha
        self.me = me
        self.kappa = (1-self.gamma**2)**0.5
        self.tu_e_bar = -self.L*self.ng/(c*log(sqrt(1-self.kappa**2)))/t0
        self.tu_o_bar = -self.L*self.ng/(c*log(alpha))/t0
        self.tu_t_bar = (1/self.tu_e_bar+1/self.tu_o_bar)**(-1)
        self.f_res_bar = 51*c/(self.neff*self.L)*t0
        self.lambda0 = c*t0/self.f_res_bar
        self.FCA_coeff_ratio = FCA_coeff_ratio
        self.Q = np.real( ( (1/self.tu_e_bar + 1/self.tu_o_bar  )/(self.f_res_bar*np.pi) )**(-1) )

    def scan_frequency(self,wl_min:float,wl_max:float,lambda_incident,time):
        """
        Call this function when you wants to perform frequency scan of the ring.
        I will perform resonant frequency scan, instead of incident frequency scan
        input:
        wl_min: minimum scanning wavelength (um)
        wl_max: maximum scanning wavelength (um)
        lambda_incident: incident laser wavelength (um)
        """
        self.time = time
        self.f_max_bar = c/wl_min*t0
        self.f_min_bar = c/wl_max*t0
        print("Q factor of ring is ",self.Q)
        print("scanning frequency range for better accuracy is ",self.f_res_bar/self.Q)
        print("Current frequency range is ",(self.f_max_bar-self.f_min_bar))
        if  self.Q >self.f_res_bar/(self.f_max_bar-self.f_min_bar):
            print("Warning : This Transfer function result may not be accurate.")

        self.f_in_bar = c/lambda_incident*t0
        self.f_res = np.zeros(len(self.time.t_total))

    def w_res(self,t):
        """
        return the resonant frequency according to specified time t
        """
        if isinstance(t,np.ndarray):
            ans = -(self.f_max_bar-self.f_min_bar)/self.time.t_max*(t-self.time.buffer) + self.f_res_bar+ \
                (self.f_in_bar-self.f_min_bar)
            return np.where(t<self.time.buffer, self.f_res_bar+(self.f_in_bar-self.f_min_bar), ans)
        else:
            if t<(self.time.buffer):
                return self.f_res_bar+(self.f_in_bar-self.f_min_bar)
            else:
                return -(self.f_max_bar-self.f_min_bar)/self.time.t_max*(t-self.time.buffer)+ self.f_res_bar+ \
                    (self.f_in_bar-self.f_min_bar)
