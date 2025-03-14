from utility import *
import numpy as np
from cmath import *
from time_class import *
from sim import simulation

class ring():
    def __init__(self,L:float, 
                 neff:float,
                 ng:float,
                 gamma:float,
                 alpha:float,
                 me:float,
                 FCA_coeff_ratio = 0):
        """
        input argments:
        radius : radius of ring (micron)
        neff : effective phase index of waveguide in ring
        ng: group index of waveguide in ring
        gamma:couple through coefficient of the bus waveguide of ring
        alpha:round trip amplitude attenuation coefficient (one for no loss)
        me : modulation efficiency (pm/reverse voltage)

        """
        self.L=L
        self.neff=neff
        self.ng=ng
        self.gamma=gamma
        self.alpha = alpha
        self.me = me
        self.kappa = (1-self.gamma**2)**0.5
        self.tu_e_bar = -self.L*self.ng/(c*log(sqrt(1-self.kappa**2)))/t0
        self.tu_o_bar = -self.L*self.ng/(c*log(alpha))/t0
        self.tu_t_bar = (1/self.tu_e_bar+1/self.tu_o_bar)**(-1)
        # self.f_res_bar = 51*c/(self.neff*self.L)*t0
        self.f_res_bar = self.find_res_frequency()
        self.lambda0 = c*t0/self.f_res_bar
        self.FCA_coeff_ratio = FCA_coeff_ratio
        self.Q = np.real( ( (1/self.tu_e_bar + 1/self.tu_o_bar  )/(self.f_res_bar*np.pi) )**(-1) )

    def scan_frequency(self,wl_start:float,wl_end:float,lambda_incident,time):
        """
        Call this function when you wants to perform frequency scan of the ring.
        I will perform resonant frequency scan, instead of incident frequency scan
        input:
        wl_min: start scanning wavelength (um)
        wl_max: Ending scanning wavelength (um)
        lambda_incident: incident laser wavelength (um)
        """
        self.time = time
        self.f_start_bar = c/wl_start*t0
        # print(self.f_start_bar)
        self.f_end_bar = c/wl_end*t0
        # print(self.f_end_bar)
        self.f_in_bar = c/lambda_incident*t0
        self.f_res = np.zeros(len(self.time.t_total))

    def w_res(self,t):
        """
        return the resonant frequency according to specified time t
        """
        if isinstance(t,np.ndarray):
            ans = -(self.f_end_bar-self.f_start_bar)/self.time.t_max*(t-self.time.buffer) + self.f_res_bar+ (self.f_in_bar-self.f_start_bar)
            return np.where(t<self.time.buffer, self.f_res_bar+(self.f_in_bar-self.f_start_bar), ans)
        else:
            if t<(self.time.buffer):
                return self.f_res_bar+(self.f_in_bar-self.f_start_bar)
            else:
                return -(self.f_end_bar-self.f_start_bar)/self.time.t_max*(t-self.time.buffer)+ self.f_res_bar+ (self.f_in_bar-self.f_start_bar)

    def find_res_frequency(self):
        f_1550 = c/1.55*t0
        m = self.neff*self.L/1.55
        a = np.ceil(m)*c/self.neff/self.L*t0 - f_1550
        b = f_1550 - np.floor(m)*c/self.neff/self.L*t0
       
        if a<b:
            m= np.ceil(m)
        if a>b:
            m=np.floor(m)
        return m*c/self.neff/self.L*t0

class Transfer_function():
    def __init__(self,ring, time):
        self.ring = ring
        self.time = time
    def mapping(self,data):
        """
        Convert the frequency scanning time domain signal to transfer function of ring
        """
        i = int( np.argwhere(self.time.t_total==self.time.buffer)[0] )
        f_res_bar = self.ring.w_res(self.time.t_total)
        L = len(f_res_bar)
        wl = c/(self.ring.f_res_bar + self.ring.f_in_bar -f_res_bar)*t0
        wl = wl[i:L-1]
        data = data[i:L-1]
        return wl, data