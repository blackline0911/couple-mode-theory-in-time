from utility import *
import numpy as np
from cmath import *
from time_class import *
from sim import simulation

class ring(simulation):
    id='ring'
    def __init__(self,L:float, 
                 alpha:float,
                 me:float,
                 cross_section:float,
                 lambda_incident:float,
                 *gamma,
                 neff=None,
                 band = "C",
                 lambda0=None,
                 FSR = None,
                 ng = None,
                 FSR_shift = 0,
                 beta_TPA = 5,  #1e-13 (cm/mW)
                 tau_eff = 20,  # ns
                 sigma_FCA = 1.04,  #1e-17 cm^2
                 FCA_fit_factor = 1,
                 TPA_fit_factor = 1,
                 input_port = 1,
                ):
        """
        input argments:\n
        \tL : Length of cavity (micron)\n
        \tneff : effective phase index of waveguide in ring\n
        \tng: group index of waveguide in ring\n
        \tgamma:amplitude couple through coefficient of the bus waveguide of ring\n
        \t(gamma can input multiple values for multi-port ring (or add-drop ring), \n
        \tthe first gamma is port of input Laser by default (can specified by input port ))\n
        \talpha:round trip amplitude attenuation coefficient (one for no loss)\n
        \tme : modulation efficiency (pm/reverse voltage)\n
        \tband : specify operating wavelength, "C" for 1550nm, "O" for 1310nm\n
        \tcross_section : mode area in waveguide (um^2)\n
        \tlambda0 : specify resonant wavelength ,calculate automatically if ungiven\n
        \tFSR_shift : specify the peak of simulation shifted by the peak near 1310 nm or 1550 nm\n
        \t(positive means simulate longer resonant wavelength, negative means simulate shorter resonant wavelength)

        \tNonLinear absorption parameter:\n
        \tbeta_TPA : Two Photon Absorption coefficient (cm/mW, normalized by 1e-13)\n
        \ttau_eff : effective carrier life time (ns)\n
        \tsigma_FCA :Free carrier absorption area (cm^2, normalized by 1e-17)\n
        """
        self.L=L
        self.ng=ng
        self.gamma=np.real(gamma)
        assert (not np.imag(g)==0.0 for g in gamma) , "\nCoupling coefficient shall not over one\n"
        self.alpha = alpha
        self.me = me
        self.band = band
        self.support_band = {"C","O"}
        assert self.band in self.support_band, f"The simulation only support {self.support_band} bands \n"
        if (lambda0==None):
            self.neff=neff
            m = self.find_reference_res_wavelength()
        else:
            self.lambda0 = lambda0
            self.lambda0_reference = lambda0
        if FSR==None:
            self.ng = ng
            self.vg = c/self.ng*t0
        elif ng==None:
            self.FSR = FSR
            self.ng = self.lambda0_reference**2/(self.FSR*self.L)
            self.vg = c/self.ng*t0
        else:
            assert False, "\nDo not specify FSR and ng at the same time \nYou have to specify one of them.\n"
        
        if (lambda0==None):
            self.FSR_shift = FSR_shift
            self.f_res_bar = self.find_accurate_res_frequency(m)
            self.lambda0 = c*t0/self.f_res_bar
        else:
            self.f_res_bar = c*t0/self.lambda0

        self.round_trip_time = self.L*self.ng/(c*t0)
        self.cross_section =  cross_section #um^2
        self.kappa = (1-self.gamma**2)**0.5
        self.tu_e_bar = -self.L*self.ng/(c*t0*np.log((1-self.kappa**2)**0.5))
        self.tu_o_bar = -self.L*self.ng/(c*t0*np.log(alpha))
        self.kappa_total = 0
        self.tu_e_bar_total_inv = 0
        for tu_e_bar in self.tu_e_bar:
            self.kappa_total += (2/tu_e_bar)**0.5 
            self.tu_e_bar_total_inv += 1/tu_e_bar
        assert input_port>0 , "\nInput port should start from one\n"
        
        self.input_kappa = (2/self.tu_e_bar[input_port-1])**0.5 
        self.alpha_linear = np.real(1/(self.vg*1e-4*self.tu_o_bar))
        self.tu_t_bar = (self.tu_e_bar_total_inv+1/self.tu_o_bar)**(-1)
        self.beta_TPA = beta_TPA  #1e-13 (cm/mW)
        self.tau_eff = tau_eff      #1e-9  s,
        self.sigma_FCA = sigma_FCA  #1e-17 cm^2
        super().__init__()
        self.lambda_incident = lambda_incident
        self.photon_energy = h*c/self.lambda_incident*1000/t0  #mJ
        self.FCA_fit_factor = FCA_fit_factor
        self.TPA_fit_factor = TPA_fit_factor
        # Note. photon energy is in unit. mJ, and normalized by t0

        
        self.Q = np.real( ( (self.tu_e_bar_total_inv + 1/self.tu_o_bar  )/(self.f_res_bar*np.pi) )**(-1) )
        # self.Q = np.real( ( (1/self.tu_e_bar + 1/self.tu_o_bar  )/(self.f_res_bar*np.pi) )**(-1) )
        

        self.handle_nonlinear()

    def scan_frequency(self,wl_start:float,wl_end:float,time):
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
        self.f_in_bar = c/self.lambda_incident*t0
        self.f_res = np.zeros(len(self.time.t_total))
    
    def handle_nonlinear(self):
        #unit of TPA coeff : 1/(cm* mJ) 
        #  beta_TPA in 1e-9 order
        self.TPA_coeff = self.beta_TPA/(self.round_trip_time*self.cross_section)
        self.TPA_coeff_order = np.log10(1e-13/(1e-12*1e-8)) + np.log10(self.TPA_coeff)//1
        self.TPA_coeff = self.TPA_coeff / ( 10**(np.log10(self.TPA_coeff)//1) )
        self.TPA_coeff = self.TPA_coeff*10**( self.TPA_coeff_order)
        self.TPA_coeff = self.TPA_coeff*self.TPA_fit_factor
        self.TPA_ratio = self.TPA_coeff/self.alpha_linear
        # Note. The unit TPA_coeff here is  1/(cm*mJ), not 1/cm
     
        self.FCA_coeff = self.beta_TPA*self.tau_eff*self.sigma_FCA/ \
            (2*self.photon_energy*self.round_trip_time**2*self.cross_section**2) 
        self.FCA_coeff_order = np.log10(1e-13*1e-9*1e-17/( 1e-12*(1e-12)**2*(1e-4)**4 ) ) + np.log10(self.FCA_coeff)//1
        self.FCA_coeff = self.FCA_coeff / ( 10**(np.log10(self.FCA_coeff)//1) )
        self.FCA_coeff = self.FCA_coeff * 10**(self.FCA_coeff_order)
        self.FCA_coeff = self.FCA_coeff*self.FCA_fit_factor
        self.FCA_ratio = self.FCA_coeff/self.alpha_linear
        # Note. The unit FCA_coeff here is  1/(cm*mJ^2), not 1/cm
        

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

    def find_reference_res_wavelength(self):
        if self.band=="C":
            f_band = c/1.55*t0
            m = self.neff*self.L/1.55
        if self.band=="O":
            f_band = c/1.31*t0
            m = self.neff*self.L/1.31
        a = np.ceil(m)*c/self.neff/self.L*t0 - f_band
        b = f_band - np.floor(m)*c/self.neff/self.L*t0
        if a<b:
            m= np.ceil(m)
        if a>b:
            m=np.floor(m)
        self.lambda0_reference = (self.neff*self.L)/m
        
        return m
    
    def find_accurate_res_frequency(self,m):
        f0_v0 = c*t0/self.lambda0_reference
        m_pround = m - self.FSR_shift
        # Assume the specified neff (self.neff) is index of resonant wavelength near 1550 or 1310 nm 
        # solving index at shifted resonance frequency
        A = 1
        B = -(2*self.neff-self.ng)
        C = -(self.ng-self.neff)/(f0_v0*self.L)*m_pround*c*t0
        n = np.real((-B+sqrt(B**2-4*A*C))/(2*A))
        assert np.imag((-B+sqrt(B**2-4*A*C))/(2*A))==0.0 , "\nimaginary part of index is supposed to be zero\n"
        
        return m_pround*c/n/self.L*t0

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