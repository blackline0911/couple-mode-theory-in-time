from utility import *
import numpy as np
from cmath import *
from time_class import *
from sim import simulation
from scipy.optimize import curve_fit

class ring(simulation):
    id='ring'
    def __init__(self,L:float, 
                 L_active:float,
                 alpha:float,
                #  me:float,
                 gamma:np.ndarray,
                 cross_section:float,
                 lambda_incident:float,
                 neff:float,
                 band = "C",
                 lambda0=None,
                 FSR = None,
                 ng = None,
                 FSR_shift = 0,
                 beta_TPA = 5,  #1e-13 (cm/mW)
                 tau_eff = 20,  # ns
                 sigma_FCA = 1.04,  #1e-17 cm^2
                 eta_h = 2.99,
                 HE = 254.3,
                 Akerr = 0.204, # um^2
                 Atpa = 0.1289, # um^2
                 Afca = 0.116, # um^2
                 n2 = 11e-9, #um^2/mW
                 FCA_fit_factor = 1,
                 TPA_fit_factor = 1,
                 SPM_fit_factor = 1,
                 input_port = 1,
                ):
        """
        input argments:\n
        \tL : Length of cavity (micron)\n
        \tneff : effective phase index of waveguide in ring\n
        \tng: group index of waveguide in ring\n
        \tgamma:AMPLITUDE couple through coefficient of the bus waveguide of ring\n
        \t(gamma can input multiple values for multi-port ring (or add-drop ring), \n
        \tthe first gamma is port of input Laser by default (can specified by input port ))\n
        \talpha:round trip AMPLITUDE attenuation coefficient (one for no loss)\n
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
        \teta_h : Normalized heater temperture (oC/mW) (Temperture in waveguide when applying 1mW electric power on Heater)\n
        \tHE : Heater efficiency of ring (pm/mW)
        """
        self.L=L
        self.L_active = L_active
        self.ng=ng
        self.gamma=np.real(gamma)
        assert (not np.imag(g)==0.0 for g in gamma) , "\nCoupling coefficient shall not over one\n"
        self.alpha = alpha
        # self.me = me
        self.band = band
        self.support_band = {"C","O"}

        # ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        # ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        self.eta_h = eta_h
        self.HE = HE 
        self.dlambda_dT = self.eta_h / (self.HE*1e-6)
        
        # //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        # /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        assert self.band in self.support_band, f"The simulation only support {self.support_band} bands \n"
        if (lambda0==None):
            self.neff=neff
            m = self.find_reference_res_wavelength()
        else:
            self.lambda0 = lambda0
            self.lambda0_reference = lambda0
        if FSR==None:
            self.ng = ng
            self.vg = c/self.ng*t0  #um/ps
        elif ng==None:
            self.FSR = FSR
            self.ng = self.lambda0_reference**2/(self.FSR*self.L)
            self.ng = 2.588
            self.vg = c/self.ng*t0  #um/ps
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
        self.kappa = (1-self.gamma**2)**0.5     #AMPLITUDE Coupling ratio

        #AMPLITUDE photon life time
        self.tu_e_bar = -self.L*self.ng/(c*t0*np.log(self.gamma))
        self.tu_o_bar = self.ng/((c*1e-4)*t0*alpha(0)/2)

        self.kappa_total = 0
        self.tu_e_bar_total_inv = 0
        for tu_e_bar in self.tu_e_bar:
            self.kappa_total += (1/tu_e_bar)**0.5 
            self.tu_e_bar_total_inv += 1/tu_e_bar
        assert input_port>0 , "\nInput port should start from one\n"
        
        # 為了保持微分方程不動，alpha_linear和gamma皆為Amplitude的loss(coupling through ratio)
        self.input_kappa = (2/self.tu_e_bar[input_port-1])**0.5 
        self.alpha_linear = np.real(2/(self.vg*1e-4*self.tu_o_bar)) #Energy absorption (1/cm)
        self.tu_t_bar = (self.tu_e_bar_total_inv+1/self.tu_o_bar)**(-1)
        self.beta_TPA = beta_TPA  #1e-13 (cm/mW)
        self.tau_eff = tau_eff      #1e-9  s,
        self.sigma_FCA = sigma_FCA  #1e-17 cm^2
        self.Atpa = Atpa
        self.Afca = Afca
        self.Akerr = Akerr
        self.n2 = n2
        super().__init__()
        self.lambda_incident = lambda_incident
        
        self.FCA_fit_factor = FCA_fit_factor
        self.TPA_fit_factor = TPA_fit_factor
        self.SPM_fit_factor = SPM_fit_factor
        # normalized dw/dlambda
        self.D_bar = -2*np.pi*c/self.lambda_incident**2 * t0
        self.vg_in_cm = self.vg*1e-4
        # Note. photon energy is in unit. mJ, and normalized by t0

        
        self.Q = np.real( ( (2*self.tu_e_bar_total_inv + 2/self.tu_o_bar  )/(self.f_res_bar*2*np.pi) )**(-1) )
        # self.Q = np.real( ( (1/self.tu_e_bar + 1/self.tu_o_bar  )/(self.f_res_bar*np.pi) )**(-1) )
        
        self.renew()
        

    def renew(self):
        self.photon_energy = h*c/self.lambda_incident*1000/t0  #mJ
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
        A_TPA = self.Atpa
        self.TPA_coeff = self.beta_TPA/(self.round_trip_time*A_TPA)
        self.TPA_coeff_order = np.log10(1e-13/(1e-12*1e-8)) + np.log10(self.TPA_coeff)//1
        self.TPA_coeff = self.TPA_coeff / ( 10**(np.log10(self.TPA_coeff)//1) )
        self.TPA_coeff = self.TPA_coeff*10**( self.TPA_coeff_order)
        self.TPA_coeff = self.TPA_coeff*self.TPA_fit_factor
        # Note. The unit TPA_coeff here is  1/(cm*mJ), not 1/cm

        A_FCA = self.Afca
        self.FCA_coeff = self.beta_TPA*self.tau_eff*self.sigma_FCA/ \
            (2*self.photon_energy*self.round_trip_time**2*A_FCA**2) 
        self.FCA_coeff_order = np.log10(1e-13*1e-9*1e-17/( 1e-12*(1e-12)**2*(1e-4)**4 ) ) + np.log10(self.FCA_coeff)//1
        self.FCA_coeff = self.FCA_coeff / ( 10**(np.log10(self.FCA_coeff)//1) )
        self.FCA_coeff = self.FCA_coeff * 10**(self.FCA_coeff_order)
        self.FCA_coeff = self.FCA_coeff*self.FCA_fit_factor
        # Note. The unit FCA_coeff here is  1/(cm*mJ^2), not 1/cm

        # Self Phase Modulation
        dn_SPM = self.n2/(self.round_trip_time*1e-12*self.Akerr)
        self.df_SPM_coeff = -dn_SPM*self.f_res_bar/self.ng *self.SPM_fit_factor  # 1/(ps*mJ)

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
            m = self.neff(0)*self.L/1.55
        if self.band=="O":
            f_band = c/1.31*t0
            m = self.neff(0)*self.L/1.31
        a = np.ceil(m)*c/self.neff(0)/self.L*t0 - f_band
        b = f_band - np.floor(m)*c/self.neff(0)/self.L*t0
        if a<b:
            m= np.ceil(m)
        if a>b:
            m=np.floor(m)
        self.lambda0_reference = (self.neff(0)*self.L)/m
        
        return m
    
    def find_accurate_res_frequency(self,m):
        f0_v0 = c*t0/self.lambda0_reference
        m_pround = m - self.FSR_shift
        # Assume the specified neff (self.neff) is index of resonant wavelength near 1550 or 1310 nm 
        # solving index at shifted resonance frequency
        A = 1
        B = -(2*self.neff(0)-self.ng)
        C = -(self.ng-self.neff(0))/(f0_v0*self.L)*m_pround*c*t0
        n = np.real((-B+sqrt(B**2-4*A*C))/(2*A))
        assert np.imag((-B+sqrt(B**2-4*A*C))/(2*A))==0.0 , "\nimaginary part of index is supposed to be zero\n"
        return m_pround*c/n/self.L*t0
    
    def CMT(self,f_pround_bar,b_bar,N_bar,delta_T,f_res_bar,alpha_linear,TPA,SPM,T_args,dlambda,Heater):
        da_dt = 1j*2*np.pi*(f_res_bar-f_pround_bar + SPM*abs(b_bar)**2 + \
        \
        (-self.f_res_bar/self.ng)*T_args[0]*delta_T )*b_bar \
        \
        - (self.tu_e_bar_total_inv + \
        \
        self.vg_in_cm*alpha_linear/2 +\
        self.vg_in_cm*TPA/2*abs(b_bar)**2 \
        + self.vg_in_cm*N_bar*1e-5/2  ) *b_bar + \
        \
        self.input_kappa *1 + \
        \
        1j*self.D_bar*dlambda*b_bar + \
        \
        1j*self.D_bar*( self.dlambda_dT*(Heater.T_surround-300) + (self.HE*1e-6)*Heater.P)*b_bar
        return da_dt

    def FC_rate_equation(self,b_bar,N_bar,FCA,tau_eff):

        dN_dt = -t0*N_bar/(tau_eff*1e-9) + FCA*abs(b_bar)**4
        
        return dN_dt
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
    
class alpha_fit():
    a = 0
    b = 0
    c = 0
    """
    alpha is AMPLITUDE absorption coefficient in 1/cm
    """
    @staticmethod
    def func(v, a, b,c):
            return a*v/(abs(v)+b)**0.5 +c
    def renew(self):
        if self.input_mode=="amp":
            self.alpha_data = -2/(self.L*1e-4)*np.log(self.RoundTripLoss)
        if self.input_mode=="energy":
            self.alpha_data = -1/(self.L*1e-4)*np.log(self.RoundTripLoss)
        V = np.array([0,-0.5,-1,-1.5,-2])
        if self.fit_mode=="linear":
            self.m, self.b = np.polyfit(V, self.alpha_data, 1)
        if self.fit_mode == "func":
            self.popt, pcov = curve_fit(self.func, V, self.alpha_data)
    def __init__(self,RoundTripLoss:np.ndarray,L,input,fit_mode = "linear"):
        """input: specify which physic parameter of RoundTripLoss is using, energy or amplitude."""
        """Now, alpha_pdk is Loss of Energy"""
        self.L = L
        self.fit_mode = fit_mode
        self.RoundTripLoss = RoundTripLoss
        self.input_mode = input 
        self.renew()
    def alpha_V(self,V):
        if self.fit_mode=="linear":
            return self.m*(V)+self.b
        if self.fit_mode=="func":
            return self.func(V,self.popt[0],self.popt[1],self.popt[2])
    
class neff_fit(ring):
    c = 0
    b = 0
    a = 0
    @staticmethod
    def func(v, a, b,c):
            return a*v/(abs(v)+b)**0.5 +c
    def renew(self):
        V = np.array([0.5,0,-0.5,-1,-1.5,-2])
        if self.fit_mode=="linear":
            self.m, self.b = np.polyfit(V, self.neff_data, 1)
        if self.fit_mode == "func":
            self.popt, pcov = curve_fit(self.func, V, self.neff_data)
    def __init__(self,neff_data:np.ndarray,fit_mode = "linear"):
        self.fit_mode = fit_mode
        self.neff_data =neff_data
        self.renew()
    def neff_V(self,V):
        if self.fit_mode=="linear":
            return self.m*(V)+self.b
        if self.fit_mode=="func":
            return self.func(V,self.popt[0],self.popt[1],self.popt[2])