from scipy.integrate import *
import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import math
from scipy.signal import *
from scipy import signal
from random import *
from utility import *

class time():
    def __init__(self,N:int,
                 lambda_incident:float,
                 ring,
                 driver):
        self.N=N
        w_pround = 2*np.pi*c/(lambda_incident)*t0
        df_max =  (c/ring.lambda0**2)*abs( ring.me*1e-12/1e-6 )*( driver.vpp/2 + abs(driver.v_bias))
        if df_max < driver.f_drive*1e9:
            self.dt = 1/(driver.f_drive*1e9)/10
        self.dt = (1/df_max/100)
        inte = (math.log10(self.dt))//1
        self.dt = 10**(inte-1)

        self.t_max = math.ceil( 1/(driver.f_drive*1e9)*N/t0 )
        self.T_normalized = 1/(driver.f_drive*1e9)/t0
        
        self.t_total = np.array([])
        self.t_all_segment = np.zeros((self.N,1)).tolist()

        # recording length in each time segment
        self.number_record = np.array([])
        for r in range(self.N):
            # num = math.ceil( ( (r+1)*self.T_normalized-self.dt/t0- (0+(r)*self.T_normalized ) ) / ( self.dt/t0 ))
            # t_segment = np.linspace(0+(r)*self.T_normalized , (r+1)*self.T_normalized-self.dt/t0, num )
            t_segment = np.arange( 0+r*self.T_normalized ,  (r+1)*self.T_normalized, self.dt/t0)
            self.number_record= np.append(self.number_record,len(t_segment))
            self.t_total=np.append(self.t_total,t_segment)
            self.t_all_segment[r] = t_segment
            # self.t_all_segment[r][0] = t_segment[0]
            # self.t_all_segment[r].extend(t_segment[1:])
        
        return
    
    
    # def create_array(self):
    #     return
    
class ring():
    def __init__(self,radius:float, 
                 neff:float,
                 ng:float,
                 gamma:float,
                 alpha:float,
                 me:float):
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
        self.w_res_bar = 102*np.pi*c/(self.neff*self.L)*t0
        self.lambda0 = 2*np.pi*c*t0/self.w_res_bar


class driver():
    Cj = 0
    def __init__(self,
                 f_drive,
                 v_bias,
                 vpp,
                 R:float,
                 ):
        """
        input argments:
        f_drive: driving voltage frequency (GHz)
        v_bias: bias voltage
        vpp: peak to peak voltage
        R:series resistance of PN junction 
        """
        self.f_drive = f_drive
        self.w_drive = 2*np.pi*f_drive*1e9
        self.v_bias = v_bias
        self.vpp=vpp
        self.R = R
        return 
    def create_voltage(self, 
                       time,
                       square_wave = 0,
                       sine_wave =1,
                       raise_cosine =0,
                       PRBS = 0,
                       ):
        self.time = time
        self.square_wave = square_wave
        self.sine_wave = sine_wave
        if square_wave:
            self.v = self.vpp/2*signal.square(self.w_drive*time.t_total*t0,duty=0.5)+self.v_bias
        if sine_wave:
            self.v = self.vpp/2*np.exp(1j*self.w_drive*time.t_total*t0)+self.v_bias
        # if raise_cosine:
        #     T_period_normalized = self.time.T_normalized
        #     a=0
        #     if PRBS:
        #         for i in range(self.time.N):
        #             a +=  (vpp*rcos(t, beta, T_period,(i+1)*T_period_normalized))
        #     else:
        #         for i in range(len(bit_sequence)):
        #             a +=  (vpp*rcos(t, beta, T_period,(i+1)*T_period_normalized))
        #     return a+v_bias-vpp/2  
    
        self.v_dict = dict(zip(time.t_total,self.v))
        self.Cj = self.Cj_V(self.v_bias)

        return

    def refering_v(self,t):
        if self.square_wave:
            return self.vpp/2*signal.square(self.w_drive*t*t0,duty=0.5)+self.v_bias
        if self.sine_wave:
            return self.vpp/2*np.exp(1j*self.w_drive*t*t0)+self.v_bias

    def refering_Cj(self,voltage):
        return self.Cj_V(voltage)
        
    def Cj_V(self,vol):
        return 3.7675e-14/( (2.5485-vol)**0.5 )
    
    def varying_Cj(self):
        """
        call this function when you want to analyze large signal 
        """
        self.Cj = self.Cj_V(self.v)
        self.Cj_dict = dict(zip(self.v,self.Cj))
        return 


        