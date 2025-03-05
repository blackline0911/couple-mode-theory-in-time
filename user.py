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
        if df_max < driver.f_drive:
            self.dt = 1/(driver.f_drive)/10
        self.dt = (1/df_max/100)
        inte = (math.log10(self.dt))//1
        self.dt = 10**(inte-1)

        self.t_max = math.ceil( 1/(driver.f_drive)*N/t0 )
        self.T_normalized = 1/(driver.f_drive)/t0
        
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


class driver(time) :
    Cj = 0
    def __init__(self,
                 f_drive,
                 v_bias,
                 vpp,
                 R:float,
                 square_wave = 0,
                 sine_wave =1,
                 raise_cosine =0,
                 PRBS = 0,
                 ):
        """
        input argments:
        f_drive: driving voltage frequency (GHz)
        v_bias: bias voltage
        vpp: peak to peak voltage
        R:series resistance of PN junction 
        """
        self.f_drive = f_drive*1e9
        self.w_drive = 2*np.pi*f_drive
        self.v_bias = v_bias
        self.vpp=vpp
        self.R = R

        self.square_wave = square_wave
        self.sine_wave = sine_wave
        self.raise_cosine = raise_cosine
        self.PRBS = PRBS
        return 
    def create_voltage(self, 
                       time,
                       ):
        """
        This function creates a voltage array according to the time array refered to the time object.mro
        Hence the voltage value only exists in the specified time.mro
        Since I have no idea what time will solve_ivp solver need , so this function is only working for analysis and plotting
        """
        self.time = time
        # Generating PRBS 
        if self.PRBS:
            self.prbs =np.zeros(self.time.N)
            for i in range(self.time.N):
                self.prbs[i] = randint(0,1)
        else:
            # Periodic Bit Sequence (this part only used when raise cosine is true while PRBS is false)
            # create a dummy bit sequence : 0,1,0,1,0,1...
            self.bit_sequence  = np.array([])
            T_normalized = self.time.T_normalized
            for u in range(int(self.time.t_max/T_normalized)):
                self.bit_sequence = np.append( self.bit_sequence , u%2 )

        
        if self.square_wave:
            self.v = self.vpp/2*signal.square(self.w_drive*time.t_total*t0,duty=0.5)+self.v_bias
        if self.sine_wave:
            self.v = self.vpp/2*np.exp(1j*self.w_drive*time.t_total*t0)+self.v_bias
        if self.raise_cosine:
            T_period_normalized = self.time.T_normalized
            a=0
            if self.PRBS:
                for i in range(self.time.N):
                    a +=  (self.vpp*self.rcos(self.time.t_total,(i+1)*T_period_normalized))
            else:
                for i in range(len(self.bit_sequence)):
                    a +=  (self.vpp*self.rcos(self.time.t_total, (i+1)*T_period_normalized))
            self.v = a+self.v_bias-self.vpp/2  
    
        # self.v_dict = dict(zip(time.t_total,self.v))
        self.Cj = self.Cj_V(self.v_bias)

        return

    def refering_v(self,t):
        """
        input
        t: given a arbitary time value, this function can return the corresponding voltage value
        """
        if self.square_wave:
            """Do not use square wave is better, 
            since the discontinuous in the transition can result in artifitial peak in large signal analysis"""
            return self.vpp/2*signal.square(self.w_drive*t*t0,duty=0.5)+self.v_bias
        if self.sine_wave:
            return self.vpp/2*np.exp(1j*self.w_drive*t*t0)+self.v_bias
        if self.raise_cosine:
            T_period_normalized = self.time.T_normalized
            a=0
            for i in range(self.time.N):
                a +=  (self.vpp*self.rcos(t,shift=(i+1)*T_period_normalized))
            return a+self.v_bias-self.vpp/2  
            

    def refering_Cj(self,voltage):
        """
        input
        t: given a arbitary voltage value, this function can return the corresponding capacitance value
        """
        return self.Cj_V(voltage)
        
    def Cj_V(self,vol):
        return 3.7675e-14/( (2.5485-vol)**0.5 )
    
    def varying_Cj(self):
        """
        call this function when you want to analyze large signal 
        """
        self.Cj = self.Cj_V(self.v)
        # self.Cj_dict = dict(zip(self.v,self.Cj))
        return 

    def sinc(self,t):
        if isinstance(t, np.ndarray):
            return np.where(t==0,1,np.sin(np.pi*t)/(np.pi*t))
        else:
            if t==0:
                return 1
            else:
                return np.sin(np.pi*t)/(np.pi*t)

    def rcos(self, 
             t, 
             shift, 
             beta=1):
        """
        input:
            shift:time shift of raise cosine function, must be integar of normalized T 
            T:period of a bit 
            t: giving time second(normalized)
        """
        T = self.time.T_normalized
        assert beta !=0, "The raise cosine function is degraded to sinc function !"

        if isinstance(t, np.ndarray):
            """When input is a numpy array"""
            if self.PRBS:
                ans  = np.where( ( (t-shift)==T/2/beta) | ( (t-shift)==( -T/2/beta) ),
                                 self.prbs[int(shift/T)-1]*np.pi/4*self.sinc((1/2/beta)),
                                 self.prbs[int(shift/T)-1] * self.sinc((t-shift)/T) * np.cos(np.pi*beta*(t-shift)/T) / (1-(2*beta*(t-shift)/T)**2) )
                return ans
            else:
                ans  = np.where( ( (t-shift)==T/2/beta) | ( (t-shift)==( -T/2/beta) ),
                                 self.bit_sequence[int(shift/T)-1]*np.pi/4*self.sinc((1/2/beta)),
                                 self.bit_sequence[int(shift/T)-1]*self.sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2) )
                return ans
        else:
            """When input is a float number"""
            if ((t-shift)==T/2/beta or (t-shift)==( -T/2/beta)):
                if self.PRBS:
                    return self.prbs[int(shift/T)-1]*np.pi/4*self.sinc(1/2/beta)
                else:
                    return self.bit_sequence[int(shift/T)-1]*np.pi/4*self.sinc(1/2/beta)
            else:
                if self.PRBS:
                    return self.prbs[int(shift/T)-1]*self.sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)
                else:
                    return self. bit_sequence[int(shift/T)-1]*self.sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)



        