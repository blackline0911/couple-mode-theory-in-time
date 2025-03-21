import numpy as np
from cmath import *
from scipy.signal import *
from scipy import signal
from random import *
from utility import *
from sim import simulation
from time_class import *


class driver(simulation) :
    Cj = 0
    id='driver'
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
        self.w_drive = 2*np.pi*self.f_drive
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
        This function creates a voltage array according to the time array refered to the time object.
        Hence the voltage value only exists in the specified time.
        Since I have no idea what time will solve_ivp solver need , so this function is only working for analysis and plotting
        """
        self.time = time
        if time.mode == "voltage_drive":
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
                self.v = self.vpp/2*signal.square(self.w_drive*self.time.t_total*t0,duty=0.5)+self.v_bias
                assert not (self.sine_wave or self.raise_cosine) , "Only one kind of signal should apply "
            if self.sine_wave:
                self.v = self.vpp/2*np.exp(1j*self.w_drive*self.time.t_total*t0)+self.v_bias
                assert not (self.square_wave or self.raise_cosine) , "Only one kind of signal should apply "
            if self.raise_cosine:
                assert not (self.square_wave or self.sine_wave) , "Only one kind of signal should apply "
                T_period_normalized = self.time.T_normalized
                a=0
                if self.PRBS:
                    for i in range(self.time.N):
                        a +=  (self.vpp*self.rcos(self.time.t_total,(i+1)*T_period_normalized))
                else:
                    for i in range(len(self.bit_sequence)):
                        a +=  (self.vpp*self.rcos(self.time.t_total, (i+1)*T_period_normalized))
                self.v = a+self.v_bias-self.vpp/2  
        if time.mode == "scan_frequency":
            self.v = self.v_bias*np.ones(len(time.t_total))
            # self.v_dict = dict(zip(time.t_total,self.v))
            self.Cj = self.Cj_V(self.v_bias)

            return

    def refering_v(self,t):
        if self.time.mode == "voltage_drive":
            """
            input
            t: given a arbitary time value, this function can return the corresponding voltage value
            """
            if self.square_wave:
                """Do not use square wave is better, 
                since the discontinuous in the transition can result in artifitial peak in large signal analysis"""
                assert not (self.sine_wave or self.raise_cosine) , "Only one kind of signal should apply "
                return self.vpp/2*signal.square(self.w_drive*t*t0,duty=0.5)+self.v_bias
            if self.sine_wave:
                assert not (self.square_wave or self.raise_cosine) , "Only one kind of signal should apply "
                return self.vpp/2*np.exp(1j*self.w_drive*t*t0)+self.v_bias
            if self.raise_cosine:
                assert not (self.square_wave or self.sine_wave) , "Only one kind of signal should apply "
                T_period_normalized = self.time.T_normalized
                a=0
                for i in range(self.time.N):
                    a +=  (self.vpp*self.rcos(t,shift=(i+1)*T_period_normalized))
                return a+self.v_bias-self.vpp/2  
        if self.time.mode == "scan_frequency":
            return self.v_bias
            

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
                                 self.prbs[int(shift/T)-1]*np.pi/4*sinc((1/2/beta)),
                                 self.prbs[int(shift/T)-1] * sinc((t-shift)/T) * np.cos(np.pi*beta*(t-shift)/T) / (1-(2*beta*(t-shift)/T)**2) )
                return ans
            else:
                ans  = np.where( ( (t-shift)==T/2/beta) | ( (t-shift)==( -T/2/beta) ),
                                 self.bit_sequence[int(shift/T)-1]*np.pi/4*sinc((1/2/beta)),
                                 self.bit_sequence[int(shift/T)-1]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2) )
                return ans
        else:
            """When input is a float number"""
            if ((t-shift)==T/2/beta or (t-shift)==( -T/2/beta)):
                if self.PRBS:
                    return self.prbs[int(shift/T)-1]*np.pi/4*sinc(1/2/beta)
                else:
                    return self.bit_sequence[int(shift/T)-1]*np.pi/4*sinc(1/2/beta)
            else:
                if self.PRBS:
                    return self.prbs[int(shift/T)-1]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)
                else:
                    return self. bit_sequence[int(shift/T)-1]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)



        