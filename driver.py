import numpy as np
from cmath import *
from scipy.signal import *
from scipy import signal
from random import *
from utility import *
from sim import simulation
from time_class import *
import raise_cosine
from scipy.interpolate import CubicSpline


class driver(simulation) :
    Cj = 0
    varying_cj = 0
    Cj_last = 0
    vj = 0
    id='driver'
    method = "small_signal"
    cj_normalizing = 0

    # This voltage_record is used for calculating differetial of voltage
    voltage_record = []

    # Notes : a, b are variables in Junction capacitance formula
    a = 0
    b = 0
    c = 0
    d = 0

    # shift bit is used to shift the time in raise cosine function,
    # when using two NRZ signal to create a raise cosine signal, the first NRZ signal is used as the reference, and the second NRZ signal is shifted by shift_bit bits.
    shift_bit = 100
    def renew(self):
        # Renew driver data earned by calculation, such as w_drive, Cj

        # If you want to use two NRZ signal to combine a new NRZ, don't shift bit.
        # If you want to use two NRZ signal to combine a new PAM4, please shift bit.
        if (self.level == "NRZ") and (self.segment==1):
            self.level_num = 2
            self.shift_bit = 0
        if (self.level == "NRZ") and (self.segment==2):
            self.level_num = 2
            self.shift_bit = 0
        if (self.level == "PAM4") and (self.segment==1):
            self.level_num = 4
            self.shift_bit = 0
        if (self.level == "PAM4") and (self.segment==2):
            self.level_num = 2
            self.shift_bit = 10
        self.w_drive = 2*np.pi*self.f_drive
        

        assert len(self.cj_array[0])>=2, "\n\ncj must have at least two elements\n\n"
        assert np.any( [ self.cj_array[0,0] > self.cj_array[0,1] ] ) , "\n\nYou must input cj as a decreasing array ,since Junction capacitance decrease with voltage\n\n"
        assert np.any( [ self.cj_array[self.segment-1,0] > self.cj_array[self.segment-1,1]  ] ) , "\n\nYou must input cj as a decreasing array ,since Junction capacitance decrease with voltage\n\n"
        self.  b = ( self.cj_array[0,1]**2*(-1)  ) / (self.cj_array[0,1]**2 - self.cj_array[0,0]**2)
        self.  a = self.cj_array[0,0]*(self.  b)**0.5
        if self.segment == 2:
            self.  d = ( self.cj_array[1,1]**2*(-1)  ) / (self.cj_array[1,1]**2 - self.cj_array[1,0]**2)
            self.  c = self.cj_array[1,0]*(self.  d)**0.5
            self.cj_normalizing2 = self.Cj_V(self.v_bias,self.  c,self.  d)
            assert ( not (self.  c==0) ) and ( not (self.  d==0) ), "\n\nYou must calculate it first\n\n"
            self.C = self.  c**2/(self.cj_normalizing2**2)
            self.D = self.cj_normalizing2**2 / (2*self.  c**2)

        self.cj_normalizing1 = self.Cj_V(self.v_bias,self.  a,self.  b)
        assert ( not (self.  a==0) ) and ( not (self.  b==0) ), "\n\nYou must calculate it first\n\n"
        self.A = self.  a**2/(self.cj_normalizing1**2)
        self.B = self.cj_normalizing1**2 / (2*self.  a**2)
        
    def __init__(self,
                 f_drive,
                 v_bias,
                 vpp,
                 Rs:np.ndarray,
                 Cox:np.ndarray,
                 Rsi:np.ndarray,
                 Cp:np.ndarray,
                 cj:np.ndarray,
                 level,
                 square_wave = 0,
                 sine_wave =0,
                 raise_cosine =0,
                 PRBS = 1,
                 beta = 0.6,
                 Z0 = 50,
                 method = "small_signal",
                 segment = 1
                 ):
        """
        input argments:
        f_drive: driving voltage frequency (GHz)
        v_bias: bias voltage
        vpp: peak to peak voltage
        R:series resistance of PN junction 
        cj : a decreasing array depicting Cj varying with voltage. Note first element is 0V, -1V ,and so on.
        cj must have at least two elements
        """
        self.f_drive = f_drive*1e9
        self.w_drive = 2*np.pi*self.f_drive
        self.v_bias = v_bias
        self.vpp=vpp
        self.Rs = Rs
        self.Cox = Cox
        self.Rsi = Rsi
        self.Cp = Cp
        self.segment = segment
        assert len(Rs)==segment, "\nRs should be a 1D array with length of segment\n"
        assert len(Cox)==segment, "\nCox should be a 1D array with length of segment\n"
        assert len(Rsi)==segment, "\nRsi should be a 1D array with length of segment\n"
        assert len(Cp)==segment, "\nCp should be a 1D array with length of segment\n"
        self.Z0 = Z0

        self.square_wave = square_wave
        self.sine_wave = sine_wave
        self.raise_cosine = raise_cosine
        self.PRBS = PRBS
        self.method = method
        self.beta = beta

        self.cj_array = cj
        self.cj_array = np.array(self.cj_array)
        self.cj_array = self.cj_array.reshape((self.segment,2))
        assert np.size(self.cj_array)==self.segment*2, "\n\ncj must have at least two elements\n\n"
        self.level = level
        self.renew()
        self.cj_normalizing = self.Cj_V(self.v_bias,self.  a,self.  b)
        
        return 
    def create_voltage(self, 
                       time,
                       ):
        """
        Do not call this function Externally\n

        This function creates a voltage array according to the time array refered to the time object.
        Hence the voltage value only exists in the specified time.
        Since I have no idea what time will solve_ivp solver need , so this function is only working for analysis and plotting
        """
        self.time = time
        if time.mode == "voltage_drive":
            self.Cj = self.Cj_V(self.v_bias,self.a,self.b)
            # Generating PRBS 
            if self.PRBS:
                self.prbs =(np.zeros(self.time.N))
                if (self.level == 'NRZ' and self.segment == 2) or \
                    (self.level == 'NRZ' and self.segment == 1) or \
                    (self.level == 'PAM4' and self.segment == 2)    :
                    for i in range(self.time.N):
                        self.prbs[i] = (randint(0,1))
                if self.level == "PAM4" and self.segment == 1:
                    for i in range(self.time.N):
                        self.prbs[i] = (randint(0,3))
                self.prbs[int(randint(0,time.N-1))] = (1)
            else:
                # Periodic Bit Sequence (this part only used when raise cosine is true while PRBS is false)
                # create a dummy bit sequence : 0,1,0,1,0,1...
                self.bit_sequence  = np.array([])
                T_normalized = self.time.T_normalized
                for u in range(int(self.time.t_max/T_normalized)):
                    self.bit_sequence = (np.append( self.bit_sequence , float(u%2) ))

            
            if self.square_wave:
                self.v = self.vpp/2*signal.square(self.w_drive*self.time.t_total*t0,duty=0.5)+self.v_bias
                assert not (self.sine_wave or self.raise_cosine) , "Only one kind of signal should apply "
                self.v = a + self.v_bias - self.vpp/2
                self.cubic_interp_voltage = CubicSpline(self.time.t_total,self.v)
            if self.sine_wave:
                self.v = self.vpp/2*np.exp(1j*self.w_drive*self.time.t_total*t0)+self.v_bias
                assert not (self.square_wave or self.raise_cosine) , "Only one kind of signal should apply "
                self.v = a + self.v_bias - self.vpp/2
                self.cubic_interp_voltage = CubicSpline(self.time.t_total,self.v)
            if self.raise_cosine:
                assert not (self.square_wave or self.sine_wave) , "Only one kind of signal should apply "
                if self.segment == 1:
                    if self.PRBS:
                        a = raise_cosine.create_rcos_signal(self.prbs,time.t_total,time.T_normalized,time.N,self,self.beta)
                    else:
                        a = raise_cosine.create_rcos_signal(self.bit_sequence,time.t_total,time.T_normalized,time.N,self,self.beta)
                    self.v = a + self.v_bias - self.vpp/2
                    self.cubic_interp_voltage = CubicSpline(self.time.t_total,self.v)
                if self.segment == 2:
                    if self.PRBS:
                        # vpp_temp = self.vpp
                        # self.vpp = self.vpp*2/3
                        a = raise_cosine.create_rcos_signal(self.prbs,time.t_total,time.T_normalized,time.N,self,self.beta)
                    else:
                        a = raise_cosine.create_rcos_signal(self.bit_sequence,time.t_total,time.T_normalized,time.N,self,self.beta)
                    a = a - self.vpp/2
                    length_UI = len(self.time.t_all_segment[0])
                    shift_NRZ_wobias = a [int(self.shift_bit*length_UI):]
                    if not (self.shift_bit==0):
                        a = a[:-int(self.shift_bit*length_UI)]
                        t_slice = self.time.t_total[:-int(self.shift_bit*length_UI)]
                    else:
                        t_slice = self.time.t_total
                    a_wobias = a 
                    a_wbias = a +self.v_bias
                    shift_NRZ_wbias = shift_NRZ_wobias + self.v_bias
                    self.v_LSB = a_wbias
                    self.v_MSB = shift_NRZ_wbias
                    self.v = a_wobias + shift_NRZ_wobias + self.v_bias
                    self.cubic_interp_voltage = CubicSpline(t_slice,self.v_LSB)
                    self.cubic_interp_voltage_shift = CubicSpline(t_slice,self.v_MSB )
                    # self.vpp = vpp_temp
            
        if self.method == 'small_signal':
            self.Cj = self.Cj_V(self.v_bias,self.a,self.b)
        
        return
    step =0 
    passed = False
    counter = 0
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
                v = self.vpp/2*signal.square(self.w_drive*t*t0,duty=0.5)+self.v_bias
            if self.sine_wave:
                assert not (self.square_wave or self.raise_cosine) , "Only one kind of signal should apply "
                v = self.vpp/2*np.exp(1j*self.w_drive*t*t0)+self.v_bias
            if self.raise_cosine:
                assert not (self.square_wave or self.sine_wave) , "Only one kind of signal should apply "
                v = self.cubic_interp_voltage(t)
            return v
        if self.time.mode == "scan_frequency":
            return self.v_bias
    def refering_v_shifted(self,t):
        v = self.cubic_interp_voltage_shift(t)
        return v
    def refering_dv_dt(self,v0,t):
        # v0：voltage point at time t
        if (t<2*self.time.dt/t0):
            v1 = self.cubic_interp_voltage(t+self.time.dt/t0)
            v2 = self.cubic_interp_voltage(t+2*self.time.dt/t0)
            dv_dt = (-3*v0 + 4*v1 - v2)/2/(self.time.dt/t0)
        if (t>=2*self.time.dt/t0):
            v1 = self.cubic_interp_voltage(t-self.time.dt/t0)
            v2 = self.cubic_interp_voltage(t-2*self.time.dt/t0)
            dv_dt = (3*v0 - 4*v1 + v2)/2/(self.time.dt/t0)
        return dv_dt
    
    def refering_dv_dt_shifted(self,v0,t):
        # v0：voltage point at time t
        if (t<2*self.time.dt/t0):
            v1 = self.cubic_interp_voltage_shift(t+self.time.dt/t0)
            v2 = self.cubic_interp_voltage_shift(t+2*self.time.dt/t0)
            dv_dt = (-3*v0 + 4*v1 - v2)/2/(self.time.dt/t0)
        if (t>=2*self.time.dt/t0):
            v1 = self.cubic_interp_voltage_shift(t-self.time.dt/t0)
            v2 = self.cubic_interp_voltage_shift(t-2*self.time.dt/t0)
            dv_dt = (3*v0 - 4*v1 + v2)/2/(self.time.dt/t0)
        return dv_dt

        

    def refering_Cj(self,voltage):
        """
        input
        t: given a arbitary voltage value, this function can return the corresponding capacitance value
        """
        if self.varying_cj:
            return self.Cj_V(voltage)
        else:
            return float(self.Cj)
        
    def Cj_V(self,vol,a,b):
        assert ( not (a==0) ) and ( not (b==0) ), "\n\nYou must calculate it first\n\n"
        return a/( (b-vol)**0.5 )
    
    def Q_V(self,vol,a,b):
        assert ( not (a==0) ) and ( not (b==0) ), "\n\nYou must calculate it first\n\n"
        return a/( (b-vol)**0.5 )*vol
    
    def V_Q(self,Q_bar,A,B,a,b):
        assert ( not (a==0) ) and ( not (b==0) ), "\n\nYou must calculate it first\n\n"
        return  (-Q_bar**2*B + \
                B*Q_bar*(Q_bar**2 + 4*b*A)**0.5 )
    
    def varying_Cj(self):
        """
        call this function when you want to analyze large signal 
        """
        self.varying_cj=1
        self.Cj = self.Cj_V(self.v)
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
                                 int(self.prbs[int(shift/T)])/(self.level_num-1)*np.pi/4*sinc((1/2/beta)),
                                 int(self.prbs[int(shift/T)])/(self.level_num-1) * sinc((t-shift)/T) * np.cos(np.pi*beta*(t-shift)/T) / (1-(2*beta*(t-shift)/T)**2) )
                return ans
            else:
                ans  = np.where( ( (t-shift)==T/2/beta) | ( (t-shift)==( -T/2/beta) ),
                                 self.bit_sequence[int(shift/T)]*np.pi/4*sinc((1/2/beta)),
                                 self.bit_sequence[int(shift/T)]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2) )
                return ans
        else:
            """When input is a float number"""
            if ((t-shift)==T/2/beta or (t-shift)==( -T/2/beta)):
                if self.PRBS:
                    return self.prbs[int(shift/T)]*np.pi/4*sinc(1/2/beta)
                else:
                    return self.bit_sequence[int(shift/T)]*np.pi/4*sinc(1/2/beta)
            else:
                if self.PRBS:
                    return self.prbs[int(shift/T)]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)
                else:
                    return self. bit_sequence[int(shift/T)]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)

    def TML(self,v_neg, vj, voltage, dvpos_dt, i2, Z0,Rs,Cj_bar,Cox_bar,Rsi,Cp_bar):
        dvneg_dt = ( 1/Cp_bar*1/Z0*(voltage - v_neg) \
                -dvpos_dt\
                -1/Cp_bar/Rs*(voltage + v_neg - vj)\
                -1/Cp_bar*i2)
    
        dvj_dt = 1/(Rs*Cj_bar)*(voltage + v_neg - vj)

        di2_dt = (1/Rsi*dvpos_dt + \
                1/Rsi*dvneg_dt\
                - 1/Rsi/Cox_bar*i2)
    
        return dvneg_dt, dvj_dt, di2_dt