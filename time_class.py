
import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import math
from scipy.signal import *
from scipy import signal
from random import *
from utility import *
from sim import simulation
from ring import *
from driver import *


class time():
    
    dt = 0
    t_total = np.array([])
    t_max = 0
    buffer = 0
    def __init__(self,mode):
        self.mode = mode
        modes = {"voltage_drive", "scan_frequency"}
        assert mode in modes, f"Please choose a simulation mode from {modes}"

    def main(self, ring, driver=None,N=0,t_max = 0,buffer=0):
        # 根據 mode 建立對應的子類別實例
        if self.mode == "voltage_drive":
            if N==0:
                assert False ,"please specify N"
            self.N=N
            sim_time = VoltageDriveTime(N)
            self.dt = sim_time.set_dt(ring, driver)
            self.t_total, self.t_all_segment, self.t_max= sim_time.create_time_array(N)
            self.T_normalized = driver.f_drive*t0
        elif self.mode == "scan_frequency":
            if buffer==0 or t_max==0:
                assert False ,"please specify t_max and buffer"
            self.buffer = buffer
            self.t_max = t_max
            sim_time = ScanFrequencyTime(t_max,buffer)
            self.dt = sim_time.set_dt(ring)
            self.t_total = sim_time.create_time_array()
        else:
            raise ValueError("Unknown mode!")

    

class VoltageDriveTime():
    def __init__(self,N:int):
        self.N =N
    """在 voltage_drive 模式底下的子類別"""
    def set_dt(self, ring, driver):
        self.ring = ring
        self.driver = driver
        df_max =  (c/ring.lambda0**2)*abs( ring.me*1e-12/1e-6 )*( driver.vpp/2 + abs(driver.v_bias))
        if df_max < driver.f_drive:
            self.dt = 1/(driver.f_drive)/10
        dt = (1/df_max/100)
        inte = (math.log10(dt))//1
        dt = 10**(inte-2)
        return dt
    
    def create_time_array(self, N):
        self.N = N
        t_total = np.array([])
        t_all_segment = np.zeros((self.N,1)).tolist()
        self.t_max = math.ceil( 1/(self.driver.f_drive)*N/t0 )
        self.T_normalized = 1/(self.driver.f_drive)/t0

        # recording length in each time segment
        self.number_record = np.array([])
        for r in range(self.N):
            # num = math.ceil( ( (r+1)*self.T_normalized-self.dt/t0- (0+(r)*self.T_normalized ) ) / ( self.dt/t0 ))
            # t_segment = np.linspace(0+(r)*self.T_normalized , (r+1)*self.T_normalized-self.dt/t0, num )
            t_segment = np.arange( 0+r*self.T_normalized ,  (r+1)*self.T_normalized, self.dt/t0)
            self.number_record= np.append(self.number_record,len(t_segment))
            t_total=np.append(t_total,t_segment)
            t_all_segment[r] = t_segment
            # self.t_all_segment[r][0] = t_segment[0]
            # self.t_all_segment[r].extend(t_segment[1:])
        
        return t_total, t_all_segment, self.t_max

class ScanFrequencyTime():
    dt = 0
    def __init__(self,t_max:int,buffer):
        self.t_max =t_max
        self.buffer = buffer
        # buffer stands for the time of Transient mode
    """在 scan_frequency 模式底下的子類別"""
    def set_dt(self, ring):
        self.ring = ring
        if abs(ring.f_min_bar-ring.f_res_bar)>abs(ring.f_max_bar-ring.f_res_bar):
            df_max = abs(ring.f_min_bar-ring.f_res_bar)/t0
        else:
            df_max = abs(ring.f_max_bar-ring.f_res_bar)/t0
        dt = (1/df_max/100)
        inte = (math.log10(dt))//1
        dt = 10**(inte-1)
        self.dt = dt
        return dt

    def create_time_array(self):
        """
        In frequency scanning mode, t_max must be integer (ps)
        """
        index_num = int(self.t_max/(self.dt/t0))
        buffer_index_num = int(self.buffer/(self.dt/t0))
        t_total = np.zeros(index_num)
        b = np.zeros(buffer_index_num)
        self.t_max = self.t_max
        assert (not self.dt==0), "dt is zero"
        for j in range(buffer_index_num):
            b[j]=j*(self.dt/t0)
        for i in range(index_num):
            t_total[i] = i*(self.dt/t0)+self.buffer
        t_total = np.append(b,t_total)
        return t_total

