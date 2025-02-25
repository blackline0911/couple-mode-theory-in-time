from user import *
from utility import *
from cmt_solver import *
from scipy.integrate import *
import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import math
from scipy.signal import *
from scipy import signal
from random import *

ring_mod = ring(radius=5.0, 
            neff=2.35,
            ng=4.3, 
            gamma = 0.95012, 
            alpha = 0.95246,
            me=37.2)

v = driver(f_drive=100,v_bias=-1,vpp=0,R=53.9)
t = time(N=10,
         lambda_incident=1.5488,
         ring=ring_mod,
         driver=v)
v.create_voltage(time=t)
print(t.dt)

# ploting(t.t_all_segment[0][:],t.t_all_segment[0][:],'time','time test',filename='test_time')
ploting(t.t_total*t0,v.v,'time','time test',filename='test_time')
print(v.v_dict[0.0])
v.varying_Cj()
b,Q,s_minus = solving(ring_mod,v,t,lambda_incident=1.5488,Pin=1e-3)


