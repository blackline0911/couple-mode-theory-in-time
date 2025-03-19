from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *

wl_in = 1.5488
Pin = 100  #mW
# Pin = 1  #mW
sim = simulation(wl_in)
wl_min = 1.55
wl_max = 1.548
ring_mod = ring(L=2*np.pi*5, 
            neff=2.51464,
            ng=4.3, 
            gamma = 0.95012, 
            alpha = 0.95246,
            me=37.2,
            cross_section=0.5*0.22)

v = driver(f_drive=100,
           v_bias=-1,
           vpp=0,
           R=53.9)
t = time(mode = "scan_frequency")
# t = time(mode = "voltage_drive")
ring_mod.scan_frequency(wl_min ,wl_max,wl_in,t)

t.main(ring_mod,t_max=2000,resolution=1)

# t.main(ring_mod,N=10,driver=v)
v.create_voltage(time=t)
v.varying_Cj()

wl_scan =  c/ring_mod.w_res(t.t_total)*t0

b,Q,s_minus = solving(ring_mod,v,t,wl_in,Pin=Pin)
b0 = np.real(sqrt(t0)*sqrt(Pin))
print("b0 = ",b0)
T = Transfer_function(ring_mod,t)
wl,data = T.mapping(abs(s_minus)**2/Pin)

os.chdir("./data/")
ploting(t.t_total,ring_mod.w_res(t.t_total),x_label='time',title='f res',filename='f_res')
ploting(t.t_total,abs(b)**2,x_label='time (ps)',title='b (mJ)',filename='b_test')
ploting(t.t_total,wl_scan,x_label='time',title='lambda res (um)',filename='lambda_res')
ploting(t.t_total,Q,x_label='time (ps)',title='Q',filename='Q_test')

ploting(t.t_total,abs(s_minus)**2,x_label='time (ps)',title=r"$|s_-^2|$",filename='s_minus_power')
ploting(t.t_total,180/np.pi*np.angle(s_minus),x_label='time (ps)',title='s_minus phase',filename='s_minus phase')
np.savetxt('b.txt',b/b0,fmt="%.8f", delimiter="\n")

Pin = 1
b,Q,s_minus = solving(ring_mod,v,t,wl_in,Pin=Pin)
b0 = np.real(sqrt(t0)*sqrt(Pin))
T = Transfer_function(ring_mod,t)
wl,data_v1 = T.mapping(abs(s_minus)**2/Pin)
ploting(wl*1000,data,data_v1,x_label='wavelength (nm)',title='Transfer function',filename='T')
# np.savetxt('T.txt',data_v1,fmt="%.8f", delimiter="\n")
# np.savetxt('T_NL_0.txt',data,fmt="%.8f", delimiter="\n")
