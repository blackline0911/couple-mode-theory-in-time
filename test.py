from utility import *
import numpy as np
from cmath import *
import os
from ring import *
from driver import driver
from time_class import time
from cmt_solver import *

wl_in = 1.5488
Pin = 1
sim = simulation(wl_in)
ring_mod = ring(radius=5.0, 
            neff=2.51464,
            ng=4.3, 
            gamma = 0.95012, 
            alpha = 0.95246,
            me=37.2)

v = driver(f_drive=100,
           v_bias=-1,
           vpp=0,
           R=53.9)
# t = time(mode = "scan_frequency")
t = time(mode = "voltage_drive")

ring_mod.scan_frequency(1.5485,1.5498,1.5488,t)

# t.main(ring_mod,t_max=300,buffer=80)
t.main(ring_mod,N=10,driver=v)
v.create_voltage(time=t)

v.varying_Cj()

wl_scan =  c/ring_mod.w_res(t.t_total)*t0

ploting(t.t_total,ring_mod.w_res(t.t_total),'time','f res',filename='f_res')
ploting(t.t_total,wl_scan,'time','lambda res (um)',filename='lambda_res')
b,Q,s_minus = solving(ring_mod,v,t,wl_in,Pin=Pin)


# os.chdir("./data/")
# ploting(t.t_total,np.real(b),'time (ps)','b',filename='b_test')
# os.chdir("../")

# os.chdir("./data/")
# t = time(N=10,
#          lambda_incident=1.548,
#          ring=ring_mod,
#          driver=v)
# b,Q,s_minus = solving(ring_mod,v,t,Pin=Pin)
# ploting(t.t_total,np.real(b),'time (ps)','b',filename='b_test_v2')
# os.chdir("../")

# os.chdir("./data/")
# t = time(N=10,
#          lambda_incident=1.5485,
#          ring=ring_mod,
#          driver=v)
# b,Q,s_minus = solving(ring_mod,v,t,Pin=Pin)
# ploting(t.t_total,np.real(b),'time (ps)','b',filename='b_test_v3')
# os.chdir("../")

ploting(t.t_total,Q,'time (ps)','Q',filename='Q_test')
ploting(t.t_total,abs(s_minus)**2,'time (ps)',r"$|s_-^2|$",filename='s_minus_power')
ploting(t.t_total,180/np.pi*np.angle(s_minus),'time (ps)','s_minus phase',filename='s_minus phase')

T = Transfer_function(ring_mod,t,sim)
wl,data = T.mapping(abs(s_minus)**2)

ploting(wl*1000, data, 'wavelength (nm)', 'Transfer function of ring', filename="T")

b0 = sqrt(t0)*sqrt(Pin)
np.savetxt("b_bar.txt", b/b0, fmt="%.8f", delimiter="\n")
np.savetxt("Q_bar.txt", Q/v.Cj, fmt="%.8f", delimiter="\n")
np.savetxt("s_minus_power.txt", abs(s_minus)**2, fmt="%.8f", delimiter="\n")
np.savetxt("s_minus_phase.txt", 180/np.pi*np.angle(s_minus), fmt="%.8f", delimiter="\n")

