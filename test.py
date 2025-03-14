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
wl_min = 1.548
wl_max = 1.5498
ring_mod = ring(L=2*np.pi*5, 
            neff=2.51464,
            ng=4.3, 
            gamma = 0.95012, 
            alpha = 0.95246,
            me=37.2)
# print("FSR = ",1000*ring_mod.lambda0**2/(ring_mod.ng*ring_mod.L))
v = driver(f_drive=100,
           v_bias=0,
           vpp=0,
           R=53.9)
t = time(mode = "scan_frequency")
# t = time(mode = "voltage_drive")

ring_mod.scan_frequency(wl_min ,wl_max,wl_in,t)

t.main(ring_mod,t_max=10000,resolution=1)

print( (ring_mod.f_start_bar-ring_mod.f_res_bar)*t.t_max/(ring_mod.f_start_bar-ring_mod.f_end_bar)+t.buffer)

# t.main(ring_mod,N=10,driver=v)
v.create_voltage(time=t)
v.varying_Cj()

wl_scan =  c/ring_mod.w_res(t.t_total)*t0


b,Q,s_minus = solving(ring_mod,v,t,wl_in,Pin=Pin)
b0 = sqrt(t0)*sqrt(Pin)

T = Transfer_function(ring_mod,t)
wl,data = T.mapping(abs(s_minus)**2/Pin)

os.chdir("./data/")
ploting(t.t_total,ring_mod.w_res(t.t_total),'time','f res',filename='f_res')

ploting(t.t_total,wl_scan,'time','lambda res (um)',filename='lambda_res')
ploting(t.t_total,Q,'time (ps)','Q',filename='Q_test')
ploting(t.t_total,abs(s_minus)**2,'time (ps)',r"$|s_-^2|$",filename='s_minus_power')
ploting(t.t_total,180/np.pi*np.angle(s_minus),'time (ps)','s_minus phase',filename='s_minus phase')

ploting(t.t_total,abs(b)**2,'time (ps)','b',filename='b_test')
ploting(wl*1000, data, 'wavelength (nm)', 'Transfer function of ring', filename="T")
np.savetxt("T.txt", data, fmt="%.8f", delimiter="\n")
np.savetxt("wl.txt",wl, fmt="%.8f", delimiter="\n")
np.savetxt('time',t.t_total,fmt="%.8f", delimiter="\n")
np.savetxt('wl_scan',wl_scan ,fmt="%.8f", delimiter="\n")
np.savetxt("s_minus_power.txt", abs(s_minus)**2, fmt="%.8f", delimiter="\n")
os.chdir("../")
# 249189
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ring_mod = ring(L=2*np.pi*5, 
#             neff=2.51464,
#             ng=4.3, 
#             gamma = 0.95012, 
#             alpha = 0.95246,
#             me=37.2,
#             FCA_coeff_ratio=0)
# ring_mod.scan_frequency(wl_min ,wl_max,1.5488,t)


# wl_scan =  c/ring_mod.w_res(t.t_total)*t0
# os.chdir("./data/")
# ploting(t.t_total,ring_mod.w_res(t.t_total),'time','f res',filename='f_res')
# ploting(t.t_total,wl_scan,'time','lambda res (um)',filename='lambda_res')
# b,Q,s_minus = solving(ring_mod,v,t,wl_in,Pin=Pin)

# T = Transfer_function(ring_mod,t)
# wl,data = T.mapping(abs(s_minus)**2/Pin)

# ploting(t.t_total,Q,'time (ps)','Q',filename='Q_test')
# ploting(t.t_total,abs(s_minus)**2/Pin,'time (ps)',r"$|s_-^2|$",filename='s_minus_power_NL_0')
# ploting(t.t_total,180/np.pi*np.angle(s_minus),'time (ps)','s_minus phase',filename='s_minus phase_NL_0')

# ploting(t.t_total,abs(b)**2/Pin,'time (ps)','b',filename='b_test_NL_0')
# ploting(wl*1000, data, 'wavelength (nm)', 'Transfer function of ring', filename="T_NL_0")
# os.chdir("../")





# # np.savetxt("b_bar.txt", b/b0, fmt="%.8f", delimiter="\n")
# # np.savetxt("Q_bar.txt", Q/v.Cj, fmt="%.8f", delimiter="\n")
# # np.savetxt("s_minus_power.txt", abs(s_minus)**2, fmt="%.8f", delimiter="\n")
# # np.savetxt("s_minus_phase.txt", 180/np.pi*np.angle(s_minus), fmt="%.8f", delimiter="\n")

