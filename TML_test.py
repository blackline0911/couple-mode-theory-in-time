from scipy.integrate import *
import numpy as np
from utility import *
from driver import driver


v_bias = -1
vpp = 1
f_drive = 50
R = 53.9

v = driver()
Cj = 20e-15
Rs = 59.3
Rsi = 1439.8
Cox = 34.7e-15
Cp = 6.6e-15
f0 = 1e12

Cp_bar = Cp*f0
Cox_bar = Cox*f0
Cj_bar = Cj*f0

