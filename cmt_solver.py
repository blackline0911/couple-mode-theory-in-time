# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
from scipy.integrate import *
import numpy as np
import matplotlib.pyplot as plt
from cmath import *
import math
from scipy.signal import *
from scipy import signal
from random import *
from utility import *
from user import *

# relative solver tolerance 
rtol = 1e-14
# absolute solver tolerance
atol = 1e-20
# accuracy = atol + abs(y)*rtol

# method of solving algorithm
method = 'RK45'

def solving(ring, 
            driver,
            time,
            lambda_incident,
            Pin,
            ):
    """
    ring : ring object
    driver: driver object
    time: time object
    lambda_incident: incident laser wavelength (micron)
    Pin : Power of input laser (Watt)
    """
    w_pround = 2*np.pi*c/(lambda_incident*1e-6)
    # normalize factors
    S0 = sqrt(Pin)
    b0 = sqrt(t0)*S0
    tu_t_bar = ring.tu_t/t0
    tu_e_bar = ring.tu_e/t0

    # normalized dw/dlambda
    D_bar = -2*np.pi*c/ring.lambda0**2 * t0

    # normalized frequency
    w0_bar = ring.w_res*t0
    w_pround_bar = w_pround*t0

    def CMT(t_bar,eqs):
        b_bar , Q_pround = eqs
        
        f1 = 1j*(w0_bar-w_pround_bar)*b_bar- (1/tu_t_bar)*b_bar + sqrt(2/tu_e_bar) *1 + 1j*D_bar*ring.me*Q_pround*b_bar

        f2 = (driver.v[t_bar]/(driver.R * driver.Cj_dict[driver.v[t_bar]])*t0) - (1/( driver.R*driver.Cj_dict[ driver.v[t_bar] ] ) )*Q_pround*t0
        return [ f1,f2 ]

#     """Since the length of solution function array may not be the same as t_eval argument we specified when the time in single solve_ivp is long.
#     Hence, we divide the time according to the Baud Rate, and solve coupled differential equation by each time segments. 
#     """
    b_record = np.array([])
    Q_record = np.array([])
    b_init=0+1j*0
    Q_init=0
    for i in range(time.N):
        if i==0:
            sol =  solve_ivp(CMT ,[0+time.T_normalized*i ,time.T_normalized*(i+1)], [b_init, Q_init],method=method,t_eval = time.t_all_segment[i] ,atol = atol,rtol = rtol)
            b = sol.y[0]
            b_record = np.append(b_record,b)
            q = sol.y[1]
            Q_record = np.append(Q_record,q)
        else:
            sol =  solve_ivp(CMT ,[0 ,time.T_normalized*(i+1)], [b_init, Q_init],method=method,t_eval = np.append( np.array( [time.t_all_segment[i][0]-time.dt/t0]), time.t_all_segment[i]  ),atol = atol,rtol = rtol)
            b = sol.y[0]
            b_record = np.append(b_record,b[1::])
            q = sol.y[1]
            Q_record = np.append(Q_record,q[1::])
        b_init = b[-1]
        Q_init = q[-1]
    # sol = solve_ivp(CMT ,[0 ,t_max], [ring.b_init, driver.Q_init],method=method,t_eval = t_total,atol = atol,rtol = rtol)
    b_bar = (sol.y[0])
    Q_bar = (sol.y[1])
    s_minus_bar = (1-sqrt(2/tu_e_bar)*b_bar)
    return b_bar*b0, Q_bar*driver.Cj , s_minus_bar*S0