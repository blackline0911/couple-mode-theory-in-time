import numpy as np
cimport numpy as cnp
from cpython.array cimport array
from scipy.integrate import *

cnp.import_array()

DTYPE = np.float64
ctypedef cnp.float64_t float_data_type
ctypedef cnp.int64_t int_data_type
ctypedef cnp.ndarray array_type

cpdef object sinc(x):
    cdef float_data_type p = np.pi
    if isinstance(x, cnp.ndarray):
        return np.where(x==0,1,np.sin(p*x)/(p*x))
    else:
        if x==0:
            return 1
        else:
            return np.sin(p*x)/(p*x)

cpdef raise_cos(driver,cnp.ndarray[float_data_type, ndim=1] BS,t,shift,beta,T):
    if isinstance(t,cnp.ndarray):
        #cdef cnp.ndarray[float_data_type, ndim=1] ans = np.zeros(len(t),dtype=DTYPE)
        ans = np.where( ( (t-shift)==T/2/beta) | ( (t-shift)==( -T/2/beta) ),
                        driver.vpp*int(BS[int(shift/T)])/(driver.level_num-1)*np.pi/4*sinc((1/2/ beta)),
                        driver.vpp*int(BS[int(shift/T)])/(driver.level_num-1) * sinc((t-shift)/T) * np.cos(np.pi*beta*(t-shift)/T) / (1-(2*beta*(t-shift)/T)**2) )
        return ans
    else:
        if ( (t-shift)==T/2/beta) | ( (t-shift)==( -T/2/beta) ):
            s = driver.vpp*BS[int(shift/T)]*np.pi/4*sinc((1/2/beta))
        else:
            s = driver.vpp*BS[int(shift/T)]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)
        return s
cpdef create_rcos_signal(cnp.ndarray[float_data_type, ndim=1] BS,
                        cnp.ndarray[float_data_type, ndim=1] t,
                        T_normalized,
                        N,
                        driver,
                        bet):
    cdef float_data_type beta = bet
    cdef cnp.ndarray[float_data_type, ndim=1] rcos_signal = np.zeros(len(t))
    for i in range(N):
        #passed_T_num = int(t[i]//T_normalized)
        #for j in range(N):
        rcos_signal += raise_cos(driver,BS,t,(i)*T_normalized,beta,T_normalized)
    return rcos_signal

cpdef refering_rcos_signal(cnp.ndarray[float_data_type, ndim=1] BS,
                        t,
                        T_normalized,
                        N,
                        driver,
                        bet):
    cdef float_data_type beta = bet
    voltage = 0
    for i in range(N):
        voltage += raise_cos(driver,BS,t,(i)*T_normalized,beta,T_normalized)
    return voltage