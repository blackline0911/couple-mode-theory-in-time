import numpy as np
cimport numpy as cnp
from cpython.array cimport array

cnp.import_array()

DTYPE = np.float64
ctypedef cnp.float64_t data_type

cpdef object sinc(x):
    cdef data_type p = np.pi
    if isinstance(x, cnp.ndarray):
        return np.where(x==0,1,np.sin(p*x)/(p*x))
    else:
        if x==0:
            return 1
        else:
            return np.sin(p*x)/(p*x)

cpdef data_type raise_cos(BS,t,shift,beta,T):
    cdef array ans = 0
    if isinstance(t,cnp.ndarray):
        ans = np.where( ( (t-shift)==T/2/beta) | ( (t-shift)==( -T/2/beta) ),
                        BS[int(shift/T)-1]*np.pi/4*sinc((1/2/beta)),
                        BS[int(shift/T)-1] * sinc((t-shift)/T) * np.cos(np.pi*beta*(t-shift)/T) / (1-(2*beta*(t-shift)/T)**2) )
        return ans
    else:
        return BS[int(shift/T)-1]*sinc((t-shift)/T)*np.cos(np.pi*beta*(t-shift)/T)/(1-(2*beta*(t-shift)/T)**2)

cpdef 