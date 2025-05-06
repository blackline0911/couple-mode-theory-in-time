import numpy as np
from utility import ploting
from matplotlib.pyplot import *
def sinc(t):
    if isinstance(t,np.ndarray):
        return np.where(t==0, 1, np.sin(np.pi*t)/(np.pi*t))
    else:
        if t==0:
            return 1
        else:
            return np.sin(np.pi*t)/(np.pi*t)

def rcos(t,T,beta=1):
    if isinstance(t,np.ndarray):
        return np.where( ((t==T/2/beta) | (t==-T/2/beta)), \
                        np.pi/4/T*sinc(1/2/beta), \
                        1/T*sinc(t/T)*np.cos(np.pi*beta*t/T)/(1-(2**beta*t/T)**2) )
    else:
        if (t==T/2/beta) | (t==-T/2/beta):
            return np.pi/4/T*sinc(1/2/beta)
        else:
            return 1/T*sinc(t/T)*np.cos(np.pi*beta*t/T)/(1-(2*beta*t/T)**2)

figure()
tmax = 10
dt = 1e-5
t0 = np.arange(0,tmax,dt)
f_precise = rcos(t0,1)
# plot(t0,f_precise)

dt = 1e-2
t1 = np.arange(0,tmax,dt)
# scatter(t1,rcos(t1,1),marker='o',c='r')
t2 = np.arange(0,tmax,1e-3)
linear_interp = np.interp(t2,t1,rcos(t1,1))
# scatter(t2,linear_interp,marker='o',c='b')
# xlabel("time")
# title("raise cosine")
# show()


# error = np.zeros(len(t2))
# for i in range(len(t2)):
#     error[i] = abs( linear_interp[i] - f_precise[np.argmin(abs(t0-t2[i]))]) 
# plot(t2,error,label='Linear interpolation')

from scipy.interpolate import CubicSpline
cubic_interp = CubicSpline(t1,rcos(t1,1))
error = np.zeros(len(t2))
for i in range(len(t2)):
    error[i] = abs( cubic_interp(t2[i]) - f_precise[np.argmin(abs(t0-t2[i]))]) 
plot(t2,error,label='Cubic interpolation')


from scipy.interpolate import make_interp_spline
k=3
make_interp_spline_interp = make_interp_spline(t1,rcos(t1,1),k=k)
error = np.zeros(len(t2))
for i in range(len(t2)):
    error[i] = abs(make_interp_spline_interp(t2[i]) - f_precise[np.argmin(abs(t0-t2[i]))]) 
plot(t2,error,label='('+str(k)+'-1)th derivative ')
legend()
show()

