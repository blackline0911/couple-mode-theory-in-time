import numpy as np
from scipy.integrate import *


w0 = 2*np.pi*1
t = np.linspace(0,10,10000)
def f(t,x):
    # print("t = ",t)
    return 1j*w0*x - x + 1 - abs(x)**2*x + 1j*abs(x)**2*x
# def f(t,x):
#     return 1j*w0*x - x + 1

x_ans = 1/(1j*w0-1)*np.exp((1j*w0-1)*t)-1/(1j*w0-1)


import matplotlib.pyplot as plt
plt.figure() 

plt.plot(t,x_ans,label='Analytic solution')
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "RK45")
plt.plot(t,sol.y[0],label='RK45')
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "RK23")
plt.plot(t,sol.y[0],label='RK23')
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "DOP853")
plt.plot(t,sol.y[0],label='DOP853')
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "BDF")
plt.plot(t,sol.y[0],label='BDF')
plt.legend()
plt.show()

plt.figure()
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "RK45")
plt.plot(t,abs(x_ans-sol.y[0]),label='RK45')
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "RK23")
plt.plot(t,abs(x_ans-sol.y[0]),label='RK23')
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "DOP853")
plt.plot(t,abs(x_ans-sol.y[0]),label='DOP853')
sol =solve_ivp(f,t_span=[0,10],y0=[0+1j*0],t_eval=t,method = "BDF")
plt.plot(t,abs(x_ans-sol.y[0]),label='BDF')
plt.legend()
plt.title('Error')
plt.show()
