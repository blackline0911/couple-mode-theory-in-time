from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

c = 299792458*1e6
t0 = 1e-12
rtol = 1e-13
atol = 1e-13
# rtol = 1e-15
# atol = 1e-20
h=6.626e-34

def show(filename):
    """
    Showing matplotplib.pyplot figure saved by picke.dump()
    """
    f = open(filename,'rb')
    pickle.load(f)
    plt.show()
    f.close()
    return 

def ploting(x,*arg,x_label, title,filename='',figcolor='w',line_color='b',
            grid_color='g',grid_style='--',grid_alpha=0.5,leg=['']):
       """
       input argments:
              x:x
              y:data you want to plot
              x_label:x_label
              title: figure title
              figcolor: set the color in figure
              line_color:set the color of line you plot
              grid_color: set the color of grid
              grid_style: set the style of grid
              grid_alpha: set opacity of grid
              filename: set the file name of the figure you plot. 
                        If unspecified, the figure will not be saved.
       """
       fig=plt.figure()
       n=0
       for i in arg:
            if (not leg==['']):
                plt.plot(x,i,label=leg[n])
                plt.legend()
            else:
                plt.plot(x,i)
            n+=1
       plt.xlabel(x_label)
       plt.title(title)
       plt.grid(color=grid_color,linestyle=grid_style, alpha=grid_alpha)
       ax = plt.gca()
       ax.set_facecolor(figcolor)
       if(filename!=''):
            plt.savefig(filename)
            # with open(filename, "wb") as f:
            #     pickle.dump(fig, f)
       plt.show()
       return

def sinc(t):
        if isinstance(t, np.ndarray):
            return np.where(t==0,1,np.sin(np.pi*t)/(np.pi*t))
        else:
            if t==0:
                return 1
            else:
                return np.sin(np.pi*t)/(np.pi*t)

def dB(x):
     return 10*np.log10(x)

def dB_inv(x):
    return 10**(-x/10)

def TP(low_level,high_level):
    return -10*np.log10( abs( high_level - low_level )/2 )

def ER(time,low_level,high_level):
    ER = np.zeros( int(len(time.t_total)-time.buffer*t0/time.dt-1))
    N = len(ER)
    for i in range(N):
        if high_level>low_level:
                ER[i] = -10*np.log10(low_level/high_level)
        else:
            ER[i] = -10*np.log10(high_level/low_level)
    return ER

def FDM(dt,signal):
    N = len(signal)
    FDM = np.zeros(N)
    FDM[0] = (-3*signal[0] + 4*signal[1] - signal[2])/2/dt
    FDM[-1] = (3*signal[-1] - 4*signal[-2] + signal[-3])/2/dt
    for i in range(1,N-1):
        FDM[i] = (signal[i+1]-signal[i-1])/2/dt
    return FDM