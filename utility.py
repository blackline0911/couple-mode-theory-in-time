from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import os

c = 299792458*1e6
t0 = 1e-12
rtol = 1e-15
atol = 1e-20
h=6.626e-34

def ploting(x,*arg,x_label, title,filename='',figcolor='w',line_color='b',
            grid_color='g',grid_style='--',grid_alpha=0.5):
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
       plt.figure()
       for i in arg:
           plt.plot(x,i)
       plt.xlabel(x_label)
       plt.title(title)
       plt.grid(color=grid_color,linestyle=grid_style, alpha=grid_alpha)
       ax = plt.gca()
       ax.set_facecolor(figcolor)
       if(filename!=''):
              plt.savefig(filename)
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

if __name__=='__main__':
    x = np.linspace(0,10,1000)
    y = x**2
    ploting(x,y,r"$\alpha$",r'$|\alpha|^2$','testing.png')