from utility import ploting
import numpy as np
from scipy.interpolate import interp1d
import scipy.io
import subprocess
import importlib
import argparse

lumapi = importlib.machinery.SourceFileLoader('Lumapi','/opt/lumerical/v212/api/python/lumapi.py').load_module()

fdtd = lumapi.FDTD(hide=True)

def main(radius,
         gap,
         z_span=0.22,
         wg_width=0.45,
         from_wg=1,
         from_ring=0):
    run = False
    scale = 1e-6
    pad = 2
    res = 5
    Si = 'Si (Silicon) - Palik'
    SiO2 = 'SiO2 (Glass) - Palik'

    fdtd_x_span = (2*radius+wg_width+2*pad)*scale
    fdtd_y_min = (-gap/2-wg_width-pad)*scale
    fdtd_y_max = (gap/2+wg_width+radius+2*pad)*scale
    fdtd_z_span = (3)*scale
    simulation_time = 1000e-15

    fdtd.addring(name="new_ring")
    fdtd.set("x",0)
    fdtd.set("y",(radius+wg_width/2+gap/2)*scale)
    fdtd.set("inner radius",(radius - wg_width/2)*scale)
    fdtd.set("outer radius",(radius + wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",z_span*scale)
    fdtd.set("theta start",180)
    fdtd.set("theta stop",360)
    fdtd.set('material',Si)

    fdtd.addrect(name="bus")
    fdtd.set("x",0);
    fdtd.set("y",(-gap/2-wg_width/2)*scale);
    fdtd.set("z",0);
    fdtd.set("x span",fdtd_x_span+2*pad*scale);
    fdtd.set("y span",wg_width*scale);
    fdtd.set("z span",z_span*scale);
    fdtd.set('material',Si)

    fdtd.addfdtd(simulation_time=simulation_time)
    fdtd.set('dimension','3D')
    fdtd.set('x',0)
    fdtd.set('x span',fdtd_x_span)
    fdtd.set('y',0)
    fdtd.set('y min',fdtd_y_min)
    fdtd.set('y max',fdtd_y_max)
    fdtd.set('z',0)
    fdtd.set('z span',fdtd_z_span)
    fdtd.set('mesh accuracy',res)
    fdtd.set("background material",SiO2)
    
    fdtd.addport ()
    fdtd.set("name","bottom_left")
    fdtd.set("direction","forward")
    fdtd.set("mode selection","user select")
    fdtd.set("number of field profile samples",3)
    fdtd.set("bent waveguide",True)
    fdtd.set("bend radius",radius*scale)
    fdtd.set("theta",-45)
    fdtd.updateportmodes(1) 
    fdtd.set("x",-(radius/(2)**0.5)*scale)
    fdtd.set("y",((radius+wg_width/2+gap/2)-radius/(2)**0.5)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.set("x span",(0)*scale)
    fdtd.select("FDTD::ports") 

    fdtd.addport ()
    fdtd.set("name","bottom_right")
    fdtd.set("direction","backward")
    fdtd.set("mode selection","user select")
    fdtd.set("number of field profile samples",3)
    fdtd.set("bent waveguide",True)
    fdtd.set("bend radius",radius*scale)
    fdtd.set("theta",45)
    fdtd.updateportmodes(1)
    fdtd.set("x",(radius/(2)**0.5)*scale)
    fdtd.set("y",((radius+wg_width/2+gap/2)-radius/(2)**0.5)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.set("x span",(0)*scale)
    fdtd.select("FDTD::ports")
    
    fdtd.addport ()
    fdtd.set("name","top_left")
    fdtd.set("direction","forward")
    fdtd.set("mode selection","user select")
    fdtd.set("number of field profile samples",3)
    fdtd.updateportmodes(1) 
    fdtd.set("x",-(radius+wg_width/2)*scale)
    fdtd.set("y",(-gap/2-wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.set("x span",(0)*scale)
    fdtd.select("FDTD::ports") 

    fdtd.addport ()
    fdtd.set("name","top_right")
    fdtd.set("direction","backward")
    fdtd.set("mode selection","user select")
    fdtd.set("number of field profile samples",3)
    fdtd.updateportmodes(1) 
    fdtd.set("x",(radius/2+wg_width/2)*scale)
    fdtd.set("y",(-gap/2-wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.set("x span",(0)*scale)
    fdtd.select("FDTD::ports") 
    
    if from_ring:
        fdtd.set("source port","bottom_left")
    if from_wg:
        fdtd.set("source port","top_left")
    fdtd.set("source mode","mode 1")

    fdtd.save("coupler_fsp")

    if run:
        print('complete build\nrunning...........\n')
        process0 = subprocess.Popen(['/opt/lumerical/v212/mpich2/nemesis/bin/mpiexec --hostfile host_file /opt/lumerical/v212/bin/fdtd-engine-mpich2nem '+filename],shell=True)
        process0.wait()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("radius",help="specify radius of ring",type=np.float64)
    parser.add_argument("gap",help="specify gap between ring and bus",type=np.float64)
    arg = parser.parse_args()
    main(arg.radius,arg.gap)