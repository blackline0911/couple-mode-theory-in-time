import numpy as np
import importlib
import argparse

lumapi = importlib.machinery.SourceFileLoader('Lumapi','/opt/lumerical/v212/api/python/lumapi.py').load_module()

fdtd = lumapi.FDTD(hide=True)

def main(radius):
    scale = 1e-6
    WG_height = 0.2114
    wg_width = 0.5
    pad = 2
    Si = 'Si (Silicon) - Palik'
    SiO2 = 'SiO2 (Glass) - Palik'
    SKT_height  = 0.0622
    L = 5
    res = 4

    fdtd.setglobalsource("wavelength start",1.5e-6)
    fdtd.setglobalsource("wavelength start",1.6e-6)
    fdtd.setglobalmonitor("frequency points",1000)

    fdtd_x_span = (radius+wg_width+6)*scale
    fdtd_y_span = (radius+wg_width+L+2)*scale
    fdtd_z_span = (3)*scale
    fdtd_x = radius/2/2**0.5*scale+3*scale
    fdtd_y = radius/2/2**0.5*scale+2*scale
    simulation_time = 1000e-15

    fdtd.addrect(name="slab")
    fdtd.set("x",fdtd_x)
    fdtd.set("y",fdtd_y)
    fdtd.set("z",0)
    fdtd.set("x span",fdtd_x_span+2*pad*scale)
    fdtd.set("y span",(fdtd_y_span)+2*pad*scale)
    fdtd.set("z",(-WG_height/2+SKT_height/2)*scale)
    fdtd.set("z span",SKT_height*scale)
    fdtd.set('material',Si)

    fdtd.addrect(name="input wg")
    fdtd.set("x",-L/2*scale)
    fdtd.set("y",radius*scale)
    fdtd.set("z",0)
    fdtd.set("x span",L*scale)
    fdtd.set("y span",wg_width*scale)
    fdtd.set("z span",WG_height*scale)
    fdtd.set('material',Si)

    fdtd.addrect(name="output wg")
    fdtd.set("x",radius*scale)
    fdtd.set("y",-L/2*scale)
    fdtd.set("z",0)
    fdtd.set("x span",wg_width*scale)
    fdtd.set("y span",L*scale)
    fdtd.set("z",0)
    fdtd.set("z span",WG_height*scale)
    fdtd.set('material',Si)

    fdtd.addring(name="new_ring")
    fdtd.set("x",0)
    fdtd.set("y",0)
    fdtd.set("inner radius",(radius - wg_width/2)*scale)
    fdtd.set("outer radius",(radius + wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",WG_height*scale)
    fdtd.set("theta start",0)
    fdtd.set("theta stop",90)
    fdtd.set('material',Si)

    fdtd.addmode()
    fdtd.set("name","mode source")
    fdtd.set("injection axis","x")
    fdtd.set("x",0)
    fdtd.set("y",radius*scale)
    fdtd.set("z",0)
    fdtd.set("y span",3*scale)
    fdtd.set("z span",2*scale)
    fdtd.set("mode selection","user select")
    fdtd.set("bent waveguide",True)
    fdtd.set("bend radius",radius*scale)
    fdtd.set("center wavelength",1.55*scale)
    fdtd.set("wavelength span",0.1*scale)
    fdtd.set("mode selection","user select")
    fdtd.set('mode selection','user select')
    fdtd.set('selected mode number',1)
    
    fdtd.addpower()
    fdtd.set("name","field xy")
    fdtd.set("monitor type","2D Z-normal")
    fdtd.set("x",fdtd_x)
    fdtd.set("y",fdtd_y)
    fdtd.set("z",0)
    fdtd.set("x span",fdtd_x_span)
    fdtd.set("y span",fdtd_y_span)

    fdtd.addpower()
    fdtd.set("name","dft_output")
    fdtd.set("monitor type","2D Y-normal")
    fdtd.set("x",radius*scale)
    fdtd.set("y",0)
    fdtd.set("z",0)
    fdtd.set("x span",3*scale)
    fdtd.set("z span",2*scale)
    

    fdtd.addmodeexpansion()
    fdtd.set("mode selection","user select")
    fdtd.set("name","mode_output")
    fdtd.set("monitor type","2D Y-normal")
    fdtd.set("x",radius*scale)
    fdtd.set("y",0)
    fdtd.set("z",0)
    fdtd.set("x span",3*scale)
    fdtd.set("z span",2*scale)
    fdtd.set("bent waveguide",True)
    fdtd.set("bend radius",radius*scale)
    fdtd.set("frequency points",3)
    fdtd.setexpansion('through port','dft_output')
    fdtd.set('mode selection','fundamental TE mode')

    fdtd.addfdtd(simulation_time=simulation_time)
    fdtd.set('dimension','3D')
    fdtd.set("x",fdtd_x)
    fdtd.set("y",fdtd_y)
    fdtd.set('z',0)
    fdtd.set('x span',fdtd_x_span)
    fdtd.set('y span',fdtd_y_span)
    fdtd.set('z span',fdtd_z_span)
    fdtd.set('mesh accuracy',res)
    fdtd.set("background material",SiO2)

    fdtd.save("bend_loss_radius"+str(radius))
    fdtd.close()
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("radius",help="specify radius of ring",type=np.float64)
    arg = parser.parse_args()
    main(arg.radius)
