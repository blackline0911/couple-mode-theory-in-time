import numpy as np
import importlib
import argparse
import subprocess
import matplotlib.pyplot as plt

lumapi = importlib.machinery.SourceFileLoader('Lumapi','/opt/lumerical/v212/api/python/lumapi.py').load_module()


def main(radius,
         gap,
         wg_width=0.5,
         bus_wg_width=0.479,
         WG_height=0.2114,
         SKT_height=0.0622
         ):
    fdtd = lumapi.FDTD(hide=True)
    scale = 1e-6
    pad = 2
    res = 4
    Si = 'Si (Silicon) - Palik'
    SiO2 = 'SiO2 (Glass) - Palik'
    fdtd.switchtolayout()
    fdtd.setglobalsource("wavelength start",1.5e-6)
    fdtd.setglobalsource("wavelength start",1.6e-6)
    fdtd.setglobalmonitor("frequency points",500)

    fdtd_x_span = (2*radius+wg_width+2*pad)*scale
    fdtd_y_min = (-gap/2-wg_width-pad)*scale
    fdtd_y_max = (gap/2+radius-1)*scale
    fdtd_z_span = (3)*scale
    simulation_time = 1000e-15

    fdtd.addrect(name="slab")
    fdtd.set("x",0)
    fdtd.set("y",0)
    fdtd.set("z",0)
    fdtd.set("x span",fdtd_x_span+2*pad*scale)
    fdtd.set("y max",(gap/2+wg_width+radius+pad)*scale)
    fdtd.set("y min",(fdtd_y_min))
    fdtd.set("z",(-WG_height/2+SKT_height/2)*scale)
    fdtd.set("z span",SKT_height*scale)
    fdtd.set('material',Si)

    fdtd.addring(name="new_ring")
    fdtd.set("x",0)
    fdtd.set("y",(radius+wg_width/2+gap/2)*scale)
    fdtd.set("inner radius",(radius - wg_width/2)*scale)
    fdtd.set("outer radius",(radius + wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",WG_height*scale)
    fdtd.set("theta start",180)
    fdtd.set("theta stop",360)
    fdtd.set('material',Si)

    fdtd.addrect(name="bus")
    fdtd.set("x",0) 
    fdtd.set("y",(-gap/2-0.45/2)*scale) 
    fdtd.set("z",0) 
    fdtd.set("x span",fdtd_x_span+2*pad*scale) 
    fdtd.set("y span",bus_wg_width*scale) 
    fdtd.set("z span",WG_height*scale) 
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

    fdtd.addpower()
    fdtd.set("name","bottom_left")
    fdtd.set("monitor type","2D X-normal")
    fdtd.set("x",-(radius/(2)**0.5)*scale)
    fdtd.set("y",((radius+wg_width/2+gap/2)-radius/(2)**0.5)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)


    fdtd.addmodeexpansion()
    fdtd.set("mode selection","user select")
    fdtd.set("name","mode_output_bottom_left")
    fdtd.set("monitor type","2D X-normal")
    fdtd.set("mode selection","user select")
    fdtd.set("bent waveguide",True)
    fdtd.set("bend radius",radius*scale)
    fdtd.set("theta",-45)
    fdtd.set("frequency points",3)
    fdtd.set("x",-(radius/(2)**0.5)*scale)
    fdtd.set("y",((radius+wg_width/2+gap/2)-radius/(2)**0.5)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.setexpansion("bottom_left","bottom_left")
    fdtd.set('mode selection','user select')
    fdtd.set('selected mode numbers',1)
    fdtd.set("frequency points",3)


    fdtd.addpower()
    fdtd.set("name","bottom_right")
    fdtd.set("monitor type","2D X-normal")
    fdtd.set("x",(radius/(2)**0.5)*scale)
    fdtd.set("y",((radius+wg_width/2+gap/2)-radius/(2)**0.5)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)


    fdtd.addmodeexpansion()
    fdtd.set("mode selection","user select")
    fdtd.set("name","mode_output_bottom_right")
    fdtd.set("monitor type","2D X-normal")
    fdtd.set("mode selection","user select")
    fdtd.set("bent waveguide",True)
    fdtd.set("bend radius",radius*scale)
    fdtd.set("theta",45)
    fdtd.set("frequency points",3)
    fdtd.set("x",(radius/(2)**0.5)*scale)
    fdtd.set("y",((radius+wg_width/2+gap/2)-radius/(2)**0.5)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.setexpansion("bottom_right","bottom_right")
    fdtd.set('mode selection','user select')
    fdtd.set('selected mode numbers',1)

    fdtd.addmode()
    fdtd.set("name","mode source")
    fdtd.set("injection axis","x")
    fdtd.set("x",-(radius+wg_width/2)*scale)
    fdtd.set("y",(-gap/2-wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.set("mode selection","user select")
    fdtd.set("center wavelength",1.55*scale)
    fdtd.set("wavelength span",0.1*scale)
    fdtd.set("mode selection","user select")
    fdtd.set('selected mode number',1)


    fdtd.addpower()
    fdtd.set("name","dft_output")
    fdtd.set("monitor type","2D X-normal")
    fdtd.set("x",(radius+wg_width/2)*scale)
    fdtd.set("y",(-gap/2-wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)

    fdtd.addmodeexpansion()
    fdtd.set("mode selection","user select")
    fdtd.set("name","mode_output")
    fdtd.set("monitor type","2D X-normal")
    fdtd.set("x",(radius+wg_width/2)*scale)
    fdtd.set("y",(-gap/2-wg_width/2)*scale)
    fdtd.set("z",0)
    fdtd.set("z span",(3)*scale)
    fdtd.set("y span",(3)*scale)
    fdtd.set("mode selection","user select")
    fdtd.set("frequency points",3)
    fdtd.setexpansion('through port','dft_output')
    fdtd.set('mode selection','user select')
    fdtd.set('selected mode numbers',1)


    fdtd.save("coupler_"+str(int(radius))+"radius")


if __name__=='__main__':
    ifscan = False
    scan_number = 5
    scan = ["bus_wg_width","SKT_height","WG_height","radius","gap"]
    parser = argparse.ArgumentParser()
    parser.add_argument("radius",help="specify radius of ring",type=np.float64)
    parser.add_argument("gap",help="specify gap between ring and bus",type=np.float64)
    arg = parser.parse_args()

    print("......starting building......")
    if ifscan:
        N = 7
        if scan_number == 1:
            data = np.linspace(0.21,0.213,N)
        if scan_number == 2:
            data = np.linspace(0.052,0.072,N)
        if scan_number == 3:
            data = np.linspace(0.46,0.49,N)
        if scan_number == 4:
            data = np.linspace(5,10,N)
        if scan_number == 5:
            data = np.linspace(0.3,0.3,N)
                #data = np.linspace(0.1,0.2,N)

    print('\tcomplete building\n\nrunning simulation......\n')
    print(scan)
    if ifscan:
        plt.figure()
        for i in range(len(data)):
            print(f"\t\nRunning simulation with {scan[scan_number-1]} = {data[i]} um\n")
            if scan_number == 1:
                    main(arg.radius,arg.gap,bus_wg_width=data[i])
                    radius = arg.radius
            if scan_number == 2:
                    main(arg.radius,arg.gap,SKT_height=data[i])
                    radius = arg.radius
            if scan_number == 3:
                    main(arg.radius,arg.gap,WG_height=data[i])
                    radius = arg.radius
            if scan_number == 4:
                    main(data[i],arg.gap)
                    radius = data[i]
            if scan_number == 5:
                main(arg.radius,data[i])
                radius = arg.radius
            p1 = subprocess.Popen(['/opt/lumerical/v212/mpich2/nemesis/bin/mpiexec --hostfile host_file /opt/lumerical/v212/bin/fdtd-engine-mpich2nem '+"coupler_"+str(int(radius))+"radius.fsp"],shell=True)
            p1.wait()

            fdtd = lumapi.FDTD(hide=True)
            fdtd.load("coupler_"+str(int(radius))+"radius")
            output = fdtd.getresult("mode_output","expansion for through port")
            a = output['a']
            f = output['f']
        plt.plot(299792458*1e6/f,np.abs(a),label="%s = %.4f um"%(scan[scan_number-1],data[i]))
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Amplitude")
        plt.title("Mode Expansion Amplitude")
        plt.legend()
        plt.grid() 
        plt.savefig(f"coupler_region_{scan[scan_number-1]}.png")
        plt.show()
    else:
        main(arg.radius,arg.gap)

        
        