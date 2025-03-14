"""
This script is used to save the simulation variable that A simulation used. 
Just for Future checking or further research
"""

from test import *
import numpy as np
import ring, driver, time_class, cmt_solver ,utility

filename = 'sim.txt'

os.chdir("./data/")
with open(filename, 'w') as f:
    f.write('This script is used to save the simulation variable that A simulation used. \nJust for Future checking or further research\n\n')

    f.write("Parameter of ring modulator\n")
    f.write("\tStructure parameter : \n")
    f.write("\t\tcavity length : ")
    f.write(str(ring_mod.L)+' um')
    f.write('\n')

    f.write("\t\tneff : ")
    f.write(str(ring_mod.neff))
    f.write('\n')

    f.write("\t\tGroup index : ")
    f.write(str(ring_mod.ng))
    f.write('\n')

    f.write("\t\tcouple through coefficient : ")
    f.write(str(ring_mod.gamma))
    f.write('\n')

    f.write("\t\tLinear attenuation coefficient : ")
    f.write(str(ring_mod.alpha))
    f.write('\n')

    f.write("\t\tmodulation efficiency : ")
    f.write(str(ring_mod.me) + ' pm/V')
    f.write('\n')

    f.write("\t\tFCA coefficient ratio : ")
    f.write(str(ring_mod.FCA_coeff_ratio))
    f.write('\n')

    f.write("\tPhysical parameter : \n")
    f.write("\t\tintrinsic photon life time : ")
    f.write(str(np.real(ring_mod.tu_o_bar))+' ps')
    f.write('\n')

    f.write("\t\tcouple extrinsic photon life time : ")
    f.write(str(np.real(ring_mod.tu_e_bar))+' ps')
    f.write('\n')

    f.write("\t\tresonant frequency : ")
    f.write(str(ring_mod.f_res_bar)+' THz')
    f.write('\n')

    f.write("\t\tresonant wavelength : ")
    f.write(str(ring_mod.lambda0)+' um')
    f.write('\n')

    f.write("\t\tQ factor : ")
    f.write(str(ring_mod.Q))
    f.write('\n\n')

    m = t.mode
    f.write("Current Simulation mode is "+t.mode+'\n\n')
    if m == "voltage_drive":
        f.write("\tDriver information : \n")
        if v.sine_wave:
            f.write("\tDriving Signal is sine wave\n")
            f.write("\t\tDriving frequency is "+str(v.f_drive/1e9)+" GHz\n")
            f.write("\t\tBias voltage is "+str(v.v_bias)+" V\n")
            f.write("\t\tPeak to Peak voltage is "+str(v.vpp)+" V\n\n")
        if v.square_wave:
            f.write("\tDriving Signal is square wave\n")
            f.write("\t\tDriving frequency is "+str(v.f_drive/1e9)+" GHz\n")
            f.write("\t\tBias voltage is "+str(v.v_bias)+" V\n")
            f.write("\t\tPeak to Peak voltage is "+str(v.vpp)+" V\n\n")
        if v.raise_cosine:
            f.write("\tDriving Signal is raise_cosine signal\n")
            f.write("\t\tDriving frequency is "+str(v.f_drive/1e9)+" GHz\n")
            f.write("\t\tBias voltage is "+str(v.v_bias)+" V\n")
            f.write("\t\tPeak to Peak voltage is "+str(v.vpp)+" V\n\n")
    if m == "scan_frequency":
        f.write("\tThe start wavelength is "+str(c/ring_mod.f_start_bar*t0)+" um\n")
        f.write("\tThe end wavelength is "+str(c/ring_mod.f_end_bar*t0)+" um\n")
        f.write("\tscanning frequency range for better accuracy is "+str(ring_mod.f_res_bar/ring_mod.Q)+" THz\n")
        f.write("\tCurrent frequency range is "+str(abs(ring_mod.f_end_bar-ring_mod.f_start_bar))+" THz\n")
        f.write("\tThe time for scanning : "+str(t.t_max)+" ps \n")
        f.write("\t\tNote. This time value may affect the accuracy of Transfer function. above 10000 ps is better.  \n")
        if  ring_mod.Q >ring_mod.f_res_bar/abs(ring_mod.f_end_bar-ring_mod.f_start_bar):
            f.write("\t\tWarning : This Transfer function result may not be accurate.\n\n")

    f.write("Simulation Parameters\n")
    f.write("\ttime scaling : "+str(utility.t0)+" sec\n")
    f.write("\ttime resolution : "+str(t.dt)+" sec\n")
    f.write("\tMax time range : "+str(t.t_max)+" ps\n")

