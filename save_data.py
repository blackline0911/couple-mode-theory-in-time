"""
This script is used to save the simulation variable that A simulation used. 
Just for Future checking or further research
"""

from test import *
import numpy as np
import ring, driver, time_class, cmt_solver 

filename = 'sim.txt'

with open(filename, 'w') as f:
    f.write('This script is used to save the simulation variable that A simulation used. \nJust for Future checking or further research\n\n')

    f.write("Parameter of ring modulator\n")
    f.write("\tStructure parameter : \n")
    f.write("\t\tradius : ")
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
    f.write('\n')

    m = t.mode
    f.write("Current Simulation mode is "+t.mode+'\n\n')
    if m == "voltage_drive":
        f.write("\tDriver information : \n")
        if v.sine_wave:
            f.write("\tDriving Signal is sine wave\n")
            f.write("\t\tDriving frequency is "+str(v.f_drive/1e9)+" GHz\n")
            f.write("\t\tBias voltage is "+str(v.v_bias)+" V")
            f.write("\t\tPeak to Peak voltage is "+str(v.vpp)+" V")