import numpy as np
from utility import *
from cmath import *
import matplotlib.pyplot as plt
class simulation():
    voltage_drive = False
    scan_frequency = False
    mode = None
    discarding = 2
    method = "small_signal"
    mode_dict = {}
    id = "simulation"
    
    # Result Variables
    b = np.array([])
    N = np.array([])

    # method of solving algorithm
    algorithm = 'RK45'
    
    def __init__(self):
        pass
    def main(self,experiment_condition):
        """
        Initialize the simulation.
        mode only support "voltage_drive", "scan_frequency"
        experiment_condition : a dictionary that storing simulation mode 
                                , input laser wavelength and Power (mW)
                                eg. experiment_condition =
                                    {"mode":"scan_frequency",
                                    "lambda_incident":1.55,
                                    "Pin":1} 
        """
        modes = {"voltage_drive", "scan_frequency"}
        assert experiment_condition["mode"] in modes, f"Please choose a simulation mode from {modes}"
        match experiment_condition["mode"]:
            case "voltage_drive":
                self.voltage_drive = True
            case "scan_frequency":
                self.scan_frequency = True
        self.mode = experiment_condition["mode"]
        self.lambda_incident = experiment_condition["lambda_incident"]
        try:
            self.method = experiment_condition["method"]
        except:
            pass
         # normalized input laser frequency
        self.f_pround_bar = c/(self.lambda_incident)*t0
        self.Pin = experiment_condition["Pin"]
        self.S0 = sqrt(self.Pin)
        self.b0 = np.real(sqrt(t0)*self.S0)
    def renew(self):
        self.f_pround_bar = c/(self.lambda_incident)*t0
        self.S0 = sqrt(self.Pin)
        self.b0 = np.real(sqrt(t0)*self.S0)
        
    def set_dt(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")

    def create_time_array(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")
    
    
    
    def eye_diagram(self,time,driver,signal,title = "Output Power",filename='',plot_bit_num=2,
                    linewidth=0.15):
        assert self.mode=="voltage_drive", "\neye diagram only available at voltage driving mode\n"
        N = time.N
        # np.save("./"+filename+".npy",signal)
        cum_t_index, sig= self.discard_func(time,signal)
        step = int(cum_t_index[1]-cum_t_index[0])
        
        plt.figure()
        # discarding the first two bits
        if driver.PRBS:
            t_array = np.array(time.t_all_segment[0][:])
            for j in range(1,plot_bit_num):
                t_array = np.append(t_array,np.array(time.t_all_segment[j][:]))
            total_seg = len(sig)//(step*plot_bit_num)
            sig = sig[0:step*total_seg*plot_bit_num]
            segment = np.array_split(sig,total_seg)
            for s in (segment):
                plt.plot( t_array,s, color='blue', linewidth=linewidth)
            # for k in range(0,N-self.discarding-plot_bit_num-driver.shift_bit+1,plot_bit_num):
            #     print(k)
                
            #     sig_segment = sig[int(cum_t_index[k]-step*self.discarding):int(cum_t_index[k+plot_bit_num]-step*self.discarding)]
            #     sig_segment = sig[int(cum_t_index[k]-step*self.discarding):int(cum_t_index[k+plot_bit_num]-step*self.discarding)]
                
        else:
            t_array = np.array(time.t_all_segment[0][:])
            for j in range(1,plot_bit_num):
                    t_array = np.append(t_array,np.array(time.t_all_segment[j][:]))
            for k in range(N-self.discarding-plot_bit_num,plot_bit_num):
                sig_segment = sig[int(cum_t_index[k]-step*self.discarding):int(cum_t_index[k+plot_bit_num]-step*self.discarding)]
                plt.plot( t_array,sig_segment, color='blue', linewidth=linewidth)

        plt.grid(color='w')
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=7)
        plt.xlabel("time (second)")
        plt.title("Eye diagram of "+title)
        plt.grid(color='g',linestyle='--', alpha=0.5)
        fig = plt.gcf()
        if (not filename==''):
            plt.savefig(filename)
        
        plt.show()
    def save_data(self,*obj,file_name = 'sim.txt'):
        with open(file_name, 'w') as f:
            f.write('This script is used to save the simulation variable that A simulation used. \nJust for Future checking or further research\n\n')
            f.write("Simulation Parameters and Settings\n")
            f.write("\ttime scaling : "+str(t0)+" sec\n")
            f.write("\tSimulation mode = "+str(self.mode)+"\n")
            f.write("\tInput Laser Power : "+str(self.Pin)+" mW\n\n")
            f.write("\tInput Laser wavelength : "+str(self.lambda_incident)+" um\n\n")
            for ob in obj:
                match ob.id:
                    case "time":
                        f.write("Simulation Parameter of time \n")
                        f.write("\ttime resolution : "+str(ob.dt)+" sec\n")
                        f.write("\tMax time range : "+str(ob.t_max)+" ps\n")
                        f.write("\tNumber of bit in squence : "+str(ob.N)+"\n")
                        f.write("\n\tNote. In ScanFrequency mode ,the Max time range value may affect the accuracy of Transfer function. above 10000 or 5000 ps is better. (This value may change by case ) \n")
                    case "driver":
                        if self.mode == "voltage_drive":
                            f.write("\tDriver information : \n")

                            if ob.sine_wave:
                                f.write("\tDriving Signal is sine wave\n")
                                f.write("\t\tDriving frequency is "+str(ob.f_drive/1e9)+" GHz\n")
                                f.write("\t\tBias voltage is "+str(ob.v_bias)+" V\n")
                                f.write("\t\tPeak to Peak voltage is "+str(ob.vpp)+" V\n\n")
                            if ob.square_wave:
                                f.write("\tDriving Signal is square wave\n")
                                f.write("\t\tDriving frequency is "+str(ob.f_drive/1e9)+" GHz\n")
                                f.write("\t\tBias voltage is "+str(ob.v_bias)+" V\n")
                                f.write("\t\tPeak to Peak voltage is "+str(ob.vpp)+" V\n\n")
                            if ob.raise_cosine:
                                f.write("\tDriving Signal is raise_cosine signal\n")
                                f.write("\t\tDriving frequency is "+str(ob.f_drive/1e9)+" GHz\n")
                                f.write("\t\tBias voltage is "+str(ob.v_bias)+" V\n")
                                f.write("\t\tPeak to Peak voltage is "+str(ob.vpp)+" V\n\n")
                                f.write("\t\tbeta of raise cosine signal is "+str(ob.beta)+"\n\n")
                        
                    case "ring":
                        f.write("Parameter of ring modulator\n")
                        f.write("\tStructure parameter : \n")
                        f.write("\t\tcavity length : ")
                        f.write(str(ob.L)+' um')
                        f.write('\n')

                        if ( ob.lambda0==None):
                            f.write("\t\tneff : ")
                            f.write(str(ob.neff(0)))
                            f.write('\n')

                        f.write("\t\tGroup index : ")
                        f.write(str(ob.ng))
                        f.write('\n')

                        f.write("\t\tcouple through coefficient : ")
                        f.write(str(ob.gamma))
                        f.write('\n')

                        f.write("\t\tRound trip loss : (energy) ")
                        # f.write(str(np.real(2*ob.alpha(0))))
                        f.write(str(np.exp(-ob.L*1e-4*ob.alpha(0))))
                        f.write('\n')

                        # f.write("\t\tmodulation efficiency : ")
                        # f.write(str(ob.me) + ' pm/V (for reverse bias)')
                        # f.write('\n')

                        f.write("\t\tmode cross section : ")
                        f.write(str(ob.cross_section) + ' um^2')
                        f.write('\n')

                        f.write("\t\tLinear Energy absorption coefficient: ")
                        f.write(str(ob.alpha_linear)+" 1/cm")
                        f.write('\n')

                        f.write("\t\tEnergy coupling coefficient: ")
                        f.write(str(-2/(ob.L*1e-4)*np.log(ob.gamma))+" 1/cm")
                        f.write('\n')

                        f.write("\t\tMax Two Photon Absorption coefficient : ")
                        f.write(str(ob.TPA_coeff*(max(abs(self.b)**2)))+" 1/cm")
                        f.write('\n')
    
                        f.write("\t\tMax Free Carrier Absorption coefficient : ")
                        f.write(str(ob.FCA_coeff*(max(abs(self.b)**4)))+" 1/cm")
                        f.write('\n')

                        f.write("\t\tFree Carrier time constant : ")
                        f.write(str(ob.tau_eff)+" ns")
                        f.write('\n')

                        f.write("\t\tTwo photon absorption coefficient : ")
                        f.write(str(ob.beta_TPA*1e-13)+" cm/mW")
                        f.write('\n')

                        f.write("\t\tTwo photon absorption mode area : ")
                        f.write(str(ob.Atpa)+" um^2")
                        f.write('\n')


                        f.write("\t\tFree carrier absorption area : ")
                        f.write(str(ob.sigma_FCA*1e-17)+" cm^2")
                        f.write('\n')

                        f.write("\t\tFree carrier absorption mode area : ")
                        f.write(str(ob.Afca)+" um^2")
                        f.write('\n')

                        f.write("\t\tHeater efficiency : ")
                        f.write(str(ob.HE)+" pm/mW")
                        f.write('\n')

                        f.write("\t\tKerr effect coefficient (n2) : ")
                        f.write(str(ob.n2)+" um^2/mW")
                        f.write('\n')

                        f.write("\t\tKerr effect mode area (Akerr) : ")
                        f.write(str(ob.Akerr)+" um^2")
                        f.write('\n')

                        f.write("\t\tindex change rate per Kelvin : ")
                        f.write(str(ob.kappa_thermal)+" 1/K")
                        f.write('\n')


                        f.write("\t\tDensity of Silicon : ")
                        f.write(str(ob.ro_si)+" g/um^3")
                        f.write('\n')

                        f.write("\t\tspecific heat capacity of silicon : ")
                        f.write(str(ob.cSi)+" mJ/(g*K)")
                        f.write('\n')

                        f.write("\tPhysical parameter : \n")
                        f.write("\t\tintrinsic photon life time : ")
                        f.write(str(np.real(ob.tu_o_bar))+' ps')
                        f.write('\n')
                        for i in range(len(ob.tu_e_bar)):
                            f.write("\t\tcouple extrinsic photon life time (port "+str(i+1)+" ) : ")
                            f.write(str(np.real(ob.tu_e_bar[i]))+' ps')
                            f.write('\n')

                        f.write("\t\tresonant frequency : ")
                        f.write(str(ob.f_res_bar)+' THz')
                        f.write('\n')

                        f.write("\t\tresonant wavelength : ")
                        f.write(str(ob.lambda0)+' um')
                        f.write('\n')

                        f.write("\t\tQ factor : ")
                        f.write(str(ob.Q))
                        f.write('\n')

                        f.write("\t\tLineWidth : ")
                        f.write(str(c*(1e-12/1e-3)/ob.f_res_bar**2*ob.f_res_bar/ob.Q)+" nm")
                        f.write('\n')

                        f.write("\t\tGroup velocity : ")
                        f.write(str(ob.vg/t0)+" um/s")
                        f.write('\n')

                        f.write("\t\tRound Trip time : ")
                        f.write(str(ob.round_trip_time)+" ps")
                        f.write('\n\n')

                        f.write("\t\tInput Laser photon energy : ")
                        f.write(str(ob.photon_energy)+" fJ")
                        f.write('\n\n')

                        if self.mode == "scan_frequency":
                            f.write("\tThe start wavelength is "+str(c/ob.f_start_bar*t0)+" um\n")
                            f.write("\tThe end wavelength is "+str(c/ob.f_end_bar*t0)+" um\n")
                            f.write("\tscanning frequency range for better accuracy is "+str(ob.f_res_bar/ob.Q)+" THz\n")
                            f.write("\tCurrent frequency range is "+str(abs(ob.f_end_bar-ob.f_start_bar))+" THz\n")
                            
                            if  ob.Q >ob.f_res_bar/abs(ob.f_end_bar-ob.f_start_bar):
                                f.write("\t\tWarning : This Transfer function result may not be accurate.\n\n")

                    case "Heater":
                        f.write("\n\t\tVoltage applied to Heater is "+str(ob.V)+" V\n")
                        f.write("\t\tResistance of Heater is "+str(ob.R)+" ohm\n")
                        f.write("\t\tThermal time constant of waveguide is "+str(ob.tau_th*1e9)+" ns\n")
                        f.write("\t\tThermal mode area : ")
                        f.write(str(ob.Aeff)+" um^2")
                        f.write('\n')

    def discard_func(self,time,signal):
        N = time.N
        t = time.t_total
        assert N>self.discarding, "\nAvailable Bit number is not enough\n"
        cum_t_index = np.zeros(N)
        seg = time.t_all_segment[0]
        sig = signal
        cum_t_index[0] = (len(time.t_all_segment[0]))
        for q in range(1,N):
            cum_t_index[q] = (len(time.t_all_segment[q]) + cum_t_index[q-1])
        return cum_t_index[self.discarding-1::], sig[int(cum_t_index[self.discarding-1])::]

