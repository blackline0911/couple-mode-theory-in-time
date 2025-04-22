from ring import * 
from driver import *
from time_class import *
from cmt_solver import *
from utility import *


# Refer to paper : "O-Band Silicon Ring Modulators With Highly
# Efficient Electro- and Thermo-Optic Modulation"
# Low Q factor(5000) and Lateral PN
wl_in = 1.3055
Pin = 1
mode = "scan_frequency"
radius = 5
FSR = 0.0135
wl_res = 1.30552
me = 29.4
Rs=68.4
cj = [23.8e-15, 20.1e-15]
a = 0.96616
gamma = 0.97115
v_bias = -1
vpp = 1

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

exp_condition = {
    'mode':mode,
    'Pin' :Pin,
    'lambda_incident':wl_in
}
sim = simulation()
sim.main(experiment_condition=exp_condition)

ring_mod = ring(2*np.pi*radius,
                a,
                me,
                0.22*0.37,
                wl_in,
                gamma,
                FSR=FSR,
                lambda0=wl_res,
                )

wl_min =  1.3052
wl_max =  1.3058 

v = driver(f_drive=50,
           v_bias=v_bias,
           vpp=vpp,
           R=53.9,
           raise_cosine=1,
           sine_wave=0,
           cj = [23.6e-15, 20e-15],
           PRBS=1)

os.chdir("./imec_oband_ring/")
t = time(mode = sim.mode)
if sim.mode == "scan_frequency":
    vbias = np.array([0.5,0,-0.5])
    ring_mod.scan_frequency(wl_min ,wl_max,t)
    t.main(ring_mod,t_max=5000,resolution=1,buffer=100,driver=v)
    wl_scan =  c/ring_mod.w_res(t.t_total)*t0
    sim.save_data(ring_mod,t,v)
    TP = np.zeros( int(len(t.t_total)-t.buffer*t0/t.dt-1))
    ER = np.zeros( int(len(t.t_total)-t.buffer*t0/t.dt-1))
    data_record = np.zeros((len(vbias),int (len(t.t_total)-t.buffer*t0/t.dt-1)))
    plt.figure()
    i = 0
    for vb in vbias:
        v.v_bias = vb
        b,Q,s_minus,N = solving(sim,ring_mod,v,t)
        T = Transfer_function(ring_mod,t)
        wl,data = T.mapping((abs(s_minus)**2/sim.Pin))
        wl,data_phase = T.mapping(180/np.pi*np.angle(s_minus))
        plt.plot(wl,10*np.log10(data),
                 label=str((vb)))
        data_record[i,:] = data
        i+=1
    TP = -10*np.log10( abs( data_record[2] - data_record[0] )/2 )
    for i in range(len(data)):
        if data_record[2,i]>data_record[0,i]:
            ER[i] = -10*np.log10(data_record[0,i]/data_record[2,i])
        else:
            ER[i] = -10*np.log10(data_record[2,i]/data_record[0,i])

    plt.grid(color='g',linestyle='--', alpha=0.5)
    plt.xlabel('wavelength(um)')
    plt.title('Transfer function')
    plt.legend()
    plt.savefig("Tranfer_function")
    plt.show()

    ploting(wl, data_record[0,:], data_record[2,:] ,x_label="wavelength (um)", title="Tranfer_functio")
    ploting(wl, TP,ER, x_label="wavelength (um)", title="TP and ER", filename="TP_and_ER",leg=['TP','ER'])

    # v.v_bias=-1.5
    # b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    # T = Transfer_function(ring_mod,t)
    # wl,data_v1 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    # wl,data_phase_v1 = T.mapping(180/np.pi*np.angle(s_minus))
    # v.v_bias=0
    # b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    # T = Transfer_function(ring_mod,t)
    # wl,data_v2 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    # wl,data_phase_v2 = T.mapping(180/np.pi*np.angle(s_minus))
    # v.v_bias=-0.5
    # b,Q,s_minus,N = solving(sim,ring_mod,v,t)
    # T = Transfer_function(ring_mod,t)
    # wl,data_v3 = T.mapping(10*np.log10(abs(s_minus)**2/sim.Pin))
    # wl,data_phase_v3 = T.mapping(180/np.pi*np.angle(s_minus))