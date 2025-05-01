from driver import *
from time_class import *
from ring import *
from utility import ploting, dB

def FFT(time,driver):
    F = 1/len(time.t_total)*np.fft.fft(driver.v)
    F  = np.fft.fftshift(F)
    F_pos = 2*F[(len(F)//2)+1:(len(F)//2)+2000]
    f = np.fft.fftfreq(len(t.t_total),d = t.dt)
    f = np.fft.fftshift(f)
    f = f[(len(F)//2)+1:(len(F)//2)+2000]/1e9
    return f,F_pos

# Some Dummy components
# //////////////////////////////////////////////////////////////////////////////////////////////////////////
experiment_condition ={"mode":"voltage_drive",
                        "lambda_incident":1.55,
                        "Pin":1} 
sim = simulation()
sim.main(experiment_condition=experiment_condition)

ring_mod = ring(2*np.pi*5,
                0.99,
                29.5,
                0.5*0.22,
                1.55,
                0.99,
                lambda0=1.55,
                FSR = 0.019
                )
# ///////////////////////////////////////////////////////////////////////////////////////////////

v = driver(v_bias=-1,
           vpp=1,
           f_drive=50,
           R = 53.9,
           cj = [23.6*1e-15,20*1e-15],
           raise_cosine=1,
           PRBS=1,
           level='PAM4',)

t = time(mode = sim.mode)
t.main(ring_mod,v,N=1000)
sim.eye_diagram(t,v,signal=v.v,plot_bit_num=3,filename="pam4_test")
f_pam4, pam4_FT = FFT(t,v)

v.level = "NRZ"
v.renew()
t.main(ring_mod,v,N=1000)
sim.eye_diagram(t,v,signal=v.v,plot_bit_num=3,filename="nrz_test")

f, nrz_FT = FFT(t,v)
# ploting(time.t_total, driver.v,x_label="time",title="PAM4 signal")
# ploting(time.t_total, f,x_label="time",title="PAM4 signal")
ploting(f, dB(abs(nrz_FT)), dB(abs(pam4_FT)),x_label="frequency",title="FT of PAM4 and NRZ")
