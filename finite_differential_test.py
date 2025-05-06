from driver import *
from time_class import *
from ring import *
from utility import ploting, dB

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

# t = time(mode = sim.mode)
# t.main(ring_mod,v,N=100)
# sim.eye_diagram(t,v,signal=v.v,plot_bit_num=3)

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Testing Finite_differential Method

def FDM(time_array,signal,accuracy_level = 2):
    N = len(signal)
    dt = time_array[1]-time_array[0]
    print("dt in GDM function = ",dt)
    FDM = np.zeros(N)
    if accuracy_level==2:
        FDM[0] = (-3*signal[0] + 4*signal[1] - signal[2])/2/dt
        FDM[-1] = (3*signal[-1] - 4*signal[-2] + signal[-3])/2/dt
        for i in range(1,N-1):
            FDM[i] = (signal[i+1]-signal[i-1])/2/dt

    return FDM


import time 
dt = 0.00001
t = np.arange(0,10,dt)
y = np.cos(2*np.pi/2*t)
ploting(t,y,x_label="time",title='original signal')
dy_dt = -2*np.pi/2*np.sin(2*np.pi/2*t)
ploting(t,dy_dt,x_label="time",title='analytical signal')
t1 = time.time()
test_1 = FDM(t,y)
t2 = time.time()
t3 = time.time()
test_2 = np.gradient(y,dt)
t4 = time.time()
ploting(t,test_1,test_2,x_label="time",title='analytical signal',leg=['My code','np.gradient'])
ploting(t,dy_dt-test_2,x_label="time",title='Error of numpy.gradient')
ploting(t,dy_dt-test_1,x_label="time",title='Error of iterating FDM')
print("time for iterating FDM = ",t2-t1," second")
print("time for numpy.gradient = ",t4-t3," second")


# Result：
# Error of numpy.gradient is about 1e-7, and 5e-4 at boundary (for dt = 1e-3)
# While Error of Mycode is about 5e-6 (for dt = 1e-3)
# time for iterating FDM =  0.4827697277069092  second (for dt = 1e-5)
# time for numpy.gradient =  0.0  second (for dt = 1e-5)
