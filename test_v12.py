import numpy as np
import matplotlib.pyplot as plt

# t = np.arange(0,1,0.01)
# voltage = np.exp(1j*2*np.pi*1*t)
# output = np.cos(2*np.pi*1*t)
# # print(output/(voltage))
# plt.plot(t,voltage)
# plt.show()
# plt.plot(t,output)
# plt.show()
# plt.plot(t,output/(voltage))
# plt.show()
# print(abs(np.sum(output/(voltage)/1*0.01)))
# print((np.sum(output/(voltage)/1*0.01)))
# # plt.show()

import time
t1 = time.time()
x1 = np.arange(0,1000000,0.0001)
t2 = time.time()
print("x1 spent time= ",t2-t1)
t1 = time.time()
x2 = np.arange(0,1000,0.0001)
t2 = time.time()
print("x2 spent time= ",t2-t1)