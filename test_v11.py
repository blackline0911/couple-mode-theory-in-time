import cython_test
import time as timer
import numpy as np
def factorial(n):
    ans = 1
    for i in range(1,n+1):
        ans*=i
    return ans

t1 = timer.time()
b=0
for i in range(100):
    b+= (cython_test.factorial(5,1))
t2 = timer.time()
print("time spent run by cython = ",t2-t1)

# t1 = timer.time()
# a=0
# for i in range(1000000000):
#     a+=(factorial(5))
# t2 = timer.time()
# print("time spent run by python = ",t2-t1)

"""
time spent run by cython =  137.29582738876343
time spent run by python =  374.89023184776306
"""

N=10.0
f = np.arange(N*N).reshape((int(N),int(N)))
g = np.arange(9.0).reshape((3,3))

print(cython_test.naive_convolve(f,g))