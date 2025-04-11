import cython_test
import time as timer
import numpy as np
import raise_cosine
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

N=1000.0
f = np.arange(N*N).reshape((int(N),int(N)))
g = np.arange(9.0).reshape((3,3))
print(type(f))
print(type(g))
t1 = timer.time()
ans = cython_test.naive_convolve(f,g)
t2 = timer.time()
print(type(ans))
print("Spent time = ",t2-t1)


def test(t):
    b = np.zeros((1,len(t)))
    for i in range(len(t)):
        b[0,i] = t[i]**0.5

    return b
t = np.arange(0,10000,0.0001)
t1 = timer.time()
b = cython_test.test(t)
# b = test(t)
t2 = timer.time()
print(type(t))
print(type(b))
print("Spent time = ",t2-t1)

