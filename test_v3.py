import numpy as np
import time 

"""
This file is testing which method handling function for array input is faster.
That can be applied to method of faster way of dealing the voltage array generation, which may cost lots of time iterating by for loop.  
"""
def my_function_v1(x):
    if x > 0:
        return x**2 + np.sin(x)
    else:
        return np.cos(x) - x**3

def my_function_v2(x):
    if isinstance(x, np.ndarray):
        return np.where(x > 0, x**2 + np.sin(x), np.cos(x) - x**3)  # NumPy 陣列操作
    else:
        return x**2 + np.sin(x) if x > 0 else np.cos(x) - x**3  # 單個 float 處理

# 測試

a = time.time()
x = np.linspace(-10,10,100000)
for i in range(len(x)):
    my_function_v1(x[i])
b = time.time()
print("Operation time of iterating by for loop is : ",b-a)

a = time.time()
vectorized_function = np.vectorize(my_function_v1)
vectorized_function(x)
b = time.time()
print("Operation time of using vectorizaed function is : ",b-a)

a = time.time()
my_function_v2(x)
b = time.time()
print("Operation time of using isinstance() function is : ",b-a)


"""
Result:
Operation time of iterating by for loop is :  0.13918566703796387
Operation time of using vectorizaed function is :  0.13659191131591797
Operation time of using isinstance() function is :  0.006720542907714844

=> Using isinstance() is 20 times faster than others , since this method does not involve for loop in it.
"""

def sinc(t):
        if isinstance(t, np.ndarray):
            # return np.where(t==0.0,1,np.sin(np.pi*t)/(np.pi*t))
            a = np.where(np.isclose(t,0.0),1,np.sin(np.pi*t)/(np.pi*t))
            return a
        else:
            if t==0.0:
                return 1
            else:
                return np.sin(np.pi*t)/(np.pi*t)
x = np.array([-1,0,1])
print(x)
print(np.where(np.isclose(x,0.0),1,np.sin(np.pi*x)/(np.pi*x)))
print(sinc(x))