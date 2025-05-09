import matplotlib.pyplot as plt
import numpy as np

# 獲取待擬合數據
x = np.linspace(1, 50, 50)
f = np.poly1d([2, 5, 10])
y = f(x)

# 拟合 返回值就是多项式的系数 从左到右对应次数从高到低
# deg指的是多项式的最高次数
param = np.polyfit(x, y, deg=2)
print(param)

# 利用拟合得到系数计算函数数值
# 也可以使用以下代码计算函数数值：f = np.poly1d(param) y = f(x)
z = np.polyval(param, x)
# plt.plot(x, y, marker='o')
plt.plot(x, param[0]*x**2+param[1]*x+param[2])
plt.plot(x, z)

plt.show()


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

xdata = np.linspace(0, 4, 50)
y = func(xdata, 2.5, 1.3, 0.5)
rng = np.random.default_rng()
y_noise = 0.2 * rng.normal(size=xdata.size)
ydata = y + y_noise
plt.plot(xdata, ydata, 'b-', label='data')

popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 1., 0.5]))
print(popt) # fitting parameters
print(pcov) # error
plt.plot(xdata, func(xdata, *popt), 'g--',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
