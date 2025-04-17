import numpy as np
from functools import partial
from scipy.integrate import solve_ivp

# 先定義好各種模式下的微分方程
def mode1(t, y, alpha, beta):
    print("alpha = ",alpha)
    print("beta = ",beta)
    """
    dy/dt = alpha * y - beta * y^2
    """
    return alpha * y - beta * y**2

def mode2(t, y, omega):
    """
    dy/dt = omega * sin(t)
    """
    return omega * np.sin(t)

def mode3(t, y):
    """
    dy/dt = -y
    """
    return -y

# 用 dictionary 對應不同模式
# 注意：若不同模式需要的參數不一樣，可以用 functools.partial 帶入
mode_dict = {
    'mode1': mode1,
    'mode2': mode2,
    'mode3': mode3,
}

def simulate_ode(mode, t_span, y0, **kwargs):
    """
    mode:   選擇要模擬的模式
    t_span: (t_start, t_end)
    y0:     初始條件 (list 或 np.array)
    kwargs: 給各個微分方程對應的參數
    """
    # 用 partial 把參數帶進對應的方程式
    ode_func = partial(mode_dict[mode], 1)

    # 執行 solve_ivp
    sol = solve_ivp(
        ode_func,             # 微分方程函式
        t_span,               # 時間範圍
        y0,                   # 初始條件
        t_eval=np.linspace(t_span[0], t_span[1], 100),  # 產生 100 個取樣點
        args=(1,2)
    )
    return sol

if __name__ == "__main__":
    # 模擬模式 1，並指定 alpha=1.0, beta=0.2
    sol1 = simulate_ode('mode1', (0, 10), [1.0], alpha=1.0, beta=0.2)
    # print("Mode1 solution y(t):", sol1.y)
    
    # 模擬模式 2，指定 omega=2.0
    sol2 = simulate_ode('mode2', (0, 10), [0.0], omega=2.0)
    # print("Mode2 solution y(t):", sol2.y)
    
    # 模擬模式 3，不用額外參數
    sol3 = simulate_ode('mode3', (0, 10), [1.0])
    # print("Mode3 solution y(t):", sol3.y)

    import matplotlib.pyplot as plt 
    plt.figure()
    plt.plot(sol1.y[0])
    plt.plot(sol2.y[0])
    plt.plot(sol3.y[0])
    plt.show()
