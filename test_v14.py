import numpy as np
from scipy.integrate import solve_ivp

# 先定義好各種模式下的微分方程
def mode1(t, y, alpha, beta):
    print("alpha = ",alpha)
    print("beta = ",beta)
    """
    dy/dt = alpha * y - beta * y^2
    """
    return alpha * y - beta * y**2

sol = solve_ivp(
        mode1,             # 微分方程函式
        (0, 10),               # 時間範圍
        [1.0],                   # 初始條件
        t_eval=np.linspace(0, 10, 10),  # 產生 100 個取樣點
        args=(1,2)
    )