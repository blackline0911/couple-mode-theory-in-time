from scipy.integrate import solve_ivp
import numpy as np

def simple_diff(t, z) :
    x, y = z
    return [1 - 2*x*y, 2*x - y]

t_range = (0, 500*3)
xy_init = [0, 0]
t_eval = np.arange(500*3+1)

sol = solve_ivp(simple_diff, t_range, xy_init,t_eval=t_eval)

print(len(t_eval))
print(len(sol.y[0]))
print(len(sol.y[1]))
import sys

# 打開文件，設置寫入模式
with open("output.txt", "w") as f:
    # 重定向標準輸出到文件
    sys.stdout = f

    # 這些 print() 的內容會寫入 output.txt
    print("Hello, World!")
    print("Python is amazing!")
    print(42)

# 恢復標準輸出
sys.stdout = sys.__stdout__

# 確認輸出到終端
print("輸出已寫入 output.txt！")

