from scipy.integrate import *
import numpy as np
# Define a function for the right side
def dU_dx_new(x, U):
    y , y_pround = U
    """Right side of the differential equation to be solved.
    U is a two-component vector with y=U[0] and z=U[1]. 
    Thus this function should return [y', z']
    """
    dy_dx = y_pround

    dyp_dx = -2*dy_dx - 50*y
    return [dy_dx, dyp_dx] 
    # return [U[1], -2*U[1] - 2*U[0] + np.cos(2*x)]

# initial condition U_0 = [y(0)=0, z(0)=y'(0)=0]
U_0 = [2., 0.]

x_pts = np.linspace(0, 15, 2000)  # Set up the mesh of x points
result = solve_ivp(dU_dx_new, (0, 15), U_0, t_eval=x_pts)
y_pts = result.y[0,:]   # Ok, this is tricky.  For each x, result.y has two 
                        #  components.  We want the first component for all
                        #  x, which is y(x).  The 0 means the first index and 
                        #  the : means all of the x values.

import matplotlib.pyplot as plt
plt.grid(color='g',linestyle='--', alpha=0.5)
plt.plot(x_pts, result.y[0])
plt.plot(x_pts, result.y[1])
plt.show()