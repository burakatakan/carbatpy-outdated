"""
compressor model adapted from Dennis Roskosch
changed to use ODE solver instead of finite volume approach

author: Alexandra Welp
21.12.2022
"""

import numpy as np
import matplotlib.pyplot as plt
from carbatpy.fl_props_compressor import z_uv, z_ps, z_Tp, z_Tx, z_mm
from scipy.integrate import solve_bvp

Rm = 8.3145  # gas constant J/mol/K
Tu = 25. + 273.15  # ambient temperature
dTK = 273.15  # conversion °C / K
Ver0 = [34e-3, 34e-3, 2., .04]  # fit-compressor: D, H, cylinder, outer surface

def fun(x,u):
    a = 2
    return np.array([u[1],a * u[1]-x**2])

def bc(ya, yb):
    return np.array([yb[1] - ya[1], ya[0]])


x_var = np.linspace(0,10,100)
y_guess = np.ones([2,len(x_var)])
res = solve_bvp(fun, bc, x_var, y_guess)
plt.plot(x_var, res.y[1])
plt.show()