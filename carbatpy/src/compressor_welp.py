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

def set_up(T_inlet, p_inlet, p_outlet, fluid, comp, resolution):
    # initializing pZ vector
    pZ = np.zeros(7)
    pZ[0:6] = z_Tp(T_inlet, p_inlet, fluid,
                   comp)  # fl.zs_kg(['T','p'],[T_e,p_e],['T','p','v','u','h','s'],fluid) #state suction pipe
    pZ[6] = p_outlet  # pressure in pressure pipe
    # initializing pV vector
    pV = [34e-3, 34e-3, 3.5, .04, .06071, 48.916, 50., 50. / 2., 2.]  # parameter see above
    cycle_pos_var = np.linspace(0., 2 * np.pi, resolution)
    a_head = np.pi / 4. * pV[0] ** 2.  # area of cylinder head

def fun(x, u, pV, a_head, fluid, comp):
    pos_piston = -(pV[1] / 2. * (1. - np.cos(x) + pV[2] *
                    (1. - np.sqrt(1. - (1. / pV[2] * np.sin(x)) ** 2.)))) + pV[4] * pV[1] + pV[1]  # piston position, x=0 at UT
    volume_cylinder = a_head * pos_piston  # volume cylinder
    ht_surface = np.pi * pV[0] * pos_piston + 2. * a_head  # heat transfer surfaces
    vi = volume_cylinder / mi  # specific volume in cylinder, m³/kg
    zi = z_uv(u[1], vi, fluid, comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)
    return np.array([u[1],a * u[1]-x**2])

def bc(ya, yb):
    return np.array([yb[1] - ya[1], ya[0]])


x_var = np.linspace(0,10,100)
y_guess = np.ones([2,len(x_var)])
res = solve_bvp(fun, bc, x_var, y_guess)
plt.plot(x_var, res.y[1])
plt.show()