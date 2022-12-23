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

    # setting of Aeff_i, explicit function
    M = z_mm(300, 100., fluid, comp)[
        -1]  # CP.PropsSI("M",fluid) # molar mass kg/mol # AW ramdom inlet conditions, molar mass constant
    Aeff_i = 2.0415e-3 * (Rm / M) ** (-.9826) * pV[0] ** 2. / Ver0[
        0] ** 2.  # effective flow cross-section inlet, m²
    # setting of Aeff_o, implicit function relatively to average mass flow density over valve
    # at 1st iteration, the mass flow density is unknown, typical value is guessed
    Aeff_o = 1.5e-5 * pV[0] ** 2. / Ver0[0] ** 2.
    # print(pZyk)


def fun(x, y, pV, a_head, fluid, comp, pZ):
    pos_piston = -(pV[1] / 2. * (1. - np.cos(x) + pV[2] *
                    (1. - np.sqrt(1. - (1. / pV[2] * np.sin(x)) ** 2.)))) + pV[4] * pV[1] + pV[1]  # piston position, x=0 at UT
    volume_cylinder = a_head * pos_piston  # volume cylinder
    ht_surface = np.pi * pV[0] * pos_piston + 2. * a_head  # heat transfer surfaces
    vi = volume_cylinder / y[0]  # specific volume in cylinder, m³/kg
    [Ti, pi, vi, ui, hi, si] = z_uv(y[1], vi, fluid, comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)
    T_thermal_mass =
    dxdt = -pV[1] / 2 * (np.sin(x) - pV[2] * (0.5 * (1 - (1 / pV[2] * np.sin(x)) ** 2))**-0.5 * \
           (-2 * 1 / pV[2] * np.sin(x)) * 1 / pV[2] * np.cos(x))
    dVdt = -np.pi ** 2 * f * pV[0] ** 2 / 2 * dxdt
    dW_fric = - pV[5] * dVdt
    dW_rev = -pi * dVdt
    if x <= np.pi:
        if pi <= pZ[6]:
            [alp, m_dot_in, m_dot_out] = compression(pV, pos_piston, dxdt, Ti, pi)
        else:
            z_it = push_out(i, fluid, z_it, comp, pV, pZyk, pZ)
    else:
        if pi >= pZ[1]:
            z_it = expansion(i, fluid, z_it, comp, pV)
        else:
            z_it = suction(i, fluid, z_it, comp, pV, pZyk, pZ)
    dQ = alp * ht_surface * (y[2] - Ti)  # kW
    dthermal_dt = state_th_Masse(y, -Q, pV)
    dmdt = m_dot_in - m_dot_out
    dudt = (dQ + dW_fric + dW_rev - dmdt * y[1] + m_dot_out * hi - m_dot_in * pZ[4]) / y[0]  # kJ/kg
    return np.array([dmdt,dudt, dthermal_dt])

def getalp(pV, step, pos_piston, dxdt, Ti, pi):
    '''
    calculates heat transfer coefficient gas/cylinder wall
    Woschni correlation
    '''
    if step == 0 or step == 2:  # closed valves
        k = 2.28
    else:  # open valves, suction or push out
        k = 5.18
    alp = 127.93 * pV[0] ** (-.2) * (pi * 1e-2) ** .8 * \
          (Ti) ** (-.55) * (k * dxdt) ** .8
    return alp

def state_th_Masse(y, Q, pV):
    '''
    calculates temperature change of thermal mass as function of heat transfer inside (Q)
    and to environment (Q_u)
    '''
    ### mass and cv of thermal mass are in stationary state not crucial,
    ### parameter are chosen to achieve fast convergence without vibrations
    m = .0001  # kg
    cv = .502  # kJ/kg/K
    alp_a = 6.  # heat transfer coefficient to environment
    A = Ver0[3] * pV[8] / Ver0[2] * pV[0] / Ver0[0] * pV[1] / Ver0[
        1]  # Outer surface cylinder estimated via geometry related to fitting compressor
    Q_u = alp_a * A * (Tu - y[2]) # kW
    dthermal_dt = (Q + Q_u) / cv / y[0]
    return dthermal_dt

def compression(pV, pos_piston, dxdt, Ti, pi):
    step = 0
    alp = (pV, step, pos_piston, dxdt, Ti, pi)
    m_dot_in = 0.  # no mass flow over boundaries
    m_dot_out = 0.
    return alp, m_dot_in, m_dot_out


def push_out(i, fluid, z_it, comp, pV, pZyk, pZ):
    step = 1
    alp = getalp(pV, step, pos_piston, dxdt, Ti, pi)
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i - 1, 12] - z_it[i - 1, 5]) * \
        ((z_it[i, 0] - z_it[i - 1, 0]) / (2. * np.pi * pV[7])) * 1e-3  # kJ
    z_it = state_th_Masse(-Q, z_it, i, pV)
    m_dot = pZyk[1] / z_it[i - 1, 7] * np.sqrt(2. * (z_it[i - 1, 6] - pZ[6]) * \
                                               1000. * z_it[i - 1, 7])  # mass flow leaving the cylinder, kg/s
    dm = m_dot * ((z_it[i, 0] - z_it[i - 1, 0]) / (2. * np.pi * pV[7]))
    # pushed out mass, kg
    mi = z_it[i - 1, 11] - dm  # resulting mass in cylinder, kg
    vi = z_it[i, 2] / mi  # resulting specific volume in cylinder, m³/kg
    ui = (Q + W + Wr - dm * z_it[i - 1, 9] + z_it[i - 1, 11] * z_it[i - 1, 8]) / mi
    # energy balance, resulting internal energy in cylinder, kj/kg
    # property call for state i with u and v
    zi = z_uv(ui, vi, fluid, comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)

    # store of quantities in z_it
    z_it[i, 4] = step
    z_it[i, 5:11] = zi
    z_it[i, 11] = mi
    z_it[i, 14:16] = dm, Q
    return z_it


def expansion(i, fluid, z_it, comp, pV):
    step = 2
    alp = getalp(pV, step, pos_piston, dxdt, Ti, pi)
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i - 1, 12] - z_it[i - 1, 5]) * \
        ((z_it[i, 0] - z_it[i - 1, 0]) / (2. * np.pi * pV[7])) * 1e-3  # kJ

    state_th_Masse(-Q, z_it, i, pV)
    dm = 0.  # no mass flow over boundaries
    mi = z_it[i - 1, 11]  # mass in cylinder, kg
    ui = (Q + W + Wr) / z_it[i - 1, 11] + z_it[i - 1, 8]  # energy balance
    vi = z_it[i, 2] / mi  # resulting specific volume in cylinder, m³/kg
    # property call for state i with u and v
    zi = z_uv(ui, vi, fluid, comp)  # zi=fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)

    # store of quantities in z_it
    z_it[i, 4] = step
    z_it[i, 5:11] = zi
    z_it[i, 11] = mi
    z_it[i, 14:16] = dm, Q
    return z_it


def suction(i, fluid, z_it, comp, pV, pZyk, pZ):
    step = 3
    alp = getalp(pV, step, pos_piston, dxdt, Ti, pi)
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i - 1, 12] - z_it[i - 1, 5]) * \
        ((z_it[i, 0] - z_it[i - 1, 0]) / (2. * np.pi * pV[7])) * 1e-3  # kJ
    state_th_Masse(-Q, z_it, i, pV)

    m_dot = pZyk[0] / pZ[2] * np.sqrt(
        2. * (pZ[1] - z_it[i - 1, 6]) * 1000 * pZ[2])  # mass flow leaving cylinder, kg
    dm = m_dot * ((z_it[i, 0] - z_it[i - 1, 0]) / (2. * np.pi * pV[7]))  # pushed out mass, kg
    mi = z_it[i - 1, 11] + dm  # resulting mass in cylinder, kg
    vi = z_it[i, 2] / mi  # resulting specific volume in cylinder, m³/kg
    ui = (Q + W + Wr + dm * pZ[4] + z_it[i - 1, 11] * z_it[i - 1, 8]) / mi  # isenthalpic throttling to cylinder
                                                                            # pressure
    # property call for state i with u and v
    zi = z_uv(ui, vi, fluid, comp)  # zi=fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)
    # store of quantities in z_it
    z_it[i, 4] = step
    z_it[i, 5: 11] = zi
    z_it[i, 11] = mi
    z_it[i, 14:16] = dm, Q
    return z_it

def bc(ya, yb):
    return np.array([yb[1] - ya[1], ya[0]])


x_var = np.linspace(0,10,100)
y_guess = np.ones([2,len(x_var)])
res = solve_bvp(fun, bc, x_var, y_guess)
plt.plot(x_var, res.y[1])
plt.show()