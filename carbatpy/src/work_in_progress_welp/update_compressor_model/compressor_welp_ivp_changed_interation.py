"""
compressor model adapted from Dennis Roskosch
changed to use ODE solver instead of finite volume approach

author: Alexandra Welp
21.12.2022
"""

import numpy as np
import matplotlib.pyplot as plt
from carbatpy.fl_props_compressor import z_uv, z_ps, z_Tp, z_Tx, z_mm
from scipy.integrate import solve_ivp

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

    pZyk = np.zeros(2)
    pZyk[0] = Aeff_i
    pZyk[1] = Aeff_o
    resolution = 3600

    x_max = 1 / pV[7]
    y_start = np.zeros(3)
    y_start[0] = Ver0[1] * a_head / pZ[2]       # mass in cylinder
    y_start[1] = pZ[3]                          # u in cylinder
    y_start[2] = Tu                             # T of Thermal mass
    err = 10                                    # start value
    count = 0                                   # counts number of cycles until desired accuracy as achieved
    y_timetrack_m = []                          # tracks properties after every cycle
    y_timetrack_u = []
    y_timetrack_t = []
    y_timetrack_m.append(y_start[0])
    y_timetrack_u.append(y_start[1])
    y_timetrack_t.append(y_start[2])
    while err > 0.5:
        res = solve_ivp(fun, [0, x_max], y_start, method='RK45', args=[pV, a_head, fluid, comp, pZ, pZyk], max_step=1/resolution)
        err = np.sqrt((res.y[0, -1] - y_start[0]) ** 2) + np.sqrt((res.y[1, -1] - y_start[1]) ** 2) + np.sqrt((res.y[2, -1] - y_start[2]) ** 2)
        pZyk[1] = help_variable
        if pZyk[1] == Aeff_o:
            (f"Aeffout not updatet, count = {count}")
        y_timetrack_m.append(res.y[0,-1])
        y_timetrack_u.append(res.y[1, -1])
        y_timetrack_t.append(res.y[2, -1])
        y_start[0] = res.y[0, -1]
        y_start[1] = res.y[1, -1]
        y_start[2] = res.y[2, -1]
        count += 1
    print(f"Anzahl der Zyklen: {count}")

    iss_eff, degree_delivery = cal_efficiency_delivery(res, pV, pZ, pZyk, a_head, fluid, comp)

    #TODO implementieren von dhdt, damit h auch output ist


    return res, count, y_timetrack_m, y_timetrack_u, y_timetrack_t, iss_eff, degree_delivery

def cal_efficiency_delivery(res, pV, pZ, pZyk, a_head,fluid,comp ):
    # calculation of isentropic efficiency and delivery ratio
    m_max = np.array(res.y[0]).max()
    m_min = np.array(res.y[0]).min()
    m_aus = m_max - m_min  # overall pushed out mass
    m0 = np.pi * pV[0] ** 2. * pV[1] / pZ[2] / 4.  # sucked-in mass ideal compressor
    degree_delivery = m_aus / m0  # degree of delivery

    theta = res.t * (2 * np.pi * pV[7])
    pos_piston = -(pV[1] / 2. * (1. - np.cos(theta) + pV[2] *
                                 (1. - np.sqrt(1. - (1. / pV[2] * np.sin(theta)) ** 2.)))) + pV[4] * pV[1] + pV[
                     1]  # piston position, x=0 at UT
    volume_cylinder = a_head * pos_piston  # volume cylinder
    vi = volume_cylinder / res.y[0]  # specific volume in cylinder, m³/kg
    dxdtheta = -pV[1] / 2 * np.sin(theta) * (
                1 + 1 / pV[2] * np.cos(theta) * (1 - (1 / pV[2] * np.sin(theta)) ** 2) ** -0.5)
    dxdt = (2 * np.pi * pV[7]) * dxdtheta
    m_out = []
    h_out = []
    for i, t in enumerate(res.t):
        [Ti, pi, vi[i], ui, hi, si] = z_uv(res.y[1, i], vi[i], fluid,
                                    comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)
        if Ti == -9999990.:
            raise ValueError("invalid properties")
        if theta[i] <= np.pi and pi >= pZ[6]:
            [alp, m_dot_in, m_dot_out] = push_out(pV, pos_piston[i], dxdt[i], Ti, pi, pZ, pZyk, vi[i])
            mass_flow_density = m_dot_out / pZyk[1]
            pZyk[1] = 5.1109e-4 * mass_flow_density ** -.486 * pV[0] ** 2 / Ver0[0] ** 2  # Aeff_o neu
            array_parts = len(theta)
            h_out.append(m_dot_out * 1/pV[7] / array_parts * hi)
            m_out.append(m_dot_out* 1/pV[7] / array_parts)


    h_aus = np.sum(h_out)/ np.sum(m_out)  # average push out enthalpy
    h_aus_s = z_ps(pZ[6], pZ[5], fluid, comp)[
        4]  # fl.zs_kg(['p','s'],[pZ[6],pZ[5]],['h'],fluid)[0]  # isentropic outlet enthalpy
    is_eff = (h_aus_s - pZ[4]) / (h_aus - pZ[4])  # isentropic efficiency
    return is_eff, degree_delivery


def fun(x, y, pV, a_head, fluid, comp, pZ, pZyk):
    if y[0] < 0:
        y[0] = 0.0002

    theta = x * (2 * np.pi * pV[7])
    pos_piston = -(pV[1] / 2. * (1. - np.cos(theta) + pV[2] *
                    (1. - np.sqrt(1. - (1. / pV[2] * np.sin(theta)) ** 2.)))) + pV[4] * pV[1] + pV[1]  # piston position, x=0 at UT
    volume_cylinder = a_head * pos_piston  # volume cylinder
    ht_surface = np.pi * pV[0] * pos_piston + 2. * a_head  # heat transfer surfaces
    vi = volume_cylinder / y[0]  # specific volume in cylinder, m³/kg
    dxdtheta = -pV[1] / 2 * np.sin(theta) * (1 + 1/pV[2] * np.cos(theta) * (1 - (1/pV[2] * np.sin(theta))**2)**-0.5)
    dxdt = (2 * np.pi * pV[7]) * dxdtheta
    dVdt = a_head * dxdt
    [Ti, pi, vi, ui, hi, si] = z_uv(y[1], vi, fluid, comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)
    if Ti == -9999990.:
            raise ValueError("invalid properties")
    if theta <= np.pi:
        dW_fric = - pV[5] * dVdt
        if pi <= pZ[6]:
            [alp, m_dot_in, m_dot_out] = compression(pV, pos_piston, dxdt, Ti, pi)
        else:
            [alp, m_dot_in, m_dot_out] = push_out(pV, pos_piston, dxdt, Ti, pi, pZ, pZyk, vi)
            mass_flow_density = m_dot_out / pZyk[1]
            pZyk[1] = 5.1109e-4 * mass_flow_density ** -.486 * pV[0] ** 2 / Ver0[0] ** 2  # Aeff_o neu
            global help_variable
            help_variable = pZyk[1]
    else:
        dW_fric = pV[5] * dVdt
        if pi >= pZ[1]:
                [alp, m_dot_in, m_dot_out] = expansion(pV, pos_piston, dxdt, Ti, pi)
        else:
                [alp, m_dot_in, m_dot_out] = suction(pV, pos_piston, dxdt, Ti, pi, pZ, pZyk)



    dW_rev = -np.multiply(pi, dVdt)
    m_dot_in = m_dot_in #* stepwidth #/ (2 * np.pi * pV[7]) #/10
    m_dot_out = m_dot_out #* stepwidth #/ (2 * np.pi * pV[7]) #/10

    dQ = alp * ht_surface * (y[2] - Ti) * 1e-3 # kW
    dthermal_dt = state_th_Masse(y, -dQ, pV)
    dmdt = m_dot_in - m_dot_out
    dudt = (dQ + dW_fric + dW_rev - dmdt * ui - m_dot_out * hi + m_dot_in * pZ[4]) / y[0]  # kJ/kg
    return np.array([dmdt, dudt, dthermal_dt])

def getalp(pV, step, dxdt, Ti, pi):
    '''
    calculates heat transfer coefficient gas/cylinder wall
    Woschni correlation
    '''
    if step == 0 or step == 2:  # closed valves
        k = 2.28
    else:  # open valves, suction or push out
        k = 5.18
    alp = 127.93 * pV[0] ** (-.2) * (pi * 1e-2) ** .8 * (Ti) ** (-.55) * (k * abs(dxdt)) ** .8
    return alp

def state_th_Masse(y, Q, pV):
    '''
    calculates temperature change of thermal mass as function of heat transfer inside (Q)
    and to environment (Q_u)
    '''
    ### mass and cv of thermal mass are in stationary state not crucial,
    ### parameter are chosen to achieve fast convergence without vibrations
    m = .001  # kg
    cv = .502  # kJ/kg/K
    alp_a = 6.  # heat transfer coefficient to environment
    A = Ver0[3] * pV[8] / Ver0[2] * pV[0] / Ver0[0] * pV[1] / Ver0[
        1]  # Outer surface cylinder estimated via geometry related to fitting compressor
    Q_u = alp_a * A * (Tu - y[2]) * 1e-3 # kW
    dthermal_dt = (Q + Q_u) / cv / m
    return dthermal_dt

def compression(pV, pos_piston, dxdt, Ti, pi):
    step = 0
    alp = getalp(pV, step, dxdt, Ti, pi)
    m_dot_in = 0.  # no mass flow over boundaries
    m_dot_out = 0.
    return alp, m_dot_in, m_dot_out


def push_out(pV, pos_piston, dxdt, Ti, pi, pZ, pZyk, vi):
    step = 1
    alp = getalp(pV, step, dxdt, Ti, pi)
    m_dot_out = pZyk[1] / vi * np.sqrt(2. * (pi - pZ[6]) * \
                                               1000. * vi)  # mass flow leaving the cylinder, kg/s
    m_dot_in = 0.
    return alp, m_dot_in, m_dot_out


def expansion(pV, pos_piston, dxdt, Ti, pi):
    step = 2
    alp = getalp(pV, step, dxdt, Ti, pi)
    m_dot_in = 0.  # no mass flow over boundaries
    m_dot_out = 0.
    return alp, m_dot_in, m_dot_out


def suction(pV, pos_piston, dxdt, Ti, pi, pZ, pZyk):
    step = 3
    alp = getalp(pV, step, dxdt, Ti, pi)
    m_dot_in = pZyk[0] / pZ[2] * np.sqrt(
        2. * (pZ[1] - pi) * 1000 * pZ[2])  # mass flow entering cylinder, kg
    m_dot_out = 0
    return alp, m_dot_in, m_dot_out

def bc(ya, yb):
    return np.array([100000*yb[0] - 100000*ya[0], yb[1] - ya[1], yb[2] - ya[2]])

if __name__ == "__main__":
    fluid = 'Propane * Butane'
    comp = [1.0, 0.]
    p_in = z_Tx(263, 0, fluid, comp)[1]  # fl.zs_kg(['T','q'],[0.,0.],['p'],fluid)[0]
    p_out = z_Tx(355, 0, fluid, comp)[1]  # fl.zs_kg(['T','q'],[35.,0.],['p'],fluid)[0]
    T_in = 9.5 + 273.15
    resolution = 3600
    result, count, y_timetrack_m, y_timetrack_u, y_timetrack_t, is_eff, degree_delivery = set_up(T_in, p_in, p_out, fluid, comp, resolution)
    print(f"isentropic effciency: {is_eff}")
    print(f"degree of delivery: {degree_delivery}")
    time_red = result.t.reshape(-1,2).mean(axis=1)
    mass_red = result.y[0].reshape(-1,2).mean(axis=1)
    inner_energy_red = result.y[1].reshape(-1, 2).mean(axis=1)
    temperature_red = result.y[2].reshape(-1, 2).mean(axis=1)
    plt.figure(4)
    plt.plot(time_red, mass_red)
    plt.ylabel("m")
    plt.figure(5)
    plt.plot(time_red, inner_energy_red)
    plt.ylabel("u")
    plt.figure(6)
    plt.plot(time_red, temperature_red)
    plt.ylabel("T thermal")
    plt.figure(7)
    x_values = np.linspace(0,count,count+1)
    plt.plot(x_values, y_timetrack_m)
    plt.title("Timetrack mass")
    plt.figure(8)
    plt.plot(x_values, y_timetrack_u)
    plt.title("Timetrack inner energy")
    plt.figure(9)
    plt.plot(x_values, y_timetrack_t)
    plt.title("timetrack Temperature")
    print(result.message)
    #print(result.y)
    #plt.plot(np.linspace(0, 2* np.pi, resolution), result.y[1])
    plt.show()