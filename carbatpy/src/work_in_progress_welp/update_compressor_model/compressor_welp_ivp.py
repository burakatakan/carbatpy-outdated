"""
compressor model adapted from Dennis Roskosch
changed to use ODE solver instead of finite volume approach

author: Alexandra Welp
21.12.2022
"""

import numpy as np
import matplotlib.pyplot as plt
from fl_props_compressor import z_uv, z_ps, z_Tp, z_Tx, z_mm
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
    number_of_cycles = 20
    x_max = number_of_cycles / pV[7]
    y_start = np.zeros(3)
    y_start[0] = Ver0[1] * a_head / pZ[2]       #m
    y_start[1] = pZ[3]                          #u
    y_start[2] = Tu                             #T
    #Ti = np.zeros(len(x_var))
    #pi = np.zeros(len(x_var))
    #hi = np.zeros(len(x_var))
    #si = np.zeros(len(x_var))
    #m_dot_in = np.zeros(len(x_var))
    #m_dot_out = np.zeros(len(x_var))
    #alp = np.zeros(len(x_var))
    #Ti, pi, hi, si, m_dot_in, m_dot_out, alp
    res = solve_ivp(fun, [0, x_max], y_start, method='RK45', args=[pV, a_head, fluid, comp, pZ, pZyk], max_step=1/(pV[7]*3600))
    return res

def initialization(x_length):
    si = np.zeros(x_length)
    Ti = np.zeros(x_length)
    alp = np.zeros(x_length)
    hi = np.zeros(x_length)
    m_dot_in = np.zeros(x_length)
    m_dot_out = np.zeros(x_length)
    pi = np.zeros(x_length)
    ui = np.zeros(x_length)
    return si, Ti, alp, hi, m_dot_in, m_dot_out, pi, ui

def fun(x, y, pV, a_head, fluid, comp, pZ, pZyk):
    #si, Ti, alp, hi, m_dot_in, m_dot_out, pi, ui = initialization(len(x))
    #for i in range(0, len(x)):
    if y[0] < 0:
        y[0] = 0.0002
        #if y[0, i] > 0.0003:
            #y[0, i] = 0.0002

        #if y[2, i] < 0:
            #y[2, i] = Tu
    theta = x * (2 * np.pi * pV[7])
    while theta >= 2 * np.pi:
        theta = theta - 2 * np.pi
    pos_piston = -(pV[1] / 2. * (1. - np.cos(theta) + pV[2] *
                    (1. - np.sqrt(1. - (1. / pV[2] * np.sin(theta)) ** 2.)))) + pV[4] * pV[1] + pV[1]  # piston position, x=0 at UT
    volume_cylinder = a_head * pos_piston  # volume cylinder
    ht_surface = np.pi * pV[0] * pos_piston + 2. * a_head  # heat transfer surfaces
    vi = volume_cylinder / y[0]  # specific volume in cylinder, m³/kg
    dxdtheta = -pV[1] / 2 * np.sin(theta) * (1 + 1/pV[2] * np.cos(theta) * (1 - (1/pV[2] * np.sin(theta))**2)**-0.5)
    dxdt = (2 * np.pi * pV[7]) * dxdtheta
    dVdt = a_head * dxdt
    #for i in range(0, len(x)):
    [Ti, pi, vi, ui, hi, si] = z_uv(y[1], vi, fluid, comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)
    if Ti == -9999990.:
            raise ValueError("invalid properties")
    if theta <= np.pi:
        dW_fric = - pV[5] * dVdt
        if pi <= pZ[6]:
            [alp, m_dot_in, m_dot_out] = compression(pV, pos_piston, dxdt, Ti, pi)
        else:
            [alp, m_dot_in, m_dot_out] = push_out(pV, pos_piston, dxdt, Ti, pi, pZ, pZyk, vi)
                #plt.figure(5)
                #plt.plot(i, alp[i],"*r", label="alpha")
                #plt.figure(6)
                #plt.plot(i, m_dot_out[i],"*b", label="m_dot_out")
                #plt.plot(i, dxdt[i],"*c", label="dxdt")
                #plt.plot(i, vi[i],"*y", label="vi")
                #plt.figure(7)
                #plt.plot(i, pi[i],"*m",label="pi")
                #if i == 1790:
                    #plt.legend()
            mass_flow_density = m_dot_out / pZyk[1]
            pZyk[1] = 5.1109e-4 * mass_flow_density ** -.486 * pV[0] ** 2 / Ver0[0] ** 2  # Aeff_o neu
    else:
        dW_fric = pV[5] * dVdt
        if pi >= pZ[1]:
                [alp, m_dot_in, m_dot_out] = expansion(pV, pos_piston, dxdt, Ti, pi)
        else:
                [alp, m_dot_in, m_dot_out] = suction(pV, pos_piston, dxdt, Ti, pi, pZ, pZyk)



    dW_rev = -np.multiply(pi, dVdt)
    #stepwidth = x[1] - x[0]
    #print(f"m_aus {sum(m_dot_out)*stepwidth / (2 * np.pi * pV[7])}")
    #print(f"m_in {sum(m_dot_in)*stepwidth / (2 * np.pi * pV[7])}")
    m_dot_in = m_dot_in #* stepwidth #/ (2 * np.pi * pV[7]) #/10
    m_dot_out = m_dot_out #* stepwidth #/ (2 * np.pi * pV[7]) #/10

    dQ = alp * ht_surface * (y[2] - Ti) * 1e-3 # kW
    dthermal_dt = state_th_Masse(y, -dQ, pV)
    dmdt = m_dot_in - m_dot_out
    dudt = (dQ + dW_fric + dW_rev - dmdt * ui - m_dot_out * hi + m_dot_in * pZ[4]) / y[0]  # kJ/kg
    plt.figure(1)
    plt.plot(x,dudt)
    plt.title("dudt")
    plt.figure(2)
    plt.plot(x,dmdt)
    plt.title("dmdt")
    plt.figure(3)
    plt.plot(x,dthermal_dt)
    plt.title("dTdt")
    #plt.show()
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
    m = 1000.  # kg
    cv = .502  # kJ/kg/K
    alp_a = 6.  # heat transfer coefficient to environment
    A = Ver0[3] * pV[8] / Ver0[2] * pV[0] / Ver0[0] * pV[1] / Ver0[
        1]  # Outer surface cylinder estimated via geometry related to fitting compressor
    Q_u = alp_a * A * (Tu - y[2]) * 1e-3 # kW
    dthermal_dt = (-Q + Q_u) / cv / m
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
    result = set_up(T_in, p_in, p_out, fluid, comp, resolution)
    plt.figure(4)
    plt.plot(result.t[3600*19::], result.y[1,3600*19::])
    plt.ylabel("u")
    plt.figure(5)
    plt.plot(result.t[3600*19::], result.y[0,3600*19::])
    plt.ylabel("m")
    plt.figure(6)
    plt.plot(result.t[3600*19::], result.y[2,3600*19::])
    plt.ylabel("T thermal")

    print(result.message)
    print(result.y)
    #plt.plot(np.linspace(0, 2* np.pi, resolution), result.y[1])
    plt.show()