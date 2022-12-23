# -*- coding: utf-8 -*-
"""
Compressor model of Dennis Roskosch
changed to work with the actual Refprop module 10.12 2022 (B. Atakan)
It uses extensivly global variables (fixed 19.12.2022 Welp) This should be changed and rewritten
to use numba and an ODE-solver

Created on Thu Jan 31 17:46:25 2019

@author: roskosch
"""
import numpy as np
import matplotlib.pyplot as plt
from carbatpy.fl_props_compressor import z_uv, z_ps, z_Tp, z_Tx, z_mm

Rm = 8.3145  # gas constant J/mol/K
Tu = 25. + 273.15  # ambient temperature
dTK = 273.15  # conversion °C / K
Ver0 = [34e-3, 34e-3, 2., .04]  # fit-compressor: D, H, cylinder, outer surface
# TODO what is fit-compressor? duplication with pV!


###########################################
#  compressor specific parameter
# pV:    0 - D, piston diameter, m
#        1 - H, stroke, m
#        2   ratio of crankshaft radius and connecting rod length
#        3   A_a, surface compressor, m²
#        4 - c1, dead volume relative to displacement
#        5 - pfric, friction pressure, kPa
#        6 - f, electrical compressor frequency Hz
#        7 - compressor speed Hz
#        8 - number of cylinders


# pZ:
# Index S: suction side
# Index D: pressure side
#       0 - TS, °C  # BA: seems to be K # AW: is initialized in line getETA as K
#       1 - pS, kPa
#       2 - vS, m³/kg
#       3 - uS, kJ/kg
#       4 - hS, kJ/kg
#       5 - sS, J/kg/K
#       6 - pD, kPa

# pZyk: quantities that remain constant during one cycle or are iterated over different cycles respectively
#       0 - Aeff_i, effective flow cross-section inlet, m²
#       1 - Aeff_o, effective flow cross-section outlet, m²


# z_it: stores quantities for every iteration step
#       0 - tet, crank angle
#       1 - x, piston position
#       2 - V, cylinder volume
#       3 - surface heat transfer cylinder wall
#       4 - cycle step, 0=compression, 1=push out, 2=expansion, 3=suction
#       5 - T, temperature in cylinder °C  # BA: seems to be K
#       6 - p, pressure in cylinder kPa
#       7 - v, specific volume in cylinder m³/kg
#       8 - u, internal energy in cylinder kJ/kg
#       9 - h, enthalpy in cylinder kJ/kg
#      10 - s, entropy in cylinder J/kg/K
#      11 - m, mass in cylinder kg
#      12 - T_th, temperature of thermal mass
#      13 - alpha, inner heat transfer coefficient
#      14 - dm, inflowing or outflowing mass in iteration step kg
#      15 - transferred heat


def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res


def getalp(z_it, i, pV):
    '''
    calculates heat transfer coefficient gas/cylinder wall
    Woschni correlation
    '''
    if z_it[i, 4] == 0 or z_it[i, 4] == 2:  # closed valves
        k = 2.28
    else:  # open valves, suction or push out
        k = 5.18
    v_p = np.abs(z_it[i, 1] - z_it[i - 1, 1]) / ((z_it[i, 0] - z_it[i - 1, 0]) /
                                                 (2. * np.pi * pV[7]))  # dX/dt
    alp = 127.93 * pV[0] ** (-.2) * (z_it[i - 1, 6] * 1e-2) ** .8 * \
          (z_it[i - 1, 5]) ** (-.55) * (k * v_p) ** .8
    z_it[i, 13] = alp
    return z_it


def state_th_Masse(Q, z_it, i, pV):
    '''
    calculates temperature change of thermal mass as function of heat transfer inside (Q)
    and to environment (Q_u)
    '''
    ### mass and cv of thermal mass are in stationary state not crucial,
    ### parameter are chosen to achieve fast convergance without vibrations
    m = .0001  # kg
    cv = .502  # kJ/kg/K
    alp_a = 6.  # heat transfer coefficient to environment
    A = Ver0[3] * pV[8] / Ver0[2] * pV[0] / Ver0[0] * pV[1] / Ver0[
        1]  # Outer surface cylinder estimated via geometry related to fitting compressor
    Q_u = alp_a * A * (Tu - z_it[i - 1, 12]) * ((z_it[i, 0] - z_it[i - 1, 0]) /
                                                (2. * np.pi * pV[7])) * 1e-3  # kJ
    z_it[i, 12] = (Q + Q_u) / cv / m + z_it[i - 1, 12]
    return z_it


def compression(i, fluid, z_it, comp, pV):
    step = 0
    W = -z_it[i - 1, 6] * (z_it[i, 2] - z_it[i - 1, 2])  # compression work, kJ
    Wr = -pV[5] * (z_it[i, 2] - z_it[i - 1, 2])  # friction work, kJ
    z_it = getalp(z_it, i, pV)
    Q = z_it[i, 13] * z_it[i, 3] * (z_it[i - 1, 12] - z_it[i - 1, 5]) * \
        ((z_it[i, 0] - z_it[i - 1, 0]) / (2. * np.pi * pV[7])) * 1e-3  # kJ
    z_it = state_th_Masse(-Q, z_it, i, pV)
    dm = 0.  # no mass flow over boundaries
    mi = z_it[i - 1, 11]  # mass inside cylinder, kg
    ui = (Q + W + Wr) / mi + z_it[i - 1, 8]  # kJ/kg
    vi = z_it[i, 2] / mi  # specific volume in cylinder, m³/kg
    # property call state i with u and v, return: T, p, v, u, h, s
    zi = z_uv(ui, vi, fluid, comp)  # fl.zs_kg(['u','v'],[ui,vi],['T','p','v','u','h','s'],fluid)

    # store of quantities in z_it
    z_it[i, 4] = step
    z_it[i, 5:11] = zi
    z_it[i, 11] = mi
    z_it[i, 14:16] = dm, Q
    return z_it


def push_out(i, fluid, z_it, comp, pV, pZyk, pZ):
    step = 1
    W = -z_it[i - 1, 6] * (z_it[i, 2] - z_it[i - 1, 2])  # compression work, kJ
    Wr = -pV[5] * (z_it[i, 2] - z_it[i - 1, 2])  # friction work, kJ
    z_it = getalp(z_it, i, pV)
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
    W = -z_it[i - 1, 6] * (z_it[i, 2] - z_it[i - 1, 2])  # compression work, kJ
    Wr = pV[5] * (z_it[i, 2] - z_it[i - 1, 2])  # friction work, kJ
    getalp(z_it, i, pV)
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
    W = -z_it[i - 1, 6] * (z_it[i, 2] - z_it[i - 1, 2])  # compression work, kJ
    Wr = pV[5] * (z_it[i, 2] - z_it[i - 1, 2])  # friction work # , kJ
    getalp(z_it, i, pV)
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


def process_iteration(fluid, pZyk, z_it, IS, IS0, comp, pV, pZ):
    # setting of Aeff_i, explicit function
    M = z_mm(300, 100., fluid, comp)[-1]  # CP.PropsSI("M",fluid) # molar mass kg/mol # AW ramdom inlet conditions, molar mass contant
    pZyk[0] = 2.0415e-3 * (Rm / M) ** (-.9826) * pV[0] ** 2. / Ver0[
        0] ** 2.  # effective flow cross-section inlet, m²
    # setting of Aeff_o, implicit function relatively to average mass flow density over valve
    # at 1st iteration, the mass flow density is unknown, typical value is guessed
    pZyk[1] = 1.5e-5 * pV[0] ** 2. / Ver0[0] ** 2.
    # print(pZyk)
    count = 0

    while 1:
        count += 1
        for i in range(1, IS):
            if z_it[i, 0] <= np.pi:
                if z_it[i - 1, 6] <= pZ[6]:
                    z_it = compression(i, fluid, z_it, comp, pV)
                else:
                    z_it = push_out(i, fluid, z_it, comp, pV, pZyk, pZ)
            else:
                if z_it[i - 1, 6] >= pZ[1]:
                    z_it = expansion(i, fluid, z_it, comp, pV)
                else:
                    z_it = suction(i, fluid, z_it, comp, pV, pZyk, pZ)

        # error square sum T, p, T_th_average
        error = np.sqrt((z_it[-1, 5] - z_it[0, 5]) ** 2.) + np.sqrt((z_it[-1, 6]
                                                                     - z_it[0, 6]) ** 2.) + np.sqrt(
            (np.average(z_it[-1, 12])
             - np.average(z_it[0, 12])) ** 2.)
        # print(IS,count)

        if error < .01:  # two-step: first many calculations with low resolution (BA)
            if IS == IS0:
                IS = 10 * IS0  # then higher resolution
                z0_ = z_it[:, :]
                z_it = np.zeros([IS, 16])
                z_it[:IS0, :] = z0_
                z_it[IS0:, :] = z0_[-1, :]
                geometry(pV, pZ, z_it, fluid, IS)  # calculate angle correctly
            else:
                break
        else:
            cell_push_out = find(z_it[:, 4] == 1)
            m_aus = np.sum(z_it[cell_push_out, 14])  # overall pushed out mass
            t_aus = (z_it[cell_push_out[-1], 0] -
                     z_it[cell_push_out[0], 0]) / (2. * np.pi * pV[7])  # time of push out
            m_dichte = m_aus / t_aus / pZyk[1]  # mass flow density kg/s/m²
            # print("BA:",m_dichte,pV[0], error)
            pZyk[1] = 5.1109e-4 * (m_dichte) ** (-.486) * pV[0] ** 2. / Ver0[0] ** 2.  # Aeff_o neu
            z_it[0, 5:14] = z_it[-1, 5:14]  # End values of last cycle = Start values of next cycle

    # Efficiency evaluation
    cell_push_out = find(z_it[:, 4] == 1)  # BA find()
    m_aus = np.sum(z_it[cell_push_out, 14])  # overall pushed out mass
    m0 = np.pi * pV[0] ** 2. * pV[1] / pZ[2] / 4.  # sucked-in mass ideal compressor
    degree_delivery = m_aus / m0  # degree of delivery

    h_aus = np.sum(z_it[cell_push_out, 9] * z_it[cell_push_out, 14]) \
            / m_aus  # average push out enthalpy
    h_aus_s = z_ps(pZ[6], pZ[5], fluid, comp)[
        4]  # fl.zs_kg(['p','s'],[pZ[6],pZ[5]],['h'],fluid)[0]  # isentropic outlet enthalpy
    is_eff = (h_aus_s - pZ[4]) / (h_aus - pZ[4])  # isentropic efficiency

    return is_eff, degree_delivery


def getETA(T_e, p_e, p_a, fluid_in, comp, pV, pZ, z_it, IS, pZyk, IS0):
    fluid = fluid_in
    comp = comp

    ###############################     Parametersatz spezifisch für Verdichter   ###################################

    pV = [34e-3, 34e-3, 3.5, .04, .06071, 48.916, 50., 50. / 2., 2.]  # parameter see above

    #################################################################################################################
    pZ[0:6] = z_Tp(T_e, p_e, fluid,
                   comp)  # fl.zs_kg(['T','p'],[T_e,p_e],['T','p','v','u','h','s'],fluid) #state suction pipe
    pZ[6] = p_a  # pressure in pressure pipe
    print(pZ)
    ############### set geometry ##################################
    z_it[:, 0] = np.linspace(0., 2 * np.pi, IS)
    z_it[:, 1] = -(pV[1] / 2. * (1. - np.cos(z_it[:, 0]) + pV[2] *
                                 (1. - np.sqrt(1. - (1. / pV[2] * np.sin(z_it[:, 0])) ** 2.)))) + \
                 pV[4] * pV[1] + pV[1]  # piston position, x=0 at UT
    a_head = np.pi / 4. * pV[0] ** 2.   # area of cylinder head
    z_it[:, 2] = a_head * z_it[:, 1]  # volume cylinder
    z_it[:, 3] = np.pi * pV[0] * z_it[:, 1] + 2. * a_head  # heat transfer surfaces
    #plt.plot(z_it[:,0], z_it[:, 1]) #AW
    #plt.show() #AW
    ########## set start conditions in the cylinder ##################
    z_it[0, 5:11] = pZ[0:6]
    z_it[0, 11] = z_it[0, 2] / z_it[0, 7]  # V/v, Cylinder completely filled with suction gas
    ##### start temperature thermal mass, only one temperature per cylce
    ##### Start value freely selectable, significantly influences iteration time
    z_it[:, 12] = 42. + 273
    is_eff, degree_delivery = process_iteration(fluid, pZyk, z_it, IS, IS0, comp, pV, pZ)
    return np.array((is_eff, degree_delivery))


def geometry(pV, pZ, z_it, fluid, IS):
    z_it[:, 0] = np.linspace(0., 2 * np.pi, IS)
    z_it[:, 1] = -(pV[1] / 2. * (1. - np.cos(z_it[:, 0]) + pV[2] *
                                 (1. - np.sqrt(1. - (1. / pV[2] * np.sin(z_it[:, 0])) ** 2.)))) + \
                 pV[4] * pV[1] + pV[1]  # piston position, x=0 bei UT
    A_kopf = np.pi / 4. * pV[0] ** 2.
    z_it[:, 2] = A_kopf * z_it[:, 1]  # cylinder volume
    z_it[:, 3] = np.pi * pV[0] * z_it[:, 1] + 2. * A_kopf  # heat transfer surfaces


# Beispiel #################################################
if __name__ == "__main__":
    IS = 360  # number of differential steps for one cycle
    IS0 = IS
    pV = np.zeros(8, float)
    pZ = np.zeros(7, float)
    pZyk = np.zeros(2, float)
    z_it = np.zeros([IS, 16])
    # fluid = []
    # comp = [1.0]  # must be checked BA

    fluid = 'Propane * Butane'
    comp = [1.0, 0.]
    pe = z_Tx(263, 0, fluid, comp)[1]  # fl.zs_kg(['T','q'],[0.,0.],['p'],fluid)[0]
    pa = z_Tx(355, 0, fluid, comp)[1]  # fl.zs_kg(['T','q'],[35.,0.],['p'],fluid)[0]

    fo = open("Daten.txtx", "w")
    print("Drücke %2.2f kPa %2.2f kPa" % (pe, pa))
    dt_all = np.linspace(9.5, 20.5, 3)
    out = []
    for dt in dt_all:
        o1 = getETA(dt + 273.15, pe, pa, fluid, comp, pV, pZ, z_it, IS, pZyk, IS0)
        # o1.append((np.max(z_it[:,11]) - np.min(z_it[:,11]) * pV[7]))  # mass flow
        out.append(o1)
        print(dt, o1)
    out = np.array(out)
    plt.plot(dt_all, out)
    plt.show()
    fo.write(str(out))
    fo.close()
