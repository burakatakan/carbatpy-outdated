# -*- coding: utf-8 -*-
"""
Heat pump cycle incl. heat transfer.

Calculate the pressures and the mass flow rate to fulfill the heat transfer
canditions for given (constant!) temperatures of the heat sink and source
x=1 before the compressor with isentropic efficienc eta and x=0 entering
the throttle.
This is done by root finding for a fluid, taking the properties from CoolProp
Will be used for sensitivity analysis and optimization (hopefully).


Created on Tue Feb 23 17:15:50 2021

@author: atakan
"""

import CoolProp.CoolProp as CP

import fluid_properties_rp as fprop
import numpy as np
# import ht
import scipy.optimize as opti
import matplotlib.pyplot as plt
import os
# p_0 = 1.013e5
# temp_0 = 300
    
fluid_a = "Propane"
composition = [1.0]
temp = np.array([279, 330])  # source and sink Temperature
eta = .65
p = np.array([6, 15]) * 1e5


FLUIDMODEL = "REFPROP"
if FLUIDMODEL == "REFPROP":
    os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
    from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
    _props ="REFPROP"  # or "CoolProp"
    working_fluid = fluid_a
    if _props =="REFPROP":
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        _units = RP.GETENUMdll(0, "MASS BASE SI").iEnum # be careful pressure is in Pa!
    else: _units = 21

elif FLUIDMODEL == "CoolProp":
    # h_0 = CP.PropsSI('H', 'P', p_0, 'T', temp_0, fluid_a)
    working_fluid = CP.AbstractState("HEOS", fluid_a)
    



def prop_pt(p, _temp, fluid):
    """
    Properties needed for integration at given p and h, single phase.

    Parameters
    ----------
    p : float
        pressure in Pa.
    _temp : float
        Temperature in K.
    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

    """
    fluid.update(CP.PT_INPUTS, p, _temp)
    reihe = [CP.iHmass, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
             CP.iconductivity, CP.iCpmass, CP.iviscosity,
             CP.iisobaric_expansion_coefficient, ]
    props = [fluid.keyed_output(k) for k in reihe]
    h, x, s, rho, prandtl, phase, lambda_s, cp, mu, beta = props[:]
    alle = [_temp, p, x, h,  s, rho, mu,
            cp, lambda_s, phase, prandtl, beta * _temp]
    return alle


def prop_ps(p, s, fluid):
    """
    Properties needed for integration at given p and h, single phase.

    Parameters
    ----------
    p : float
        pressure in Pa.
    s : float
        specific entropy in J/kg/K.
    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

    """
    fluid.update(CP.PSmass_INPUTS, p, s)
    reihe = [CP.iHmass, CP.iQ, CP.iT, CP.iDmass, CP.iPrandtl, CP.iPhase,
             CP.iconductivity, CP.iCpmass, CP.iviscosity,
             CP.iisobaric_expansion_coefficient, ]
    props = [fluid.keyed_output(k) for k in reihe]
    h, x, _temp, rho, prandtl, phase, lambda_s, cp, mu, beta = props[:]
    alle = [_temp, p, x, h,  s, rho, mu,
            cp, lambda_s, phase, prandtl]
    return alle


def prop_ph(p, h, fluid):
    """
    Properties needed for integration at given p and h, single phase.

    Parameters
    ----------
    p : float
        pressure in Pa.
    h : float
        specific enthalpy in J/kg.
    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

    """
    fluid.update(CP.HmassP_INPUTS, h, p)
    reihe = [CP.iT, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
             CP.iconductivity, CP.iCpmass, CP.iviscosity]
    props = [fluid.keyed_output(k) for k in reihe]
    _temp, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]
    alle = [_temp, p, x, h,  s, rho, mu,
            cp, lambda_s, phase, prandtl]

    return alle


def prop_pq(p, q, fluid):
    """
    Properties needed for integration at given p and h, single phase.

    Parameters
    ----------
    p : float
        pressure in Pa.
    q : float
        quality (0-1).
    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

    """
    fluid.update(CP.PQ_INPUTS, p, q)
    reihe = [CP.iT, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
             CP.iconductivity, CP.iCpmass, CP.iviscosity, CP.iHmass]
    props = [fluid.keyed_output(k) for k in reihe]
    _temp, x, s, rho, prandtl, phase, lambda_s, cp, mu, h = props[:]
    alle = [_temp, p, x, h,  s, rho, mu,
            cp, lambda_s, phase, prandtl]

    return alle


def check_bound(p, bounds):
    for i in range(len(p)):
        if not((p[i] >= bounds[i][0]) and (p[i] <= bounds[i][1])):
            return False
    return True



def heat_pump_opti(para, eta, U, T_s, working_fluid, bounds):
    p = para[:3]
    areas = para[3:]
    #print("p,a",p, areas)
    loes = opti.root(heat_pump_ht, p,
                     args=(eta, U, areas, T_s, working_fluid,True, False),
                     options={"factor": .25})
    # print(loes)
    if loes.success:
        p = loes.x[:3]
        A = areas
        # print(loes.x, p,A)
        rat_eff = heat_pump_ht(p, eta, U, A, T_s, working_fluid, 
                               rootfind =False, optim=True)
        return -rat_eff
    else:
        return 1e6
        
    



def heat_pump_ht(p, eta, U, A, T_s, working_fluid, rootfind =True, optim=False,
                 bounds=None):
    """
    Heat pump cycle with isothermal source and sink, x=0 after condenser.

    x=1 entering the compressor. For optimizing the two working fluid pressures
    or calculating states, COPs etc

    Parameters
    ----------
    p : numpy.array, 3 floats
         with the pressures of evaporator and condenser and  the
        mass flow rate (all SI Pa, kg/s.
    eta : float
        compressor efficienvy (isentropic).
    U : numpy.array, 3floats
        mean overall heat ransfer coefficients for evaporator, condenser(gas)
        condenser two-phase (W/m2/K).
    A : numpy.array, 2 floats
        areas of evaporatorand condenser (m2).
    T_s : numpy.array, 2 floats
        secondary fluid temperatures low and high in K.
    working_fluid : coolprop-fluid
        DESCRIPTION.
    rootfind: Boolean, optional.
        are we finding the root? (T-differences according to heat exchanger size)
    optim : Boolean, optional
        if true: output for optimization
        if false: output of several properties. The default is True.
    bounds: not yet implemented

    Returns
    -------
    if rootfind == True:
        array with three differences of energy balances
    if optim == True:
        rational efficiency (for optimization)
    if False:
          state: array with all calculated states (T,p,x,h,s,rho ...)
          q_w: all specific energies transfered (J/kg)
          A_hg: heat exchanger Area until condensation
          cop: COP of the cycle
          cop_carnot: Carnot-COP for the tempereatures of condenser and
              evaporator (superheating is neglected)
          eta_rational: ratio of the latter two COPs
          p_compr: Compressor Power (W)
          q_h_dot: Heat flow rate at high T (W).

    """
    problem = 0
    if bounds!= None:
        inside = check_bound(p, bounds)
    else:
        inside = True
    if inside:
        p = np.abs(p)
        m_dot = p[2]
        if _props == "REFPROP":
            state = np.zeros((8, 6))
        else:
            state = np.zeros((8, 11))
        diff = np.zeros((3))
        q_w = np. zeros((4))
        try:
            state[0, :] = fprop.prop_pq(p[0], 1, working_fluid, composition, 
                                    option=1, units =_units, props=_props)  # vor dem Kompressor
            if state[0, 0] >= T_s[0]: problem =1
        except:
            state[0, :] = fprop.prop_pq(1.3e5, 1, working_fluid, composition, 
                                    option=1, units =_units, props=_props)  # vor dem Kompressor
            print("Fehler!:", p,eta, U, A, T_s, optim)
            
        # Kompressor
        state[1, :] = fprop.sp(state[0, 4], p[1], working_fluid, composition, 
                                option=1, units =_units, props=_props)  # isentropic
        q_w[0] = (state[1, 2] - state[0, 2]) / eta
        
        # kondensatoreintritt:
        state[2, :] = fprop.hp(state[0, 2] + q_w[0], p[1],  working_fluid, composition, 
                                option=1, units =_units, props=_props)
        if state[2, 0] < T_s[1]: problem = 2
        # print ("Problem:", problem, p)
        state[6, :] = fprop.prop_pq(p[1], 1, working_fluid, composition, 
                                option=1, units =_units, props=_props)  # Taupunkt
        if state[6, 0] <= T_s[1]: problem = 3
        q_w[1] = (state[6, 2] - state[2, 2])  # Waermeuebtragung bis dahin

        state[7, :] = fprop.prop_pq(p[1], 0, working_fluid, composition, 
                                option=1, units =_units, props=_props)  # gerade kondensiert, Wunsch
        q_w[2] = state[7, 2] - state[6, 2]  # isotherme Kondensation
        state[3, :] = fprop.hp(state[7, 2], p[0], working_fluid, composition, 
                                option=1, units =_units, props=_props)  # hinter Drossel
        q_w[3] = state[0, 2] - state[3, 2]  # Verdampfer
        # - heat transfer:
        q_v = U[0] * A[0] * (T_s[0] - state[0, 0]) / m_dot  # isotherm Verdampfer
        A_hg = q_w[1] / U[1] * m_dot  # Gas-abkühlung, Fläche
        A_hg = -A_hg / ((state[2, 0] - state[6, 0]) /
                      np.log(np.abs((T_s[1] - state[2, 0]) /
                             (T_s[1] - state[6, 0]))))
        q_kond = U[2] * (A[1] - A_hg) * (T_s[1] - state[7, 0]) / m_dot  # isotherme Kondensation

        diff[0] = q_w[3] - q_v  # passt der Druck im Verdampfer zum Waermestrom?
        diff[1] = (state[7, 2] - state[2, 2]) - (q_kond + q_w[1])  # Kondensation incl. Gas
        diff[2] = state[7, 2] - state[6, 2] - q_kond  # nur Kondensation
        # print (q_w, q_kond, q_v)
        if problem > 0:
            diff[:] = 1e7
        cop = -q_w[1:3].sum()/q_w[0]
        cop_carnot = T_s[1] / (T_s [1] - T_s[0])
        eta_rational = cop / cop_carnot
        if rootfind:
            return diff
        elif optim:
            print(eta_rational, m_dot, p,A ,"opti" )
            return eta_rational
        else:
            
            p_compr = q_w[0] * m_dot
            q_h_dot = q_w[1:3].sum() * m_dot
            state =np.delete(state, [1,4,5],0) #  so far unused pointes deleted

            return state, q_w, A_hg, cop, cop_carnot, eta_rational,\
                p_compr,  q_h_dot
    else:
        if rootfind:
            return np.ones(3) * 1e5
        else:
            return np.zeros(14)

if __name__ == "__main__":
    diameter = 15e-3  # tube diameter in m
    schleif = False
    optimierung = True
    U = np.array([1300, 250, 1300])
    areas = np.array([(12 * np.pi * diameter), (28 * np.pi * diameter)]) * 1.3 # sehr nichtlinear!
    p0 = np.array([2e5, 25e5, 0.018])  # Propan np.array([4e5, 19e5, 0.008])

    loes = opti.root(heat_pump_ht, p0,
                     args=(eta, U, areas, temp, working_fluid),
                     options={"factor": .25})
    if loes.success:
    
        st, qw, A_hg, cop, cop_carnot, eta_rat, p_compr, q_h_dot = \
            heat_pump_ht(loes.x, eta, U, areas, temp, working_fluid, False)
        sortierung =[0, 1, 3, 4, 2, 0]  # order of the states along the cycle
        plt.axhline(temp[0])
        plt.axhline(temp[1])
        plt.plot(st[sortierung,4],st[sortierung,0], "-o")
    
        print("COP: %2.2f, COP(Carnot): %2.2f, eta(rational): %2.2f"
              % (cop, cop_carnot, eta_rat))
        print("P-Compressor/W: %3.2f, Q_dot-highT/W: %3.2f"
              % (p_compr, q_h_dot))
        
    if optimierung:
        bou = [(5e4, 4e5), (14e5, 2.5e6), (0.017, 0.023), (0.8, 3), (1, 5)]
        
        opti_loes = opti.shgo(lambda x:  heat_pump_opti(x,eta, U, temp, 
                                                        working_fluid, bou), bou)
                                    
        print("\n\n",opti_loes)
        
        
            
    # if schleif:
    #     aa = []
    #     for a_factor in np.arange(1, 100, .75):
    #         areas = np.array([(1.9 * np.pi * 10e-3), (1 * np.pi * 10e-3)]) * a_factor
    #         try:
    #             loes = opti.root(heat_pump_ht, p0,
    #                              args=(eta, U, areas, temp, working_fluid),
    #                              options={"factor": .25})
    
    #             st, qw, A_hg, cop, cop_carnot, eta_rat, p_compr, q_h_dot = \
    #                 heat_pump_ht(loes.x, eta, U, areas, temp, working_fluid, False)
    
    #             print("COP: %2.2f, COP(Carnot): %2.2f, eta(rational): %1.3f"
    #                   % (cop, cop_carnot, eta_rat))
    #             print("P-Compressor/W: %3.2f, Q_dot-highT/W: %3.2f"
    #                   % (p_compr, q_h_dot))
    #             p0 = loes.x
    #         except:
    #             pass
    #         aa.append([a_factor, eta_rat, p_compr, cop, q_h_dot, *loes.x])
    #     aa = np.array(aa)
    
    #     f, ax = plt.subplots(2, 2, sharex=True)
    #     ax[0, 0].plot(aa[:, 0], aa[:, 1], "v")
    #     ax[0, 0].set_title("$\eta_{rat}$")
    #     ax[1, 1].plot(aa[:, 0], aa[:, 2]/1000, ".-")
    #     ax[1,  1].plot(aa[:, 0], -aa[:, 4] / 1000, "-")
    #     ax[1, 1].set_title("P, $\dot Q_h$ / kW")
    #     ax[1, 0].plot(aa[:, 0], aa[:, -3]/1e5, "v-")
    #     ax[1, 0].plot(aa[:, 0], aa[:, -2]/1e6, "-")
    #     ax[1, 0].set_title("$p_l, p_h/10$ / bar")
    #     ax[0, 1].plot(aa[:, 0], aa[:, -1]*1000, "o")
    #     ax[0, 1].set_title("$\dot m$ / g/s")
    #     f.suptitle(fluid_a+", T(l):%3i K, T(h):%3i K" % (temp[0], temp[1]))
    #     for axx in ax.flat:
    #         axx.set(xlabel='tubes')  # , ylabel='y-label')
    # else:
    #     print(loes)
