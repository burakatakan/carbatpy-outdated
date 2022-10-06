# -*- coding: utf-8 -*-
"""
functions to obtain fluid properties from cool prop using the low level interface
so far for known T and P (tp()) and for known h and p (hps(). The latter
                                                       is also vectorized 
by hand (hps_v())

The funcions for the saturated states must be checked.
Created on Wed Dec  9 17:37:20 2020
changes 04.10.2022

@author: atakan
"""


import numpy as np
import CoolProp.CoolProp as CP
import os

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
_props ="REFPROP"  # or "CoolProp"
if _props =="REFPROP":
    RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    _units = RP.GETENUMdll(0, "MASS BASE SI").iEnum # be careful pressure is in Pa!
else: _units = 21


__Tenv__ = 283.15 # Temp. of the environment in K
__penv__= 1.013e5  # Pressure of the environment in Pa

def mdot_area_function(m_dot, diameter):
    area = np.pi * (diameter / 2)**2
    m_dot_area = m_dot / area
    return m_dot_area, area

def hp_exergy(h, p, fluid, T_env=__Tenv__, p_env=__penv__, props=_props, 
              composition=[1.0]):
    pr_test = False
    state_env = tp(T_env, p_env, fluid, option =1, props=props, 
                   composition=composition)
    state = hps(h, p, fluid, props=props, composition=composition) 
    dstate = state - state_env
    ex = dstate[2] - __Tenv__ * dstate[4]
    if pr_test: print(state, state_env,dstate,"\n")
    return ex

def hps(h, p, fluid, composition=[1.0], option=1, units =_units, props=_props):
    """
    Properties needed for integration at given p and h, single phase.

    Parameters
    ----------
    h : float
        specific enthalpy in J/kg.
    p : float
        pressure in Pa.
    
    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        if option = 0 it includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.
        
        if option =1: it is compatible to the high level output
        T,p,h,v,s,x.

    """
    if props == "REFPROP":
        o = RP.REFPROP2dll(fluid,"HP","T;D;S;q", units, 0, h, p, composition)
        if option == 0:
            alle =[]
        if option == 1:
            alle =[o.Output[0], p, h, 1/o.Output[1], *o.Output[2:4]]
            
    elif props =="CoolProp":
        fluid.update(CP.HmassP_INPUTS, h, p)
        reihe = [CP.iT, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                 CP.iconductivity, CP.iCpmass, CP.iviscosity]
        props = [fluid.keyed_output(k) for k in reihe]
        _temp, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]
        
        if option == 0:
            alle = [_temp, p, x, h,  s, rho, mu,
                    cp, lambda_s, phase, prandtl]
        elif option == 1:
            alle = [_temp, p,  h, 1/rho, s, x]
    return np.array(alle)
         
    
    
name_properties = [
    ["temperature", "p", "x", "h",  "s", "rho", "mu", "cp", "lambda_s", 
     "phase", "prandtl"],
    ["temperature", "p",  "h", "v", "s","x"] 
    ]

def hps_v(h, p, fluid, composition=[1.0], option=1, units =_units, props=_props):
    """ Vectorization of the single phase properties function"""
    _n = len(h)
    
    if option ==1:
        alle = np.zeros((6, _n))
    else:
        alle = np.zeros((11, _n))
    for _i in range(_n):
        if type(p) is float:
            alle[:, _i] = hps(h[_i], p, fluid, composition, option, 
                              units, props)
        else:
            alle[:, _i] = hps(h[_i], p[_i], fluid, composition, option, 
                              units, props)
    return alle

def tp(temp, p,  fluid, composition=[1.0], option=1, units =_units, 
       props=_props):
    """
    Properties needed for integration at given p and h, single phase.

    Parameters
    ----------
   
    p : float
        pressure in Pa.
    temp : float
         temperature in K.
    
    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        if option = 0 it includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.
        
        if option =1: it is compatible to the high level output
        T,p,h,v,s,x.

    """
    if props == "REFPROP":
        o = RP.REFPROP2dll(fluid,"TP","H;D;S;q", units, 0, temp, p, composition)
        if option == 0:
            alle =[]
            
        elif option == 1:
            alle =[temp, p, o.Output[0], 1 / o.Output[1], *o.Output[2:4]]
    elif props =="CoolProp":
       fluid.update(CP.PT_INPUTS, p, temp)
       reihe = [CP.iHmass, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                CP.iconductivity, CP.iCpmass, CP.iviscosity]
       props = [fluid.keyed_output(k) for k in reihe]
       h, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]
       
       if option == 0:
           alle = np.array([temp, p, x, h,  s, rho, mu,
                   cp, lambda_s, phase, prandtl])
       elif option == 1:
           return [temp, p,  h, 1/rho, s,  x]
    
    return np.array(alle)
#  below must be checked!


def ht_properties_sat(p, fluid):
    """
    Properties needed for integration at given p and h for heat transfer.

    Parameters
    ----------
    p : float
        pressure in Pa.

    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array (2,4)
        includes: tranport properies in saturated state at given pressure p
         of liquid (0,:) and vapor(1,:),
        densities,viscosities, cp, thermal conductivity

        all in SI units.

    """

    props_all = np.zeros((2,4))
    for phase in [0,1]:
        fluid.update(CP.PQ_INPUTS, p, phase)
        reihe = [CP.iDmass, CP.iCpmass, CP.iviscosity, CP.iconductivity]
        props = [fluid.keyed_output(k) for k in reihe]
        props_all[phase,:] = props[:]

    return props_all


def ht_properties_satV(p, h, fluid): # unbenutzte vektorisierung
    _n = len(p)
    alle = np.zeros((14, 2, _n))
    for _i in range(_n):
        alle[:, :, _i] = ht_properties_sat(p[_i], h[_i], fluid)
    return alle


def properties_V(p, h, fluid, option=1):
    """ Vectorization of the single phase properties function"""
    _n = len(p)
    alle = np.zeros((11, _n))
    for _i in range(_n):
        alle[:, _i] = hps(h[_i], p[_i], fluid, option=option)
    return alle

if __name__ == "__main__":
        
    # _vielPrint__ = False
    
    _props ="CoolProp"
    x_0 = 1.
    p_0 = 20e5     # Anfangsdruck Pa
    
    
    temp_sur = 283.15
    temp_0_s = 373.15
    
    
    p_sur = 1.013e5
    if _props =="CoolProp":
        fluid_a = "n-Propane"  # Working fluid
        fluid_s = "Water"
        p_c = CP.PropsSI('Pcrit', fluid_a)
        temp_0 = CP.PropsSI('T', 'P', p_0, 'Q', x_0, fluid_a) 
        
        h_0 = CP.PropsSI('H', 'P', p_0, 'T', temp_0 + 1, fluid_a)
        working_fluid = CP.AbstractState("BICUBIC&HEOS", fluid_a)
        mm = ht_properties_sat(1e6, working_fluid)
        secondary_fluid = CP.AbstractState("TTSE&HEOS", fluid_s) 
        h_0_s = tp(temp_0_s, p_sur,  secondary_fluid, props=_props)[2]
        h_end = CP.PropsSI('H', 'P', p_sur, 'T', temp_0_s, fluid_a)
    elif _props == "REFPROP":
        # Sekundärfluid --------------------------------
        fluid_s = "Propane * Pentane"
        comp =[.5, 0.5]
        #secondary_fluid = CP.AbstractState("TTSE&HEOS", fluid_s) 
        # interesting, when using "BICUBIC&HEOS" the exergy of the ambient state is 0.15!
        
        state_data = tp(temp_0_s, p_sur, fluid_s, composition=comp)
        print(state_data)
        print(hps(state_data[2],p_sur, fluid_s, composition=comp))
    h_0_s = tp(temp_0_s, p_sur,  secondary_fluid, props=_props)[2]
    h_end = CP.PropsSI('H', 'P', p_sur, 'T', temp_0_s, fluid_a)
    ex1 = hp_exergy(h_0_s, p_sur, secondary_fluid, props=_props)
    ex2 = hp_exergy(h_0, p_sur, working_fluid, props=_props)
    print( "Exergies (J/kg):", ex1, ex2)
    
