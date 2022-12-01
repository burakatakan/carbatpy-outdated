# -*- coding: utf-8 -*-
"""
Simple calculation of heat pump states for multinerary mixtures (using Refprop)
mainly for visualization of temperature and enthalpy levels for different 
pressure levels and slight superheatig after evaporation.

Created on Tue Nov 22 14:55:50 2022

@author: atakan
"""

import fluid_properties_rp as fprop

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opti

_props = "REFPROP"
Tevap = 273.15 - 20
Tcond = 273.15+60
Tin = 273.15 + 20
Tout = 273.15 + 131.5

fluid_s = "Propane * Pentane * Heptane"
x0=.695
x1 = 0.304
comp = [x0, x1, 1-x0-x1]# [.4, 0.5 , 0.1]

def hp_calc(pev, pcond, comp =[1.0], dT_cin=10, fluid_s =fluid_s, eta_s =0.65):
    """
    simple heat pump calculation (also for plotting the state changes)
    

    Parameters
    ----------
    pev : float
        evaporator pressure in Pa.
    pcond : float
        condenser pressure in Pa.
    comp : list of floats, optional
        molar composition, must sum to 1. The default is [1.0].
    dT_cin : float, optional
        superheating at compresor entrance. The default is 10.
    fluid_s : depends on property module(REFPROP/CoolProp), optional
        working fluid compounds (list). The default is fluid_s.
    eta_s : float, optional
        isentropic efficiency compressor. The default is 0.65.

    Returns
    -------
    numpy array (npo+5, 6)
        states along the cyle, for each T,p,h,v,s,q.

    """
   
    evap = fprop.p_prop_sat(pev, fluid_s, composition = comp, option=1) # evaporator pressure
    comp_in = fprop.tp(evap[0, 0]+ dT_cin, evap[0, 1], fluid_s, composition=comp, option = 1)  # compressor entrance at 20C
    cond = fprop.p_prop_sat(pcond, fluid_s, composition = comp, option=1) # condenser pressure
    npo =10
    con =[]
    for h in np.linspace(cond[0,2],cond[1,2],npo):
        con.append(fprop.hp(h, cond[0,1],fluid_s, composition = comp, option=1))
    # isentropic state after compressor
    comp_s = fprop.sp(comp_in[ 4], cond[0, 1], fluid_s, composition=comp) 
    dhs = comp_s[2] - comp_in[2]
    
    dh = dhs/eta_s
    comp_out = fprop.hp(comp_in[2]+dh, cond[0,1], fluid_s, composition=comp) 
    
    
    return np.array([evap[0,:], comp_in, comp_out, *con,
                     evap[1],evap[0,:]])

def find_comp(x, x0, p=1.1e5, T=289.):
    """
    find the ternary composition for given compressor inlet T & p
    function for scipy.root

    Parameters
    ----------
    x : float (0<x<1)
        mole fraction to vary, the second one in the list.
    x0 : float (0<x<1)
        fixed mole fraction , the first one in the list.
    p : float, optional
        inlet pressure in Pa. The default is 1.1e5.
    T : float, optional
        inlet pressure in K. The default is 289..

    Returns
    -------
    float
        difference between wanted T and T of the mixture at p.

    """
    
    comp = [x0, x, 1- x - x0]
    if (comp[-1] > 0) and (x>0) :
        evap = fprop.p_prop_sat(p, fluid_s, composition = comp, option=1)
        print(comp, evap[0,0], evap[0,0]-T)
        return evap[0,0] - T
    else:
        return 1e6
    

if __name__ == "__main__":
    # calculate a heat pump cycle for a given composition and pressure levels
    res = hp_calc(1.1e5,1.1e6, comp= comp)
    f, ax = plt.subplots(1)
    ax.plot(res[:,2], res[:,0])
    # find the ternary composition for given compressor inlet T & p:
    loes = opti.root(find_comp, 0.2, args=(0.695,1.5e5, 290.))
    
    