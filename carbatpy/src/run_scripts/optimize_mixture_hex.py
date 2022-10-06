# -*- coding: utf-8 -*-
"""
Example: Optimize mixture composition for minimum heat transfer to a single
storage fluid running in counterflow
Created on Thu Oct  6 15:17:49 2022

@author: atakan
"""

import numpy as np
import CoolProp.CoolProp as CP
from fluid_properties_rp import tp, hps, hps_v, hp_exergy
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
props = "REFPROP"


import heat_exchanger as HE

from scipy.optimize import minimize



def mixdep_hex_entropy(x, argumente):
    """ minimize the entropy generation in the heat exchanger by varying the
    composition of one fluid (constraint: mole fractions sum to 1 and are 
    positive). The low of the other fluid is given and its initial state.
    For the fluid to optimize T is given.
    Mass flow rate and diameter/ length could be included to the optimization.
    """
    mdot1, p1, x1, x2= x
    
    props, mdot, alpha, Tin, p, diameters, length, compositions, fl, ha_in, hb_in = argumente
    mdot[0]=mdot1
    compositions[0] = [x1, x2]
    p[0] = p1
    # props = "REFPROP"  # "CoolProp"  # "REFPROP"
    # mdot=np.array((mdot1, .025)) # kg/s for both fluids
    # alpha = 500  # heat transfer coefficient through the wall (from total resistance)
    # Tin = [354, 309]  # initial fluid temperatures, assuming single phae each!
    # p = [p1, 4.e5]  # pressure for each fluid, Pa
    # diameters =[1.5e-2, 5e-2]  # m
    # length = 3.  # m
    
    # fl1 ="Isobutane * Pentane"
    # fl2 = "Water"
    # compositions =[[x1, x2],[ 1.]]
    # fl =[fl1,fl2]
    
    
    
    # #  evaluate enthalpies and maximum possible enthalpy changes:
    # ha_in = tp(Tin[0], p[0], fl[0], props=props, 
    #            composition=compositions[0])[2]  # state of fluid 1 left
    # hb_in=tp(Tin[1],p[1],fl[1], props=props, 
    #          composition =compositions[1])[2]  # state of fluide 2 right (at L)
    
    
    heat_ex = HE.counterflow_hex(fl, mdot, p, [ha_in,hb_in], 
                              length, diameters, U=alpha, no_tubes=2, 
                              props=props, compositions =compositions)  # assign parameters
    res =heat_ex.he_bvp_solve()  # solve the heat exchanger problem
    if res.success: 
        f1, f2, entropy_flowrate, dq =heat_ex.he_state(res, 1) # evaluate results (and plot)
    else:  entropy_flowrate = 1e9
    print(x, -dq / entropy_flowrate)
    return -dq / entropy_flowrate



if __name__ == "__main__":
    
    props = "REFPROP"
    mdot=np.array((.02, .025))
    alpha = 500  
    Tin = [354, 309] 
    p = [1e6, 4.e5] 
    diameters =[1.5e-2, 5e-2]  
    length = 3.
    fl1 ="Propane * Pentane"
    fl2 = "Water"
    compositions =[[.5, .5],[ 1.]]
    fl =[fl1,fl2]
    ha_in = tp(Tin[0], p[0], fl[0], props=props, 
               composition=compositions[0])[2]  
    hb_in=tp(Tin[1],p[1],fl[1], props=props, 
             composition =compositions[1])[2]
    ha_outMax = tp(Tin[1],p[0],fl[0], props=props, 
                   composition=compositions[0])[2]  # state of fluid 1 left
    hb_outMax = tp(Tin[0],p[1],fl[1], props=props, 
                   composition =compositions[1])[2]  # state of fluide 2 righ
    mdot[0] = mdot[1] * (hb_outMax -  hb_in)/(ha_in -ha_outMax)
    argumente =[props, mdot,
    alpha, Tin  , p  , diameters , length,
    compositions,
    fl, ha_in, hb_in
    ]
    
    cons = ({'type': 'eq', 'fun': lambda x:  x[2] + x[3] - 1})
    bnds =((mdot[0] / 2, mdot[0] * 2), (1e5, 2e6),(0,1),(0,1))
    x0 = [mdot[0]*1.05, 5e5, 0.5, 0.5]  # initial guess
    result = minimize(mixdep_hex_entropy, x0, args =argumente, 
                      method ='SLSQP', bounds =bnds, 
                       constraints=cons)
    print(result)
    