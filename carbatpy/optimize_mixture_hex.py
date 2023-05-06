# -*- coding: utf-8 -*-
"""
Example: Optimize mixture composition for minimum heat transfer to a single
storage fluid running in counterflow
preparation for multiobjective optimiztion.

Warning; Must be checked work in progress!
Created on Thu Oct  6 15:17:49 2022/ Changed 06.05.2023

@author: atakan
"""

import numpy as np
#import CoolProp.CoolProp as CP
from fluid_properties_rp import tp
#from scipy.integrate import solve_bvp
#import matplotlib.pyplot as plt
props = "REFPROP"


import heat_exchanger as HE

from scipy.optimize import minimize, differential_evolution



def mixdep_hex_entropy(x, argumente):
    """ minimize the entropy generation in the heat exchanger by varying the
    composition of one fluid (constraint: mole fractions sum to 1 and are 
    positive). The low of the other fluid is given and its initial state.
    For the fluid to optimize T is given. mass flow rate of the first fluid 
    is calculated to match the mass flow rate of the second for a reversible
    heat exchanger.
    Mass flow rate and diameter/ length could be included to the optimization.
    """
    p1, x1= x
    
    props, mdot, alpha, Tin, p, diameters, length, compositions, fl, ha_in, \
    hb_in, multi_objective = argumente
    
    compositions[0] = [x1, 1-x1]
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
    
    verbose = False  # if True it will not work with multiprocessing
    
    
    # #  evaluate enthalpies and maximum possible enthalpy changes:
    a_in = tp(Tin[0], p[0], fl[0], props=props, 
               composition=compositions[0])  # state of fluid 1 left
    b_in=tp(Tin[1],p[1],fl[1], props=props, 
              composition =compositions[1])  # state of fluide 2 right (at L)
    a_out = tp(Tin[1], p[0], fl[0], props=props, 
               composition=compositions[0])  # state of fluid 1 left
    hb_out=tp(Tin[0],p[1],fl[1], props=props, 
              composition =compositions[1])[2]  # state of fluide 2 right (at L)
   
    
    ha_in = a_in[2]
    hb_in = b_in[2]
    ha_out = a_out[2]
    mdot[0]= mdot[1] *(hb_out-hb_in) / (ha_in-ha_out)
    if verbose: print(*a_in, mdot[0],*x, a_out[0])
    heat_ex = HE.counterflow_hex(fl, mdot, p, [ha_in,hb_in], 
                              length, diameters, U=alpha, no_tubes=6, 
                              props=props, compositions =compositions)  # assign parameters
    res =heat_ex.he_bvp_solve()  # solve the heat exchanger problem
    if res.success: 
        f1, f2, entropy_flowrate, dq =heat_ex.he_state(res, 1) # evaluate results (and plot)
        if (entropy_flowrate < 0) or (dq<0):
            if verbose: print("thermo failed:",x, entropy_flowrate, dq)
            if multi_objective: return 1e8, 1e8, mdot[0]
            else: return 1e8
        else:
            if verbose: print(x, f"q:{-dq:-2f} ,S:{ entropy_flowrate:.3f}, \
                              m:{mdot[0]:.3f}")
            if multi_objective:
                return -dq , entropy_flowrate, mdot[0] # changed for pymoo! 2 objectives
            else:
                return entropy_flowrate
    else:  
        if verbose: print("bvp failed:",x, res.message)
        return  np.array([1e9, 1e9, mdot[0]])
    



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
    
    mdot[0] = mdot[1] * (hb_outMax -  hb_in)/(ha_in -ha_outMax)  #select a reasonable mass flow rate of the working fluid
    argumente =[props, mdot,
    alpha, Tin  , p  , diameters , length,
    compositions,
    fl, ha_in, hb_in, False
    ]
    dq_min =3e3
    #cons = ({'type': 'eq', 'fun': lambda x:  x[2] + x[3] - 1})
    bnds =( (1e5, 2e6),(0,1))
    x0 = [ 15e5, 0.75]  # initial guess
    print("running ...")
    result = differential_evolution(lambda x: mixdep_hex_entropy(x, argumente), 
                                    bnds, workers=1 )
    # result = minimize(mixdep_hex_entropy, x0, args =argumente, 
    #                   method ='SLSQP', bounds =bnds)  #, constraints=cons)
    print(result)
    mdot[0], p[0], xa=result.x
    compositions[0] =[xa,1-xa]
    heat_ex = HE.counterflow_hex(fl, mdot, p, [ha_in,hb_in], 
                              length, diameters, U=alpha, no_tubes=2, 
                              props=props, compositions =compositions)
    res =heat_ex.he_bvp_solve()
    heat_ex.he_state(res, option =6)
    