# -*- coding: utf-8 -*-
"""
calculate a few characteristic numbers for a heat pump 

running with a multicomponent mixture (5) and having 10 objectives, which are 
calculated. To be usd in combination with a pymoo-Optimization. Towards 
Carnot-battery fluid pre-selection.

Created on Thu Jun 22 17:40:55 2023

@author: atakan
"""

import fluid_props as fprop
import numpy as np


v_names = ["delt_ex_destruct", 
           "delt_ex_destruct_compr",
           "dT_high",  
           "T_liq(p_low)",
           "T_after compr irr",
           "T_glide(p_high)",
           "quality_subcooled_throttle", 
           "p_low", 
           "p_high", 
           "p_ratio"]   # , "$\Delta T$", "$UA_h$", "$UA_l$"]
OBJECTIVES = len(v_names)
MIN_MAX = np.array([-1, -1, -1, -1, -1,1, -1,  1, -1, -1]) # which objective to minimize(-1)

def dewpoints(input_all, fluid, fixed_points):
    """
    calculate a few characteristic numbers for a heat pump 

    running with a multicomponent mixture (5) and having 10 objectives, which are 
    calculated. To be usd in combination with a pymoo-Optimization.Pressures
    are calculated according to the temperatures.

    Parameters
    ----------
    input_all : array length 4
        mole fractions of all but one component,sum must be <=1.
    fluid : an fprop FLUID
        a fluid mixture including property model (Refprop).
    fixed_points : dictinary
        with characteristic temperatures of the heat pump and isentropic
        efficiency.

    Returns
    -------
    an array with the objectives, e.g.: Exergy-destruction, tiotal and in the
    throttle, pressure levels, pressure ratio, T-glide, lowest temperature.
    see v_names

    """

# two test cases condenser and evaporator:
    x1, x2, x3, x4 = input_all
    x5 = 1 - x1 - x2 -x3 - x4
    comp = np.array([x1, x2, x3, x4, x5])
    if (x5 < 0) or (comp.sum() > 1.0000000001):
        # print(comp)
        return np.ones((OBJECTIVES))*1e9
    
    try:
        myFluid = fluid
        myFluid.set_composition(comp)
        # print(comp, input_all)
        state_hvap = myFluid.set_state([fixed_points["T_hh"], 1.], "TQ")  # find high pressure
    
        state_lvap = myFluid.set_state([fixed_points["T_lh"], 1.], "TQ")  # find minimum pressure, dew
        p_low = state_lvap[1]
        p_high = state_hvap[1]
        diff1 =state_hvap - state_lvap
        ds_rel = diff1[4]
        
        state_isentropic = myFluid.set_state([state_lvap[4], p_high], "SP")
        # print(state_isentropic)
        w_s = state_isentropic[2]-state_lvap[2]
        # print(w_s)
        state_irr = myFluid.set_state([w_s / fixed_points["eta_s"] + 
                                       state_lvap[2], p_high], "HP")
        d_exergy_rel = fixed_points["T_hl"] * (state_irr[4] - 
                                               state_lvap[4]) / w_s
        state_lliq = myFluid.set_state([p_low, 0.], "PQ")  # find minimum pressure, boil
        state_hliq = myFluid.set_state([p_high,0], "PQ") 
        # diff2 =state_hliq - state_lliq
        throttle_liq = myFluid.set_state([p_low, state_hliq[2]], "PH")
        
        state_hl = myFluid.set_state([fixed_points["T_hl"], p_high], "TP")
        
        ds_rel = diff1[4]/(state_lvap[4] - state_lliq[4])
        dT_high = state_hliq[0] - fixed_points["T_hl"]
        throttle_subcooled_liq = myFluid.set_state([p_low, state_hl[2]], "PH")
        q_throttle = throttle_liq[5]
        q_throttle_sub = throttle_subcooled_liq[5]
        throttle_rev = myFluid.set_state([p_low, state_hl[4]], "PS")
        d_exergy_throttle_rel = np.abs((throttle_rev[2] - state_hl[2])/w_s)
        d_exergy = d_exergy_rel+d_exergy_throttle_rel
        p_ratio =p_high / p_low
        delta_low_p =state_lvap - state_lliq
        delta_high_p =state_hvap - state_hliq
        return np.array((d_exergy, d_exergy_rel, dT_high, state_lliq[0], 
                         state_irr[0], delta_high_p[0],
                          q_throttle_sub,p_low, p_high, p_ratio))
    
    except:
        print(input_all, "Problem")
        return np.ones((OBJECTIVES)) * 1e9
    
    
if __name__ == "__main__":
    fixed_points = {
                    "p_low": 1e5,
                    "T_hh": 370,
                    "T_hl": 290,
                    "T_lh": 290,
                    "T_ll": 276.0,
                    "eta_s": 0.65
                    }
    FLUID = "Propane * Butane * Pentane * Hexane * Isobutane"
    FLS = "Methanol"  # "Water"  #
    comp = [.20, 0.1, .4, .2,.1]
    flm = fprop.FluidModel(FLUID)
    myFluid = fprop.Fluid(flm, comp)
    new_in =[.3,.3,.1,.1]
    print (dewpoints(new_in, myFluid, fixed_points))