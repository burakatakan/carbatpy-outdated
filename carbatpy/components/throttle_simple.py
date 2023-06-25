# -*- coding: utf-8 -*-
"""
Created on Sun May 21 08:51:33 2023

@author: atakan
"""


import CoolProp.CoolProp as CP

import fluid_properties_rp as fprop
import numpy as np


def throttle(state_in, p_out, working_fluid, props="REFPROP",
               composition=[1.0], calc_type="const_h",
               name="compressor",
               units=21, WF=fprop.RP):
    """
    throttle output state calculation
    
    so far only for a constant enthalpy

    Parameters
    ----------
    state_in : array of float
        state containing [T,p,h,v,s,q].
    p_out : float
        output pressure.
    working_fluid : string
        fluid (mixture according to Refprop.
    props : string, optional
        property model. The default is "REFPROP".
    composition : array of floats, optional
        mle fractions of each component. The default is [1.0].
    calc_type : string, optional
        how to calculate, so far, only one implemented. The default is "const_h".
    name : string, optional
        name of the device. The default is "compressor".
    units : integer, optional
        to select SI (here or REFPROP). The default is 21.
    WF : instnce of refprop, optional
        DESCRIPTION. The default is fprop.RP.

    Returns
    -------
    state_out : array of float
        compressor output state containing [T,p,h,v,s,q].

    """
    if calc_type == "const_h":
        prop_input = [working_fluid, composition , 1, units, props, WF]
        
        state_out = fprop.hp(state_in[2] ,
                             p_out, *prop_input)
    else:
        raise Exception(f"The option{calc_type} is not yet implemented for compressors")
    return state_out


if __name__ == "__main__":
    

    fluid = "Propane * Butane"
    composition =[0.8,0.2]
    #wf = fprop.setRPFluid(fluid, fprop.REFPROPFunctionLibrary)
    p_in = 19e5
    T_in = 300.
    p_out = 10e5
    
    state_in = fprop.tp(T_in, p_in, fluid,composition)
    state_ina = fprop.sp(state_in[4], p_out, fluid,composition)
    state_out = throttle(state_in, p_out, fluid, composition=composition)
    print(state_in,"\n", state_out,"\n", state_out-state_in)
