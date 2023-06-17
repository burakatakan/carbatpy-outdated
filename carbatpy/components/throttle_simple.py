# -*- coding: utf-8 -*-
"""
Created on Sun May 21 08:51:33 2023

@author: atakan
"""


import CoolProp.CoolProp as CP

import fluid_props as fprop
import numpy as np


def throttle(state_in, p_out, fluid, calc_type="const_h",
               name="throttle"):
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
        state_out = fluid.set_state([fluid.properties.enthalpy, p_out],"HP")
    else:
        raise Exception(f"The option{calc_type} is not yet implemented for compressors")
    return state_out


if __name__ == "__main__":
    FLUID = "Propane * Pentane"
    comp = [.50, 0.5]
    flm = fprop.FluidModel(FLUID)
    myFluid = fprop.Fluid(flm, comp)
    state_in = myFluid.set_state([320., 19e5], "TP")

    p_out = 5e5
    

    state_out = throttle(state_in, p_out, myFluid)
    print("Throttle:\nInput:", state_in,"\nOutput:", state_out,"\nDifference", state_out-state_in)
