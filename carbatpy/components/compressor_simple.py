# -*- coding: utf-8 -*-
"""
functions for compressor and expander output state calcultations

so far/here: only for fixed isentropic efficiencies

((Part of carbatpy.))
Created on Sun May 21 08:51:33 2023

@author: atakan
"""


import CoolProp.CoolProp as CP

import fluid_properties_rp as fprop
import numpy as np


def compressor(state_in, p_out, eta_c, working_fluid, props="REFPROP",
               composition=[1.0], calc_type="const_eta",
               name="compressor",
               units=21, WF=fprop.RP):
    """
    compressor output state calculation
    
    so far only for a constant isentropic efficiency

    Parameters
    ----------
    state_in : array of float
        state containing [T,p,h,v,s,q].
    p_out : float
        output pressure.
    eta_c : float
        isentropic efficiency.
    working_fluid : string
        fluid (mixture according to Refprop.
    props : string, optional
        property model. The default is "REFPROP".
    composition : array of floats, optional
        mle fractions of each component. The default is [1.0].
    calc_type : string, optional
        how to calculate, so far, only one implemented. The default is "const_eta".
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
    if calc_type == "const_eta":
        prop_input = [working_fluid, composition , 1, units, props, WF]
        state_out_isentropic = fprop.sp( state_in[4], p_out, *prop_input)  
        enthalpy_change_real = (state_out_isentropic[2] - state_in[2]) / eta_c
        state_out = fprop.hp(state_out_isentropic[2] + enthalpy_change_real,
                             p_out, *prop_input)
    else:
        raise Exception(f"The option{calc_type} is not yet implemented for compressors")
    return state_out


def expander(state_in, p_out, eta_s, working_fluid, props="REFPROP",
               composition=[1.0], calc_type="const_eta",
               name="compressor",
               units=21, WF=fprop.RP):
    """
    expander output state calculation
    
    so far only for a constant isentropic efficiency

    Parameters
    ----------
    state_in : array of float
        state containing [T,p,h,v,s,q].
    p_out : float
        output pressure.
    eta_s : float
        isentropic efficiency.
    working_fluid : string
        fluid (mixture according to Refprop.
    props : string, optional
        property model. The default is "REFPROP".
    composition : array of floats, optional
        mle fractions of each component. The default is [1.0].
    calc_type : string, optional
        how to calculate, so far, only one implemented. The default is "const_eta".
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
    if calc_type == "const_eta":
        prop_input = [working_fluid, composition , 1, units, props, WF]
        state_out_isentropic = fprop.sp( state_in[4], p_out, *prop_input)  
        enthalpy_change_real = (state_out_isentropic[2] - state_in[2]) * eta_s
        state_out = fprop.hp(state_out_isentropic[2] + enthalpy_change_real,
                             p_out, *prop_input)
    else:
        raise Exception(f"The option{calc_type} is not yet implemented for compressors")
    return state_out

if __name__ == "__main__":
    

    fluid = "Propane * Butane"
    composition =[0.8,0.2]
    #wf = fprop.setRPFluid(fluid, fprop.REFPROPFunctionLibrary)
    p_in = 1e5
    T_in = 300.
    p_out = 10e5
    eta_c = .65
    state_in = fprop.tp(T_in, p_in, fluid, composition)
    state_ina = fprop.sp(state_in[4], p_out, fluid, composition)
    state_out = compressor(state_in, p_out, eta_c, fluid, composition=composition)
    print(state_in,"\n", state_out,"\n", state_out-state_in)
