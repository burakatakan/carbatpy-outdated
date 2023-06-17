# -*- coding: utf-8 -*-
"""
functions for compressor and expander output state calcultations

so far/here: only for fixed isentropic efficiencies

((Part of carbatpy.))
Created on Sun May 21 08:51:33 2023

@author: atakan
"""


import CoolProp.CoolProp as CP

import fluid_props as fprop
import numpy as np


def compressor(p_out, eta_s, fluid, calc_type="const_eta",
               name="compressor"):
    """
    compressor or expander output state calculation

    so far only for a constant isentropic efficiency, according to the pressure
    change an expansion or compression is detected and handled.

    Parameters
    ----------
    state_in : array of float
        state containing [T,p,h,v,s,q].
    p_out : float
        output pressure.
    eta_s : float
        isentropic efficiency.
    fluid : fprop.Fluid
        entering fluid, including properties, composition, and model.
    calc_type : string, optional
        how to calculate, so far, only one implemented. The default is "const_eta".
    name : string, optional
        name of the device. The default is "compressor".

    Returns
    -------
    state_out : array of float
        compressor output state containing [T,p,h,v,s,q].

    """
    expander = False
    if fluid.properties.pressure > p_out:
        expander = True

    if calc_type == "const_eta":
        h_in = fluid.properties.enthalpy

        state_out_isentropic = fluid.set_state(
            [fluid.properties.entropy, p_out], "SP")

        diff_enthalpy_s = fluid.properties.enthalpy-h_in

        if expander:
            diff_enthalpy = diff_enthalpy_s * eta_s
        else:
            diff_enthalpy = diff_enthalpy_s / eta_s

        state_out = fluid.set_state([h_in + diff_enthalpy, p_out], "HP")
    else:
        raise Exception(
            f"The option{calc_type} is not yet implemented for compressors")
    return state_out


if __name__ == "__main__":

    FLUID = "Propane * Pentane"
    comp = [.80, 0.2]
    flm = fprop.FluidModel(FLUID)
    myFluid = fprop.Fluid(flm, comp)
    p_low = 1e5
    state_in = myFluid.set_state([310., p_low], "TP")
    p_out = 10e5
    eta_s = .7

    state_out = compressor(p_out, eta_s, myFluid)
    print(myFluid.properties.temperature)
    print("\nCompressor", state_in, "\n", state_out, "\n", state_out-state_in)
    state_in = state_out
    state_out = compressor(p_low, eta_s, myFluid)
    print("\nExpander:", state_in, "\n", state_out, "\n", state_out-state_in)
