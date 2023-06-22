# -*- coding: utf-8 -*-
"""
Created on Mon May 22 21:09:12 2023

@author: atakan
"""

import os
# from time import time
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import numpy as np
# import CoolProp.CoolProp as CP


os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
os.environ['RPPREFIXs'] = r'C:/Program Files (x86)/REFPROP/secondCopyREFPROP'
_PROPS = "REFPROP"  # or "CoolProp"
_fl_properties_names = ("Temperature", "Pressure", "spec. Enthalpy",
                        "spec. Volume", "spec. Entropy", "quality",
                        "spec. internal Energy",
                        "viscosity", "thermal conductivity",
                        "Prandtl number", "k.viscosity", "molecular mass")
_THERMO_STRING = "T;P;H;V;S;QMASS;E"
_TRANS_STRING = _THERMO_STRING + ";VIS;TCX;PRANDTL;KV;M"
# order for coolprop,alle_0:[_temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, prandtl, phase]"
if _PROPS == "REFPROP":
    try:
        rp_instance = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        # be careful pressure is in Pa!
        _UNITS = rp_instance.GETENUMdll(0, "MASS BASE SI").iEnum
    except:
        print("print refprop not installed!")

else:
    _UNITS = 21


class FluidModel:

    def __init__(self, fluid, units=_UNITS, props=_PROPS, rp_inst=rp_instance):

        self.fluid = fluid
        self.props = props
        self.units = units
        self.rp_instance = rp_inst
        self.set_rp_fluid()

    def set_rp_fluid(self, modwf=REFPROPFunctionLibrary, name='RPPREFIX'):
        """
        A new instance of Refpropdll for the given fluid. It can then be called 
        using fluid =""

        Parameters
        ----------
        fluid : string
            fluid (mixture) name, as described in REFPROP.

        Returns
        -------
        self.rp_instance : REFPROP Instance
            for further usage.

        """

        self.rp_instance = modwf(os.environ[name])
        self.rp_instance.SETPATHdll(os.environ[name])
        ierr = self.rp_instance.SETFLUIDSdll(self.fluid)
        if ierr != 0:
            print(f"Fehler in setfluid {ierr}")
            print(self.rp_instance.ERRMSGdll(ierr))
        return self.rp_instance


class FluidState:
    def __init__(self, state):
        self.temperature = state.Output[0]
        self.pressure = state.Output[1]
        self.sp_volume = state.Output[3]
        self.enthalpy = state.Output[2]
        self.entropy = state.Output[4]
        self.quality = state.Output[5]
        self.int_energy = state.Output[6]
        self.state = state.Output[:7]
        self.prop_names = _fl_properties_names


class FluidStateTransport(FluidState):
    def __init__(self, state):
        super().__init__(state)
        # ";VIS;TCX;PRANDTL;KV;M"
        self.viscosity = state.Output[7]
        self.thermal_conductivity = state.Output[8]
        self.prandtl = state.Output[9]
        self.kin_viscosity = state.Output[10]
        self.molecular_mass = state.Output[11]
        self.transport = state.Output[7:12]


class Fluid:

    def __init__(self, fluidmodel, composition=[1.0],
                 option=1):
        self.fluidmodel = fluidmodel
        self.composition = composition
        self.option = option
        self.no_compounds = len(composition)

    def set_state(self, values, given="TP", wanted=_THERMO_STRING):
        if self.fluidmodel.props == "REFPROP":
            state = self.fluidmodel.rp_instance.REFPROP2dll(
                self.fluidmodel.fluid, given, wanted,
                self.fluidmodel.units,
                0, values[0], values[1],
                self.composition)

            if state.ierr == 0:
                if wanted == _THERMO_STRING:
                    # print(state)
                    self.properties = FluidState(state)
                elif wanted == _TRANS_STRING:
                    self.properties = FluidStateTransport(state)
                else:

                    raise Exception(f"properties{wanted} not implemented yet!")
            else:
                self.herr = state.herr
                raise Exception(f"Property-Refprop problem: {state.herr}!")
        else:
            raise Exception(
                f"Property model {self.fluidmodel.props} not implemented yet!")
        return np.array([*self.properties.state])

    def set_state_v(self, values, given="TP", wanted=_THERMO_STRING):
        dimension = np.shape(values)
        number_wanted = wanted.count(";")+1
        output = np.zeros((dimension[0], number_wanted))
        for count, value in enumerate(values):
            output[count, :] = self.set_state(value, given, wanted)
        self.state_v = output
        return output
    
    def print_state(self):
        pr =self.properties
        flm = self.fluidmodel
        print(f"\n{flm.fluid}, composition: {self.composition}")
        print(f"T:{pr.temperature:.2f} K, p: {pr.pressure/1e5 :.2f} bar,  h: {pr.enthalpy/1000: .2f} kJ/kg, s: {pr.entropy/1000:.3f} kJ/kg K\n")
        if pr.quality>=0:
            print(f"Quality: {pr.quality :.3f}")
            


if __name__ == "__main__":
    FLUID = "Propane * Pentane"
    comp = [.50, 0.5]
    flm = FluidModel(FLUID)
    myFluid = Fluid(flm, comp)
    st0 = myFluid.set_state([300., 1e5], "TP")
    st1 = myFluid.set_state([300., 1e5], "TP",
                            _TRANS_STRING)
    print(st0, st1)
    myFluid.print_state()
    # value_vec = np.array([[300, 1e5], [400, 1e5], [500, 1e5]])
    # stv = myFluid.set_state_v(value_vec, "TP")

    # print(myFluid.set_state_v(value_vec, "TP"))
    # print(myFluid.set_state([300., 1.], "TQ"))
