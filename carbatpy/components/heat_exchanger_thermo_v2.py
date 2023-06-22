# -*- coding: utf-8 -*-
"""
Created on Sun May 21 08:51:33 2023

@author: atakan
"""


import CoolProp.CoolProp as CP
import copy
import fluid_props as fprop
import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt


class static_heat_exchanger:
    """
    Class for static counter-flow heat exchanger

    means: no time dependence and not heat transfer coefficients * areas
    are used!
    Only first law and second law will be checked (the latter must be improved)

    """

    def __init__(self, fluids, dH_min, T_out_s,
                 points=30, dT_separation_min=0.5, calc_type="const",
                 name="evaporator"):
        """
        class to calculate (static/steady state) heat-exchangers
        
        includes pinch-point analysis and plotting, 
        only implemented for simple thermodynamic calculations 
        (no convection coefficients and heat exchanger areas regarded yet)
        
        Parameters
        ----------
        fluids : list of 2 fprop.Fluid
            The definition of the two fluids, as they enter the heat exchanger.
        dH_min : float
            enthalpy flow rate (W) which has to be transfered.
        T_out_s : float
            exit temperature of the secondary fluid.
        points : integer, optional
            for how many points (array) shall the minimum approach temperature
            be checked and properties be returned (for plotting etc.)
            The default is 30.
        dT_separation_min : float, optional
            Minimium approach temperature (pinch point) between the two fluids. 
            The default is 0.5.
        calc_type : string, optional
            which calculation type shall be performed; only one implemented so
            far. The default is "const".
        name : string, optional
            name of the heat exchanger. The default is "evaporator".

        Returns
        -------
        None.

        """
        self.fluids = fluids
        self.dQ = dH_min
        self.m_dot_s = 0
        self.T_out_s = T_out_s
        self.points = points
        self. dT_separation_min = dT_separation_min
        self.heating = False
        if T_out_s > fluids[1].properties.temperature:
            self.heating = True # condenser
        self.calc_type = calc_type
        self.name = name
        self.all_states = np.zeros((points, len(fprop._THERMO_STRING)))
        self.h_in_out = np.zeros((2, 4))
        
        
        
    def pinch_calc(self, m_dot_w_factor=1):
        """
        Calculate the changes in tnthalpy and temperature along the heat exchanger
        
        counter-flow hex assumed! 
        Is used to check, whether the second low is violated. The factor can 
        be used to vary the mass flow rate of the working fluid, until no 
        violation is found (done in root finding).

        Parameters
        ----------
        m_dot_w_factor : float, optional
            factor to divide the working fluid mass flow rate, as derived from 
            the enrgy balance, in order to avoid a crossing of temperature 
            curves. The default is 1.

        Raises
        ------
        Exception
            if temperatures are not consistent.

        Returns
        -------
        m_dot_w : float
            working fluid mass flow rate (kg/s.
        dTall : numpy-array
            The temperature differences along the counter-flow heat exchanger.
        w_array : array
            properties of the working fluid along the heat exchanger 
            (T,p,h, etc. see fluid class).
        s_array : array
            properties of the secondary fluid along the heat exchanger 
            (T,p,h, etc. see fluid class).

        """
        w_in = self.fluids[0]
        s_in = self.fluids[1]

        w_out = copy.copy(self.fluids[0])
        s_out = copy.copy(self.fluids[1])

        self. h_in_out[1, 0] = h_in_s = s_out.properties.enthalpy

        # fixed output temperature, secondary fluid
        s_out.set_state([self.T_out_s, s_in.properties.pressure], "TP")
        if self.heating: 
            condenser = 1
        else: 
                condenser = -1
                
        self. h_in_out[1, 1] = dh_s = (
            s_out.properties.enthalpy - s_in.properties.enthalpy) * condenser
        self.m_dot_s = self.dQ / dh_s  # fixed heat flow, determines mass flow rate
        
        
        outw = w_out.set_state([s_in.properties.temperature + 
                                self.dT_separation_min * condenser,
                                w_in.properties.pressure], "TP")
        print(outw, w_in.properties.temperature, w_out.properties.temperature)
        self. h_in_out[0, 0] = w_out.properties.enthalpy
        
        self. h_in_out[0, 1] = dh_w = (w_in.properties.enthalpy - w_out.properties.enthalpy) / \
            m_dot_w_factor * condenser # absolute value
        m_dot_w = np.abs(self.dQ / dh_w)
        

        # check pinch, first secondary fluid
        h_array = np.linspace(h_in_s, h_in_s + dh_s * condenser, self.points)
        values = np.zeros((self.points, 2))
        values[:, 0] = h_array
        values[:, 1] = s_out.properties.pressure

        s_array = s_out.set_state_v(values, given="HP")
        # print(values, m_dot_w, self.m_dot_s)

        # now working fluid:
        h_array = np.linspace(self.h_in_out[0, 0], 
                              self.h_in_out[0, 0] + self.h_in_out[0, 1] * condenser,  
                              self.points)
        values = np.zeros((self.points, 2))
        values[:, 0] = h_array.T
        values[:, 1] = w_out.properties.pressure

        w_array = w_out.set_state_v(values, "HP")
        dTall = w_array[:, 0]-s_array[:, 0]
        print(h_in_s, h_in_s + dh_s * condenser, self.h_in_out[0, 0], 
                              self.h_in_out[0, 0] + self.h_in_out[0, 1] * condenser, condenser)
        # print(w_out.properties.temperature, w_in.properties.temperature, dh_w, m_dot_w, condenser)
        return m_dot_w, dTall, w_array, s_array
    

    def pinch_root(self, factor):
        """
        function for root-finding of the minimum approach temperature
        
        a factor to divide the working fluid mass flow rate is varied. input
        for root

        Parameters
        ----------
        factor : float
            as said above.

        Returns
        -------
        float
            root tries to reach a value of 0.

        """
        mw, dTs, w, s = self.pinch_calc(factor)
        if self.heating:
            return dTs.min() - self.dT_separation_min
        else:
            return (dTs.max() - self.dT_separation_min)

    def find_pinch(self):
        if self.heating:
            start =1. 
        else: start = 1.905
        result = root(self.pinch_root, start)
        if result.success:
            return result.x
        else:
            print("root-finding problem!", result)
            return result.x

    def pinch_plot(self, factor):
        m_dot_w, dTall, w_array, s_array = self.pinch_calc(factor)
        fig, ax = plt.subplots(1, 1)
        h_w_plot_array = (
            w_array[:, 2] - self.h_in_out[0, 0]) * m_dot_w / self.m_dot_s
        ax.plot(s_array[:, 2] - self.h_in_out[1, 0], s_array[:, 0], "v")
        ax.plot(h_w_plot_array, w_array[:, 0], "o")
        ax.set_xlabel(
            "specific enthalpy flow per kg of secondary fluid / (J/kg)")
        ax.set_ylabel("temperature / (K)")
        ax.set_title("heat exchanger, simple")
        
        return m_dot_w, dTall, w_array, s_array


if __name__ == "__main__":
    FLUID = "Propane * Pentane"
    FLS = "Water"
    comp = [.50, 0.5]
    flm = fprop.FluidModel(FLUID)
    myFluid = fprop.Fluid(flm, comp)
    
    secFlm = fprop.FluidModel(FLS)
    secFluid = fprop.Fluid(secFlm, [1.])
    
    # Condenser, secondary fluid fixes all:
    sT_in = 300.0
    sT_out = 370.0
    sp_in = 2e5
    dH_min = 1e3
    
    state_sec_in = secFluid.set_state([sT_in, sp_in], "TP")
    # working fluid
    dT_super = 5.
    T_dew = sT_out
    state_in = myFluid.set_state([T_dew, 1.], "TQ") # find minimum pressure
    
    
    T_in = T_dew + dT_super 
    state_in = myFluid.set_state([myFluid.properties.pressure,
                                  myFluid.properties.temperature + dT_super], "PT")
    dT_min = 2.0
    
    hex0 = static_heat_exchanger([myFluid, secFluid], dH_min, sT_out, 
                                 dT_separation_min=dT_min)
    mw, dTall, w, s = hex0.pinch_calc(160)
    factor0 = hex0.find_pinch()
    mw0, dTall0, w0, s0 = hex0.pinch_plot(factor0)
    
    # Evaporator:
        
    
    sT_in = 300.0
    sT_out = 276
    sp_in = 2e5
    dH_min = 1e3
    
    state_sec_in = secFluid.set_state([sT_in, sp_in], "TP")
    dT_min = 2.0
    T_in = sT_out - dT_min
    state_in = myFluid.set_state([T_in, .1], "TQ")
    print("state in", state_in)
    
    
    
    hex0 = static_heat_exchanger([myFluid, secFluid], dH_min, sT_out, 
                                 dT_separation_min=dT_min)
    mw, dTall, w, s = hex0.pinch_calc(1)
    factor = hex0.find_pinch()
    mw, dTall, w, s = hex0.pinch_plot(1)

    p_out = 5e5

    # state_out = throttle(p_out, myFluid)
    # print("Throttle:\nInput:", state_in,"\nOutput:", state_out,"\nDifference", state_out-state_in)
