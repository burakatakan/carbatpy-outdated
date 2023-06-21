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
    Class for static heat exchanger

    means: no time dependence and not heat transfer coefficients * areas
    are used!
    Only first law and second law will be checked (the latter must be improved)

    """

    def __init__(self, fluids, dH_min, T_out_s,
                 points=30, dT_separation_min=0.5, calc_type="const",
                 name="evaporator"):
        self.fluids = fluids
        self.dQ = dH_min
        self.m_dot_s = 0
        self.T_out_s = T_out_s
        self.points = points
        self. dT_separation_min = dT_separation_min
        self.heating = False
        if T_out_s > fluids[1].properties.temperature:
            self.heating = True
        self.calc_type = calc_type
        self.name = name
        self.all_states = np.zeros((points, len(fprop._THERMO_STRING)))
        self.h_in_out = np.zeros((2, 4))

    def pinch_calc(self, m_dot_w_factor=1):
        w_in = self.fluids[0]
        s_in = self.fluids[1]

        w_out = copy.copy(self.fluids[0])
        s_out = copy.copy(self.fluids[1])

        
        self. h_in_out[1, 0] = h_in_s = s_out.properties.enthalpy
        
        

        s_out.set_state([self.T_out_s, s_in.properties.pressure], "TP")  # fixed output temperature
        self. h_in_out[1, 1] = dh_s = (s_out.properties.enthalpy - s_in.properties.enthalpy)
        self.m_dot_s = self.dQ / dh_s  # fixed heat flow, determines mass flow rate
        
        if self.heating:  # condenser
            if w_in.properties.temperature <= s_out.properties.temperature:
                raise Exception("Working fluid entering to cold!")
            outw = w_out.set_state([s_in.properties.temperature + self.dT_separation_min,
                                    w_in.properties.pressure], "TP")
            self. h_in_out[0, 0]  = w_out.properties.enthalpy
            self. h_in_out[0, 1] = dh_w = (w_in.properties.enthalpy - w_out.properties.enthalpy) / \
                    m_dot_w_factor # absolute value
            m_dot_w = self.dQ / dh_w

        # check pinch, first secondary fluid
        h_array = np.linspace(h_in_s, h_in_s + dh_s, self.points)
        values = np.zeros((self.points, 2))
        values[:, 0] = h_array
        values[:, 1] = s_out.properties.pressure

        s_array = s_out.set_state_v(values, given="HP")
        print(values, m_dot_w, self.m_dot_s)
        
        # now working fluid:
        h_array = np.linspace(self.h_in_out[0,0], self.h_in_out[0,:].sum(),  self.points)
        values = np.zeros((self.points, 2))
        print(np.shape(values), np.shape(h_array))
        values[:, 0] = h_array.T
        values[:, 1] = w_out.properties.pressure   
        
        w_array = w_out.set_state_v(values, "HP")
        dTall = w_array[:, 0]-s_array[:, 0]
        return m_dot_w, dTall, w_array, s_array
    
    def pinch_root(self, factor):
        mw, dTs, w, s = self.pinch_calc(factor)
        return dTs.min() - self.dT_separation_min
    
    def find_pinch(self):
        result = root(self.pinch_root,1.)
        if result.success:
            return result.x
        else: print("root-finding problem!")
    def pinch_plot(self, factor):
        m_dot_w, dTall, w_array, s_array = self.pinch_calc(factor)
        fig, ax =plt.subplots(1,1)
        h_w_plot_array = (w_array[:,2] - self.h_in_out[0,0]) * m_dot_w / self.m_dot_s
        ax.plot(s_array[:,2] - self.h_in_out[1, 0], s_array[:,0] ,"v") # 
        ax.plot(h_w_plot_array, w_array[:,0] ,"o") #
        fig.show()
        return m_dot_w, dTall, w_array, s_array
            
        
        


if __name__ == "__main__":
    FLUID = "Propane * Pentane"
    FLS = "Water"
    comp = [.50, 0.5]
    flm = fprop.FluidModel(FLUID)
    myFluid = fprop.Fluid(flm, comp)
    T_in = 370.
    state_in = myFluid.set_state([T_in, 1.], "TQ")
    dT_super = 5.
    dT_min = 2.0
    state_in = myFluid.set_state([myFluid.properties.pressure,
                                  myFluid.properties.temperature + dT_super], "PT")
    sT_in = 300.0
    sp_in = 2e5
    dH_min = 1e3
    sT_out = T_in - dT_super - dT_min
    secFlm = fprop.FluidModel(FLS)
    secFluid = fprop.Fluid(secFlm, [1.])
    state_sec_in = secFluid.set_state([sT_in, sp_in], "TP")
    hex0 = static_heat_exchanger([myFluid, secFluid], dH_min, sT_out)
    mw, dTall, w, s = hex0.pinch_calc(160)
    factor = hex0.find_pinch()
    mw, dTall, w, s = hex0.pinch_plot(factor)

    p_out = 5e5

    # state_out = throttle(p_out, myFluid)
    # print("Throttle:\nInput:", state_in,"\nOutput:", state_out,"\nDifference", state_out-state_in)
