# -*- coding: utf-8 -*-
"""
Heat exchangers
simple NTU with constant convection coefficient 
and
solving the bvp along an axial coordinate

if nothing else given: counterflow/double-pipe
if nothing else is said given: steady state

Created on Sun Oct  2 16:27:20 2022

@author: atakan
"""
import numpy as np
import CoolProp.CoolProp as CP
from fluid_properties_hl import tp, hps, xp, xT
from scipy.integrate import solve_bvp

class heat_exchanger:
    # Heat exchanger base class
    
    def __init__(self, fluids, mass_flows, pressures, enthalpies, UA=10, 
                 calc_type="const", name="HEX_0"):
        """
        

        Parameters
        ----------
        fluids : TYPE list 
            two coolProp-Fluid names.
        mass_flows : TYPE list or numpy-array length 2
           both mass flows.
        pressures : list or numpy-array length 2
            initial pressures of each fluid.
        enthalpies : list or numpy-array length 2
            initial pressures of each fluid.
        UA : float
            overall heat transfer coefficient times area, 
            as constant (for simple calculations), default = 10.
        calc_type : string
            if "const" = does not vary along heat exchanger (default)
            if "calc" = calculated from Nusselt correlation 
            along heat exchanger
        name: string.
            name of the heat-exchanger, default = "HEX_0"

        Returns
        -------
        None.

        """
        self.fluids = fluids
        self.mass_flows = mass_flows
        self.pressures = pressures
        self.enthalpies = enthalpies
        self.UA = UA
        self.calc_type = calc_type
        self.name = name
    
    def q_max(self, option=0):
        """ maximum possible heat transfer for isobaric, adiabatic 
        heat exchanger
        As it is programmed now, saturated final states are not considered"""
        
        state_variables = np.zeros((2,6))
        final_states = np.zeros((2,6))
        q_dot = np.zeros((2))
        
        for n in range(2): # get initial states (Temperatures)
            state_variables[n,:] = hps(self.enthalpies[n],self.pressures[n], 
                                       self.fluids[n])
        final_states[0,:] = tp(state_variables[1, 0],state_variables[0, 1], 
                             self.fluids[0])
        final_states[1,:] = tp(state_variables[0, 0],state_variables[1, 1], 
                             self.fluids[1])
        for n in range(2):
            q_dot[n] = self.mass_flows[n] * (final_states[n, 2] 
                                             - state_variables[n, 2])
        
        if np.abs(q_dot[0]) > np.abs(q_dot[1]):
            q_max = q_dot[1]
        else: q_max = q_dot[0]
        print (q_dot)
        if option == 0:
            return q_max
        if option == 1:
            return q_max, q_dot, final_states
            
class counterflow_hex(heat_exchanger):
    def __init__(self, fluids, mass_flows, pressures, enthalpies, length, 
                 diameters, U=10, no_tubes=1, no_points=100,
                 calc_type="const", name="HEX_0"):
        """
        

        Parameters
        ----------
        fluids : TYPE list 
            two coolProp-Fluid names.
        mass_flows : TYPE list or numpy-array length 2
           both mass flows.
        pressures : list or numpy-array length 2
            initial pressures of each fluid.
        enthalpies : list or numpy-array length 2
            initial pressures of each fluid.
        length : float
            heat exchanger length in m.
        diameters : array length 2
            inner diameters of inner tubes and outer tube unit:m.
        U : float, optional
            overall heat transfer coefficient in W /( K m2). The default is 10.
        no_tubes : integer, optional
            number of inner tubes. The default is 1.
        no_points : integer
            No of initial points along x, for solving the bvp.
            
        calc_type : string
            if "const" = does not vary along heat exchanger (default)
            if "calc" = calculated from Nusselt correlation 
            along heat exchanger
        name: string.
            name of the heat-exchanger, default = "HEX_0"

        Returns
        -------
        None.

        """
        
        self.fluids = fluids
        self.mass_flows = mass_flows
        self.pressures = pressures
        self.enthalpies = enthalpies
        self.length = length
        self. diameters = diameters
        self.U = U
        self.no_tubes = no_tubes
        self.no_points = no_points
        self.calc_type = calc_type
        self.name = name
        self.x = np.linspace(0, length, no_points)
        self.area = self.length * self.diameters[0] *self.no_tubes
        qm, qd, f_states = self.q_max()
        self.min_flow = np.where(qd == qm)[0]
        self.h_in = np.linspace(self.enthalpies[0], 
                                    self.enthalpies[0] + qm , no_points)
        self.h_out = np.linspace(self.enthalpies[1] - qm, 
                                    self.enthalpies[1] , no_points)
        
        
    def energy(self,x,h): 
        """
        I must think about it, whether the "self" is needed or helpful
        energy balance 
        
        output: both changes in enthalpy in x-direction
        
        depends on global variables: 
            fl (fluid-names)
            alpha: convection coefficient W/m2/K
            mdot: mass flow rates kg/s
            U: circumference of tube m
        function hps returns an array, the first value is temperature
        """
        T1 = hps(h[1],self.pressures[1],self.fluids[1])[0]
        T0 = hps(h[0],self.pressures[0],self.fluids[0])[0]
        q_konv = T1-T0
        
        dh0 = self.U *self.area / self.mass_flows[0]*q_konv
        dh1 = self.U *self.area / self.mass_flows[1]*q_konv
        return np.array([dh0,dh1])


    def bc(self, ha, hb): #boundary conditions
        # return np.array([ha[0]-ha_in,hb[1]-hb_in,])
        return np.array([ha[0] - self.enthalpies[0], 
                         hb[1] - self.enthalpies[1]])
    
    
    def he_bvp_solve (self):
        y=np.zeros((2,self.no_points))
        y[0,:] = self.h_in
        y[1,:] = self.h_out
        
        result=solve_bvp(self.energy,self.bc, self.x, y, tol=5e-3, 
                         max_nodes=1000)
        return result
    
if __name__  == "__main__":
    mdot=np.array((.0029711, .0351)) # kg/s for both fluids
    alpha = 500  # heat transfer coefficient through the wall (from total resistance)
    fl =["ISOBUTANE","Water"]   # which fluids?
    Tin = [354, 313]
    p = [6.545e5, 4.e5]  # pressure for each fluid, Pa
    ha_in=tp(Tin[0], p[0], fl[0])[2]  # state of fluid 1 left
    hb_in=tp(Tin[1],p[1],fl[1])[2]  # state of fluide 2 right (at L)
    ha_outMax = tp(Tin[1],p[0],fl[0])[2]  # state of fluid 1 left
    hb_outMax = tp(Tin[0],p[1],fl[1])[2]  # state of fluide 2 right (at L)
    heat_ex = heat_exchanger(fl, mdot, p, [ha_in,hb_in])
    qm = heat_ex.q_max()
    print (qm)
    pass
