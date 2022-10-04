# -*- coding: utf-8 -*-
"""
Heat exchanger base class 
and counterflow heat exchnager class are introduced
For the first, only the maximum possible heat flow rate s evaluated.
For the counterflow heat exchanger at the moment the boundary value problem 
is solved with constant overall heat transfer coefficient  along the axial 
coordinate is implemented, together with some graphical output.

Fluid properties stem from CoolProp, the low-level interface is used.

Planned: 
    -the convection coefficienets along the axial coordinate 
    shall be evaluated using appropriate 
    Nu-correlations
    - optimization will be implemented (Entropy minimization)
    - perhaps pressure drop will be implemented to the bvp



if nothing else given: counterflow/double-pipe
if nothing else is said given: steady state

Created on Sun Oct  2 16:27:20 2022

@author: atakan
"""
import numpy as np
import CoolProp.CoolProp as CP
from fluid_properties_ll import tp, hps, hps_v
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

class heat_exchanger:
    # Heat exchanger base class
    
    def __init__(self, fluids, mass_flows, pressures, enthalpies, UA=10, 
                 calc_type="const", name="HEX_0"):
        """
        Initializing the heat exchanger with fluids, mass flow rates etc.
        at the moment pressures are constant.

        Parameters
        ----------
        fluids : TYPE list 
            two coolProp-Fluid abstrct states (for low level).
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
        """ maximum possible heat transfer for an isobaric, adiabatic 
        heat exchanger
        """
        
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
        Counter flow heat exchanger class initialization, for double-pipe 
        heat exchangers
        or having no_tubes tubes inside. Now geometrical parameters and 
        a convection coefficient are needed.
        no_points is the initial division of the length for solving the bvp

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
        self.perimeter =self.diameters[0] * np.pi * self.no_tubes
        qm, qd, f_states = self.q_max(1)
        self.min_flow = np.where(qd == qm)[0]
        self.h_in = np.linspace(self.enthalpies[0],  # maximum changes in enthalpy
                                    self.enthalpies[0] + qm , no_points)
        self.h_out = np.linspace(self.enthalpies[1] - qm, 
                                    self.enthalpies[1] , no_points)
        
        
    def energy(self,x,h): 
        """
        energy balance for solving the boundary value problem
        couples the energy changes of each fluid with the convective heat 
        transfer between both fluids.
        At the moment the convection coefficient from the heat exchanger 
        is used;
        this shall be evaluated later as a function of local
        heat exchanger parameters.
        
        output: both changes in enthalpy in x-direction
        
        depends on heat-exchanger variables: 
            fluids (fluid-names)
            U: heat transfer coefficient W/m2/K
            mass_flows: mass flow rates kg/s
            perimeter: circumference of tube m
        function hps returns an array, the first value is temperature
        """
        T1 = hps_v(h[1], self.pressures[1], self.fluids[1])[0]
        T0 = hps_v(h[0], self.pressures[0], self.fluids[0])[0]
        q_konv = T1-T0
        
        dh0 = self.U *self.perimeter / self.mass_flows[0]*q_konv
        dh1 = self.U *self.perimeter / self.mass_flows[1]*q_konv
        return np.array([dh0, dh1])


    def bc(self, ha, hb): 
        """
        two boundary conditions for bvp-solver (scipy-optimize) needed
        here: the enthalpies of the inner fluid at x=0 and the enthalpy of 
        the outer fluid at the end of the heat exchanger at x=length

        Parameters
        ----------
        ha : numpy array
            enthalpies of the inner fluid along the x-coordinate.
        hb : numpy array
            enthalpies of the outer fluid along the x-coordinate.

        Returns
        -------
        numpy array of length 2
           the difference to the prescribed entrance conditions, should both 
           get zero, if succesful.

        """
       
        return np.array([ha[0] - self.enthalpies[0], 
                         hb[1] - self.enthalpies[1]])
    
    
    def he_bvp_solve (self):
        """
        solving the boundary volume problem (scipy) for the counter-flow
        heat exchanger

        Returns
        -------
        result : dictionary
           see scipy documentation, e.g. success (did it find a solution,
                                                  the x and y values)

        """
        y=np.zeros((2,self.no_points))
        y[0,:] = self.h_in
        y[1,:] = self.h_out
        
        result=solve_bvp(self.energy,self.bc, self.x, y, tol=5e-3, 
                         max_nodes=1000)
        return result
    
    
    def he_state(self, result, option = 0):
        """
        Evaluation of all temperatures, enthalpies, entropies etc. 
        along the counterflow heat exchanger and plotting T if option >1

        Parameters
        ----------
        result : TYPE
            DESCRIPTION.
        option : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        if success == True
        states_0 numpy array
           resolved state variables for the inner tube(s).
        states_1 numpy array
           resolved state variables for the outer tube.
        ds float
            entropy production per mass of the heat exchanger (J/( kg K).

        """
        if result.success:
            states_0 = hps_v(result.y[0],self.pressures[0], self.fluids[0])
            states_1 = hps_v(result.y[1], self.pressures[1], self.fluids[1])
            s0 = states_0[4]
            s1 = states_1[4]
            if option > 1:
                plt.figure()
                plt.plot(result.x, states_0[0])
                plt.plot(result.x, states_1[0])
                plt.show()
            ds = (s0[0] - s0[-1]) * self.mass_flows[0] + \
                 (s1[0] - s1[-1]) * self.mass_flows[1]
            
            if option >6:
                print("Entropieproduktion:%3.2f, relativ: %2.2f" % (ds))
            return states_0, states_1, ds
        else: 
            print("Fehler, keine Lösung!", result.message)
            return -1, -1, -1
            
    
if __name__  == "__main__":
    
    T0 = 283.  # K
    mdot=np.array((.0029711, .01)) # kg/s for both fluids
    alpha = 500  # heat transfer coefficient through the wall (from total resistance)
    # Isobutane (hot) and water (cold)
    fl1 = CP.AbstractState("BICUBIC&HEOS", "ISOBUTANE")
    fl2 = CP.AbstractState("BICUBIC&HEOS", "Water")
    fl =[fl1,fl2]   # which fluids?
    Tin = [384, 313]  # initial fluid temperatures, assuming single phae each!
    p = [7.9e5, 4.e5]  # pressure for each fluid, Pa
    
    #  evaluate enthalpies and maximum possible enthalpy changes:
    ha_in = tp(Tin[0], p[0], fl[0])[2]  # state of fluid 1 left
    hb_in=tp(Tin[1],p[1],fl[1])[2]  # state of fluide 2 right (at L)
    ha_outMax = tp(Tin[1],p[0],fl[0])[2]  # state of fluid 1 left
    hb_outMax = tp(Tin[0],p[1],fl[1])[2]  # state of fluide 2 right (at L)

    diameters =[1.5e-2, 3e-2]  # m
    length = 5.  # m
    
    heat_ex = counterflow_hex(fl, mdot, p, [ha_in,hb_in], 
                              length, diameters, U=alpha, no_tubes=2)  # assign parameters
    res =heat_ex.he_bvp_solve()  # solve the heat exchanger problem
    
    f1,f2,ds =heat_ex.he_state(res,5) # evaluate results (and plot)
    print("Entropy production rate: %2.2e W/K, exergy loss rate %2.2e W" 
          % (ds, ds * T0))
