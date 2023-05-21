# -*- coding: utf-8 -*-
"""
Class and function for static heat exchanger

Means evaporator and condenser (incl. superheating or subcooling) of a working
fluid and a secondary fluid. Both with properties from a model like
Refprop, and both can be mixtures.
Everything is calcultaed per kg of the working fluid (specific). For the 
secondary fluid the mass-flow ratio is calculated.
No time dependence, only thermodynamics, no heat transfer calculation beyond 
energy conservation (means: no heat exchanger layut or sizing).

Part of carbatpy.

Created on Fri May 19 14:57:09 2023

@author: atakan
"""


from scipy.optimize import minimize, differential_evolution
import numpy as np
# import CoolProp.CoolProp as CP
import fluid_properties_rp as flp
# from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
props = "REFPROP"


# import components.heat_exchanger as hex


class static_heat_exchanger:
    """
    Class for static heat exchanger
    
    means: no time dependence and not heat transfer coefficients * areas
    are used!
    Only first law and second law will be checked (the latter must be improved)
    
    """
    
    
    def __init__(self, fluids, temps, ps, qs=[0, -2], h_enter =[-1e9, -1e9], 
                 points=30, dT_hex=0.5,
                 dT_superh=10, heating=True, dH_min=0, props="REFPROP",
                 compositions=[[1.0], [1.0]], calc_type="const",
                 name="evaporator",
                 units=21):  
        """
        Initialization of the class static heat exchanger
        
        Will be used for two fluids flowing in counterflow, at the moment
        assuming that the high temperature fluid in the evaporator and the low-
        temperature fluid in the condenser do not go through phase changes. 
        The naming (condenser/evaporator) below, assumes a heat pump. If it is 
        used for a power cycle, one has to invert the naming. The enthalpies
        and their difference is always specific for the working fluid
        (first fluid) for the second fluid, the mass ratio is calculated.

        Parameters
        ----------
        fluids : list of two strings
            Fluid names according to Refprop each.
        temps : list of two floats
           two temperatures (K), the first is always used as the highest (saturation)
           T of the working fluid with quality of one (dT_superh is added to it), 
           the second can be used as the lowest T of the secondary
           (or storage) fluid.
        ps : list of two floats
            The prescribed pressures(Pa) of the two fluids, for the evaporator
            the first pressure can be varied by the program.
        qs : list of two floats, optional
            The quality of the two fluids, for pure liquid negative,
            for pure vapor larger 1. The default is [0, -2].
        h_enter : list of two floats, optional
            if the values are not -1e-9, the value is used as the entering 
            enthalpy for the evaporator. The default is [-1e-9, -1e-9].
        points : integer, optional
            number of state points to evaluate the states of both fluids.
            The default is 30.
        dT_hex : float, optional
            the minimum temperature difference between the two fluids, which
            is allowed for any point. The default is 0.5.
        dT_superh : float, optional
            superheating of the working fluid. The default is 10.
        heating : boolean, optional
            if True, we heat the working fluid (evaporator), if False, we cool 
            it(condenser. The default is True.
        dH_min : float, optional
            minimum specific enthalpy change of the working fluid. May be 
            prescribed by other parts of the cycle (e.g.:evaporator + compressor) 
            if it is set to zero, the enthalpy difference results from the 
            given T/p-states of the secondary fluid. The default is 0.
        props : string, optional
            either REFPROP or CoolProp (what shall be used?) Refprop is safe at
            the moment. The fluid naming conventions also depond on the choice.
            The default is "REFPROP".
        compositions : list of list of floats, optional
            the composition (mole fractions) of each of the two fluids, if they 
            are mixtures. Otherwise both: 1.0. The default is [[1.0], [1.0]].
        calc_type : string, optional
            at the moment unused, may be later used for different variations. 
            The default is "const".
        name : string, optional
            name of the heat exchanger, helps distinguishing when coupled to a 
            cycle. Could be that "evaporator" and "condenser" may get a meaning 
            ater.
            The default is "evaporator".
        units : integer, optional
            indication of the unit system of evrything in REFPROP. 21 means
            that SI (kg, J, s, W, Pa, etc.) is used throughout, The default is 21.

        Returns
        -------
        None.

        """
        
        # qs def.: <0=> liquid entering, qs >1 => vapor entering

        self.temps = temps
        self.ps = ps
        self.qs = qs
        self.h_enter = h_enter
        self.points = points
        self.dT_hex = dT_hex
        self.dT_superh = dT_superh
        self.heating = heating
        self.dH_min = dH_min
        self.compositions = compositions
        self.props = props
        self.units = units
        self.calc_type = calc_type
        self.name = name
        self.fluids = fluids

    def pinchpoint(self, verbose=False):
        """
        Calculate the state of a static heat exchanger
        
        only from thermodynamics. If it is an evaporator, the pressure may be 
        varied, such that the exiting working fluid is really at the wanted 
        superheated state. The quality of the entering working fluid can be 
        prescribed. For the heat-pump evaporator it may be better to set the 
        enthalpy of the entering working fluid. Finally the values of the 
        static_heat_exchanger instance are set.

        Parameters
        ----------
        verbose : TYPE, optional
            DESCRIPTION. The default is False.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.
        But following values are set:
        self.enthalpies : list of two arrays
            all specific enthlpies of the two fluids for the set number of points
        self.m_ratio: float
            The ratio of the mass flow between the two fluids
        self.dh : list of two floats
            the specific enthalpy chnages of the two fluids.
        self.t_all: list of two arrays
            with all specific values (T, p, h, v, s, q) of both fluids along 
            the set number of points

        """

        from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary as modwf

        from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary as modsf
        self.RP = [flp.setRPFluid(self.fluids[0], modwf, 'RPPREFIX'),
                   flp.setRPFluid(self.fluids[1], modsf, 'RPPREFIXs')]
        sat_p_initial = self.ps[0]
        flp_rest =[[self.fluids[0], self.compositions[0],1, 
                    self.units,self.props,self.RP[0]],
                   [self.fluids[1], self.compositions[1],1, 
                    self.units,self.props,self.RP[1]]]
        sat_p = flp.prop_Tq(self.temps[0], 1., *flp_rest[0])
                            
        self.dT_hex_wanted = self.dT_hex
        
        if self.heating: # the working fluid (evaporator in the heat pump)
            # self.ps[0] = sat_p[1]
            # print(
            #     f"Saturation pressure varied from {sat_p_initial:.2e} to {self.ps[0] :.2e} Pa!")
            epsilon = -100
            count =0
            while epsilon< 0.1 and count <7:
                count+=1
                sat_v = flp.prop_pq(self.ps[0], 1.,*flp_rest[0])
                                    
                if self.qs[0] < 0 or self.qs[0] > 1:
                    raise Exception(f"working fluid quality is wrong!{self.qs}!")
                    
                if self.h_enter[0] < -1e8:  # quality will be used                         
                    sat_l = flp.prop_pq(self.ps[0], self.qs[0],*flp_rest[0])
                    sf_in = flp.tp(sat_v[0]+self.dT_superh + self.dT_hex, 
                                   self.ps[1],*flp_rest[1])
                    sf_out = flp.tp(sat_l[0] + self.dT_hex, self.ps[1],
                                    *flp_rest[1])
                                        
                
                                          
                    ev_out = flp.tp(sat_v[0]+self.dT_superh, self.ps[0], 
                                *flp_rest[0])
                else: # enthalpy of the entering working fluid given.
                    sat_l = flp.hp(self.h_enter[0], self.ps[0],*flp_rest[0])
                    sf_in = flp.tp(sat_v[0]+self.dT_superh + self.dT_hex, 
                                   self.ps[1],*flp_rest[1])
                    sf_out = flp.tp(sat_l[0] + self.dT_hex, self.ps[1],
                                    *flp_rest[1])
                    ev_out = flp.tp(sat_v[0]+self.dT_superh, self.ps[0], 
                                *flp_rest[0])
                    
                if self.dH_min>0:
                    ev_out = flp.hp(sat_l[2] + self.dH_min, self.ps[0], 
                                    *flp_rest[0])
                    sf_in = flp.tp(ev_out[0] + self.dT_hex, self.ps[1],
                                    *flp_rest[1])
                    sf_out = flp.tp(sat_l[0] + self.dT_hex, self.ps[1],
                                    *flp_rest[1])
                    
    
                
    
                if verbose:
                    print(sat_v, sat_l, ev_out, "\nsf:", sf_in, sf_out)
    
                dh0 = ev_out[2] - sat_l[2]
                dh1 = sf_in[2] - sf_out[2]
                m_ratio = dh1 / dh0
                h_0 = np.linspace(sat_l[2], ev_out[2], self.points)
                h_1 = np.linspace(sf_out[2], sf_in[2], self.points)
                t0 = flp.hp_v(h_0, self.ps[0],*flp_rest[0])
                t1 = flp.hp_v(h_1, self.ps[1],*flp_rest[1])
                self.enthalpies = [h_0, h_1]
                self.m_ratio = m_ratio
                self.dh = [dh0, dh1]
                self.t_all = [t0, t1]
                self.Characteristic_points = [sat_l, sat_v, ev_out, sf_in, sf_out]
                h0_shifted = self.t_all[0][2, :] - \
                    self.t_all[0][2, 0], self.t_all[0][0, :]
                h1_shifted = (self.t_all[1][2, :] - self.t_all[1]
                              [2, 0])/self.m_ratio, self.t_all[1][0, :]
                dT_min = np.min(h1_shifted[1]-h0_shifted[1])
                epsilon =  self.dT_hex - dT_min
                if epsilon < -.010:
                    
    
                    print(f"Problem{dT_min, self.dT_hex, epsilon} heating")
                    self.dT_hex = dT_min - epsilon
        else:  # cooling the working_fluid (condenser in heat pump)
            dT_min =100
            count = 0
            while dT_min > self.dT_hex and count <10:
                count+=1
                if self.h_enter[0] < -1e8: # T, p will be used
                    con_in = flp.tp(self.temps[0]+self.dT_superh+self.dT_hex, self.ps[0],
                                    *flp_rest[0])
                else:  # enthalpie will be used
                    con_in = flp.hp(self.h_enter[0], self.ps[0],
                                    *flp_rest[0])
                con_dew = flp.prop_pq(self.ps[0], 1.,*flp_rest[0])
                dew_T = con_dew[0] # here it must be prooved whether the pinch point T-difference is ok, if not: pressure change
                sf_out = flp.tp(self.temps[0]+self.dT_superh-self.dT_hex, self.ps[1],
                                *flp_rest[1])
                sf_dew = flp.tp(dew_T-self.dT_hex, self.ps[1], *flp_rest[1])
                m_ratio =(sf_out[2]-sf_dew[2])/(con_in[2]-con_dew[2])
                
                if self.dH_min ==0:
                       
                    sf_in = flp.tp(self.temps[1]-self.dT_hex, self.ps[1],
                                   *flp_rest[1])
                    dh1 = sf_out[2] -sf_in[2]
                    dh0 = dh1 / m_ratio
                    con_out = flp.hp(con_in[2]-dh0, self.ps[0],*flp_rest[0])
                else:  # prescribed enthalpy difference will be used
                    dh_remain = np.abs(self.dH_min)-(con_in[2]-con_dew[2])
                    con_out = flp.hp(con_in[2]-self.dH_min, self.ps[0],
                                    *flp_rest[0])
                    sf_in = flp.tp(con_out[0]-self.dT_hex, self.ps[1],
                                   *flp_rest[1])
                    m_ratio = (sf_dew[2]-sf_in[2])/(dh_remain)
                    sf_out = flp.hp(sf_in[2]+self.dH_min*m_ratio, self.ps[1],
                                   *flp_rest[1])
                    dh0 =self.dH_min
                    dh1 =dh0 / m_ratio # BA hier pruegfen!!
                    
                                                
                h_0 = np.linspace(con_out[2], con_in[2], self.points)
                h_1 = np.linspace(sf_in[2], sf_out[2], self.points)
                t0 = flp.hp_v(h_0, self.ps[0],*flp_rest[0])
                t1 = flp.hp_v(h_1, self.ps[1],*flp_rest[1])
                self.enthalpies = [h_0, h_1]
                self.m_ratio = m_ratio
                self.dh = [dh0, dh1]
                self.t_all = [t0, t1]
                self.characteristic_points = [con_in, con_out, con_dew,sf_dew, sf_in, sf_out]
                h0_shifted = self.t_all[0][2, :] - \
                    self.t_all[0][2, 0], self.t_all[0][0, :]
                h1_shifted = (self.t_all[1][2, :] - self.t_all[1]
                              [2, 0])/self.m_ratio, self.t_all[1][0, :]
                dT_min = np.min(h0_shifted[1]-h1_shifted[1])
                if dT_min < self.dT_hex*.99:
                    
    
                    print(f"Problem:{dT_min, self.dT_hex} cond")
                    self.dT_hex+=dT_min

    def pp_root(self, var_in):
        if self.heating:
            if self.qs[0] >= 0 and self.qs[0] <= 1:
                pass
    def hex_plot(self, second="", fname ="hex_plot.png"):
        plt.figure()
        shift1 = self.t_all[0][2, -1] - self.t_all[0][2, 0]
        
        plt.plot(self.t_all[0][2, :] - self.t_all[0][2, 0], self.t_all[0][0, :])
        plt.plot((self.t_all[1][2, :] - self.t_all[1][2, 0]) /
                 self.m_ratio, self.t_all[1][0, :], "o")
        if second !="":
            shift1 = second.t_all[0][2, -1] - second.t_all[0][2, 0]
            
            plt.plot(second.t_all[0][2, :] - second.t_all[0][2, 0], second.t_all[0][0, :])
            plt.plot((second.t_all[1][2, :] - second.t_all[1][2, 0]) /
                     second.m_ratio, second.t_all[1][0, :], "o")
        plt.xlabel("enthalpy change / (J/(kg of working fluid))")
        plt.ylabel("temperature / (K")
        plt.savefig(fname)


if __name__ == "__main__":
    fl1 = "Propane*Pentane*Butane"
    fl2 = "Methanol"
    compositions = [[.3, .4, .3], [1.]]
    
    # Evaporator with given quality
    Ts = [290., 250.]
    # fl = fl_names
    flx = [fl1, fl2]
    ps = [1.92e5,  12e5]
    qs = [0.05, -2]
    # Evaporator

    hp0 = static_heat_exchanger(
        flx, Ts, ps, dT_superh=15, qs=qs, compositions=compositions)
    hp0.pinchpoint()

    plt.figure()
    shift1 = hp0.t_all[0][2, -1] - hp0.t_all[0][2, 0]
    
    plt.plot(hp0.t_all[0][2, :] - hp0.t_all[0][2, 0], hp0.t_all[0][0, :])
    plt.plot((hp0.t_all[1][2, :] - hp0.t_all[1][2, 0]) /
             hp0.m_ratio, hp0.t_all[1][0, :], "o")
    
    # Evporator with given enthalpy (nearly the same calculation, start later)
    h_0 = [hp0.enthalpies[0][2], -1e9]
    hp0 = static_heat_exchanger(
        flx, Ts, ps, dT_superh=15, dT_hex=1,qs=qs, h_enter=h_0,  compositions=compositions)
    hp0.pinchpoint()

    # plt.figure()
    shift1 = hp0.t_all[0][2, -1] - hp0.t_all[0][2, 0]
    
    plt.plot(hp0.t_all[0][2, :] - hp0.t_all[0][2, 0], hp0.t_all[0][0, :], "k")
    plt.plot((hp0.t_all[1][2, :] - hp0.t_all[1][2, 0]) /
             hp0.m_ratio, hp0.t_all[1][0, :], "x")
    
    # Condenser with given enthalpy difference
    print("condenser")
    Ts = [350., 290.]
    ps = [8.92e5,  12e5]
    fl2 = "Water"
    flx = [fl1, fl2]
    qs = [2, -2]
    hp0 = static_heat_exchanger(
        flx, Ts, ps, dT_superh=15, qs=qs,heating =False,dH_min=hp0.dh[0],
        compositions=compositions)
    hp0.pinchpoint()

    shift = shift1 -( hp0.t_all[0][2, -1] - hp0.t_all[0][2, 0])
    plt.plot(hp0.t_all[0][2, :] - hp0.t_all[0][2, 0]+shift, hp0.t_all[0][0, :],".")
    plt.plot((hp0.t_all[1][2, :] - hp0.t_all[1][2, 0]) /
              hp0.m_ratio +shift, hp0.t_all[1][0, :], "v")
