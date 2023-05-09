# -*- coding: utf-8 -*-
"""
Heat pump cycle incl. heat transfer.

Calculate the pressures and the mass flow rate to fulfill the heat transfer
canditions for given (constant!) temperatures of the heat sink and source
x=1 before the compressor with isentropic efficienc eta and x=0 entering
the throttle. For now, the heat transfer coefficients are constant, but
differing between superheated and saturated in the condenser.
This is done by root finding for a fluid, taking the properties from Refprop 
or CoolProp, this has to be set globally (_props) at the moment, but will perhaps be 
included to the functions later.
Will be used for sensitivity analysis and optimization (hopefully).


Created on Tue Feb 23 17:15:50 2021

@author: atakan
"""

import CoolProp.CoolProp as CP

import fluid_properties_rp as fprop
import numpy as np
# import ht
import scipy.optimize as opti
import matplotlib.pyplot as plt
import os
# p_0 = 1.013e5
# temp_0 = 300
    
fluid_a = "Propane * Pentane * Dimethylether"

x1 =.55
x2 = 0.01
composition =  [x1, x2, 1-x1-x2]  # .950, 0.050][0.9, 0.1]  #
temp = [273, 330] # source and sink Temperature
eta = .65
p = [6e5, 15e5]


FLUIDMODEL = "REFPROP"  # "REFPROP"  or"CoolProp"
_props = FLUIDMODEL
WF = ""
if FLUIDMODEL == "REFPROP":
    os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
    from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
    working_fluid = fluid_a
    WF = fprop.setRPFluid(working_fluid)
    working_fluid = ""
    if _props =="REFPROP":
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        _units = RP.GETENUMdll(0, "MASS BASE SI").iEnum # be careful pressure is in Pa!
    else: _units = 21

elif FLUIDMODEL == "CoolProp":
    # h_0 = CP.PropsSI('H', 'P', p_0, 'T', temp_0, fluid_a)
    working_fluid = CP.AbstractState("HEOS", fluid_a)
    _units = 21
    



def check_bound(p, bounds):
    """
    Check whether the parametere p are within the bounds 

    Parameters
    ----------
    p : numpy.array
        parameters to be checked.
    bounds : list of lists (len(2))
        list with lists of upper and lower limit for each bounded parameter.

    Returns
    -------
    bool
        True when within limits.

    """
    for i in range(len(p)):
        if not((p[i] >= bounds[i][0]) and (p[i] <= bounds[i][1])):
            return False
    return True



def heat_pump_opti(para, eta, U, T_s, working_fluid, bounds, WF):
    """
    Function to optimize (the COP/COP_reversible) of a heat pump with
    isothermal heat source and heat sink. For further details 
    see heat_pump_ht

    Parameters
    ----------
    para : TYPE
        DESCRIPTION.
    eta : TYPE
        DESCRIPTION.
    U : TYPE
        DESCRIPTION.
    T_s : TYPE
        DESCRIPTION.
    working_fluid : TYPE
        DESCRIPTION.
    bounds : TYPE
        DESCRIPTION.
    WF : Instance of Refprop
        fluid mxture name/file are stored inside, output from setRPFluid

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    p = para[:3]
    areas = para[3:5]
    # only binary mixture so far!
    composition = [para[-1], 1-para[-1]]
    #composition = "%1.3f, %1.3f" % (comp[0], comp[1])
    
    
    #print("p,a",p, areas)
    loes = opti.root(heat_pump_ht, p,
                     args=(eta, U, areas, T_s, working_fluid,composition, WF,
                           True, 
                           False),
                     options={"factor": .25})
    print(loes)
    if loes.success:
        p = loes.x[:3]
        A = areas
        
        print("\danach\n",loes.x)
        rat_eff = heat_pump_ht(p, eta, U, A, T_s, working_fluid, composition, WF,
                               rootfind =False, optim=True)
        return -rat_eff
    else:
        return 1e6
        
    



def heat_pump_ht(p, eta, U, A, T_s, working_fluid, composition=[1.0], WF=WF,
                 rootfind =True, optim=False,
                 bounds=None):
    """
    Heat pump cycle with isothermal source and sink, x=0 after condenser.

    x=1 entering the compressor. For optimizing the two working fluid pressures
    or calculating states, COPs etc for optimization, the COP is relative
    to the reversible COP (Carnot).

    Parameters
    ----------
    p : numpy.array, 3 floats
         with the pressures of evaporator and condenser and  the
        mass flow rate (all SI Pa, kg/s.
    eta : float
        compressor efficienvy (isentropic).
    U : numpy.array, 3floats
        mean overall heat ransfer coefficients for evaporator, condenser(gas)
        condenser two-phase (W/m2/K).
    A : numpy.array, 2 floats
        areas of evaporatorand condenser (m2).
    T_s : numpy.array, 2 floats
        secondary fluid temperatures low and high in K.
    working_fluid : coolprop or Refprop fluid
        DESCRIPTION.
    composition : list or array
        mole fractions of each compound
    WF : instance of REFPROPFluid (or empty)
        defines the fluid using REFPROP, output from setRPFluid
    rootfind: Boolean, optional.
        are we finding the root? (T-differences according to heat exchanger size)
    optim : Boolean, optional
        if true: output for optimization
        if false: output of several properties. The default is True.
    bounds: not yet implemented

    Returns
    -------
    if rootfind == True:
        array with three differences of energy balances
    if optim == True:
        rational efficiency (for optimization)
    if False:
          state: array with all calculated states (T,p,x,h,s,rho ...)
          q_w: all specific energies transfered (J/kg)
          A_hg: heat exchanger Area until condensation
          cop: COP of the cycle
          cop_carnot: Carnot-COP for the tempereatures of condenser and
          evaporator (superheating is neglected)
          eta_rational: ratio of the latter two COPs
          p_compr: Compressor Power (W)
          q_h_dot: Heat flow rate at high T (W).

    """
    problem = 0
    druck = False
    n_fluids =len(composition)
    prop_input=[]
    for ii in range(n_fluids):
        prop_input.append([working_fluid[ii], composition[ii], 1, _units, _props, WF[ii]])
        
        # caution when two fluid are used the ii must be checked!
    if n_fluids ==1:
        prop_input = prop_input[0]
    else:
        print("2 fluids are not implemented yet!")
    # ^ all often needed input for the propertties of the working fluid
    if bounds!= None:
        inside = check_bound(p, bounds)
    else:
        inside = True
    if inside:
        p = np.abs(p)
        m_dot = p[2]
        state = np.zeros((8, 6))
        diff = np.zeros((3))
        q_w = np. zeros((4))
        try:
            state[0, :] = fprop.prop_pq(p[0], 1, *prop_input)  # vor dem Kompressor
            if state[0, 0] >= T_s[0]: problem =1
        except:
            state[0, :] = fprop.prop_pq(1.3e5, 1, *prop_input)  # vor dem Kompressor
            print("Fehler!:", p,eta, U, A, T_s, optim)
            
        # Kompressor
        state[1, :] = fprop.sp(state[0, 4], p[1], *prop_input)  # isentropic
        q_w[0] = (state[1, 2] - state[0, 2]) / eta
        
        # kondensatoreintritt:
        state[2, :] = fprop.hp(state[0, 2] + q_w[0], p[1], *prop_input)
        if state[2, 0] < T_s[1]: problem = 2
        # print ("Problem:", problem, p)
        state[6, :] = fprop.prop_pq(p[1], 1, *prop_input)  # Taupunkt
        if state[6, 0] <= T_s[1]: problem = 3
        q_w[1] = (state[6, 2] - state[2, 2])  # Waermeuebtragung bis dahin

        state[7, :] = fprop.prop_pq(p[1], 0, *prop_input)  # gerade kondensiert, Wunsch
        q_w[2] = state[7, 2] - state[6, 2]  # isotherme Kondensation
        state[3, :] = fprop.hp(state[7, 2], p[0], *prop_input)  # hinter Drossel
        q_w[3] = state[0, 2] - state[3, 2]  # Verdampfer
        # - heat transfer:
        q_v = U[0] * A[0] * (T_s[0] - state[0, 0]) / m_dot  # isotherm Verdampfer
        A_hg = q_w[1] / U[1] * m_dot  # Gas-abkühlung, Fläche
        A_hg = -A_hg / ((state[2, 0] - state[6, 0]) /
                      np.log(np.abs((T_s[1] - state[2, 0]) /
                             (T_s[1] - state[6, 0]))))
        q_kond = U[2] * (A[1] - A_hg) * (T_s[1] - state[7, 0]) / m_dot  # isotherme Kondensation

        diff[0] = q_w[3] - q_v  # passt der Druck im Verdampfer zum Waermestrom?
        diff[1] = (state[7, 2] - state[2, 2]) - (q_kond + q_w[1])  # Kondensation incl. Gas
        diff[2] = state[7, 2] - state[6, 2] - q_kond  # nur Kondensation
        # print (m_dot,composition, q_w, q_kond, q_v)
        if problem > 0:
            diff[:] = 1e7
        cop = -q_w[1:3].sum()/q_w[0]
        cop_carnot = T_s[1] / (T_s [1] - T_s[0])
        eta_rational = cop / cop_carnot
        if rootfind:
            return diff
        elif optim:
            if druck: print(eta_rational, m_dot, p,A ,"opti" )
            return eta_rational
        else:
            
            p_compr = q_w[0] * m_dot
            q_h_dot = q_w[1:3].sum() * m_dot
            state =np.delete(state, [1,4,5],0) #  so far unused pointes deleted

            return state, q_w, A_hg, cop, cop_carnot, eta_rational,\
                p_compr,  q_h_dot
    else:
        print("Out of bounds", inside, p, bounds[:len(p)])
        if rootfind:
            return np.ones(3) * 1e5
        else:
            return np.zeros(14)
        
        
class heat_pump:
    
    def __init__(self, fluids, mass_flows, pressures, eta, areas, temp,
               d_in, U, props,  composition, bounds, rootfind,
               optimization=False, 
               calc_type="const", name="heatpump_0", units=21,
               time_dependent=False):
        self.fluids = fluids
        self.fluid_names = fluids.copy()
        self.mass_flows = mass_flows
        self.pressures = pressures
        self.compressor_eta = eta

        self.calc_type = calc_type
        if calc_type == "const":
            self.heat_transfer_coefficient = U
        self.name = name
        self.temperatures = temp
        self.areas = areas
        self.d_in = d_in
        self.composition = composition
        self.bounds = bounds
        self.props = props
        self.name = name
        self.units = units
        # self._abstract_state = abstract_state
        self.rootfind = rootfind
        self.opitimization = optimization
        self. time_dependence= time_dependence
        
        
    def write_yaml(self, fname="st_hex_parameters_file.yaml"):
        """
        Write a yaml file

        Parameters
        ----------
        fname : string, optional
            file name. The default is "st_hex_parameters_file.yaml".

        Returns
        -------
        None.

        """
        
        if self._abstract_state != "":
            _store = self._abstract_state
            self._abstract_state =" "

        with open(fname, mode="wt", encoding="utf-8") as file:
            yaml.dump(self, file)
        self._abstract_state = _store
            
            
    @classmethod
    def read_yaml(cls, fname="st_hex_parameters_file.yaml"):
        """
        read a yaml file

        Parameters
        ----------
        fname : string, optional
            file name. The default is "out_fileN.yaml".

        Returns
        -------
        neu : instance of st_heat_exchanger_input
            to be used as input.

        """
        with open(fname, mode="rt", encoding="utf-8") as file:
            neu = yaml.unsafe_load(file)
            return neu

    def all_out(self):
        """
        format the list output in such a way that it can be used as input 
        for the counterflow heat exchanger.
        If the latter should be modified, this has to be modified too!

        

        """
        n_fluids =len(self.fluids)
        self._abstract_state = []
        if self.props == "REFPROP":
            os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
            from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
            RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
            self.units = RP.GETENUMdll(0, "MASS BASE SI").iEnum # be careful pressure is in Pa!

            for ii in range(n_fluids):
                self._abstract_state.append(fprop.setRPFluid(self.fluids[ii]))
                self.fluids[ii] = ""
               

        elif self.props == "CoolProp":
            
            for ii in range(n_fluids):
                self._abstract_state.append(CP.AbstractState("HEOS", self.fluids[ii]))
            self.units = 21
            

        return (self.pressures, self.mass_flows, 
                self.compressor_eta, 
                  self.heat_transfer_coefficient, 
                self.areas, self.temperatures,
                self.fluids, self.composition, self._abstract_state ,
                self.rootfind, self.opitimization ,
                self.bounds)
    

class heat_pump_eval(heat_pump):
    """ A Heat pump class used for the calculations, including the abstract fluid states
    """ 
    def __init__(self, fluids, mass_flows, pressures, eta, areas, temp,
               d_in, U, props,  composition, bounds, rootfind,
               optimization=False, 
               calc_type="const", name="heatpump_0", units=21,
               time_dependent=False):
        
        super().__init__(fluids, mass_flows, pressures, eta, areas, temp,
                   d_in, U, props,  composition, bounds, rootfind,
                   optimization=False, 
                   calc_type="const", name="heatpump_0", units=21,
                   time_dependent=False)
        n_fluids =len(self.fluids)
        self._abstract_state = []
        if self.props == "REFPROP":
            os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
            from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
            RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
            self.units = RP.GETENUMdll(0, "MASS BASE SI").iEnum # be careful pressure is in Pa!

            for ii in range(n_fluids):
                self._abstract_state.append(fprop.setRPFluid(self.fluids[ii]))
                self.fluids[ii] = ""
               

        elif self.props == "CoolProp":
            
            for ii in range(n_fluids):
                self._abstract_state.append(CP.AbstractState("HEOS", self.fluids[ii]))
            self.units = 21
        
    
    

if __name__ == "__main__":
    diameter = 15e-3  # tube diameter in m
    schleif = False  # calculation for different heat exchanger areas?
    optimierung = False
    U = [1300, 250, 1300]
    areas =[float(6 * np.pi * diameter), float(14 * np.pi * diameter)]  # sehr nichtlinear!
    p0 = [1e5, 22e5, 0.0068]  # Propan np.array([4e5, 19e5, 0.008])
    bou = [(5e4, 4.9e5), (5e5, 3.5e6), (0.0041, 0.09)] #, 
        # (0.8, 9), (1, 17), (0.02, 0.99)]
    names =["pressure_evaporator", "pressure_condenser", "mass_flow_rate" ]
    import yaml
    fname ="heat_pump0.yaml"
    #neu = heat_pump_input(p0[:2], p0[2], eta, U, areas, temp, working_fluid,composition, WF,
    #      True, False, bou )
    neu = heat_pump_input([fluid_a], [p0[2]], p0[:2], eta, areas, temp, [diameter], U, 
                          _props, [composition], bou, True, False)
    
    with open(fname, mode="wt", encoding="utf-8") as file:
        yaml.dump(neu, file)
        # yaml.dump(bou, file)
    loes_input=neu.all_out()
    loes0 = opti.root(heat_pump_ht, [*loes_input[0],*loes_input[1]],
                     args=( *loes_input[2:], ),
                     options={"factor": .25})
 
    loes = opti.root(heat_pump_ht, p0,
                     args=(eta, U, areas, temp, working_fluid,composition, WF,
                           True, False, bou),
                     options={"factor": .25})
    if loes.success:
    
        st, qw, A_hg, cop, cop_carnot, eta_rat, p_compr, q_h_dot = \
            heat_pump_ht(loes.x, eta, U, areas, temp, working_fluid, 
                         composition, WF,
                         False)
        sortierung =[0, 1, 3, 4, 2, 0]  # order of the states along the cycle
        plt.axhline(temp[0])
        plt.axhline(temp[1])
        plt.plot(st[sortierung,4],st[sortierung,0], ":o")
        Vdot = loes.x[-1]*st[0,3]*3600  # m3/h
        P_comp = loes.x[-1]*(st[1,2]-st[0,2])
        print ("Vdot %2.2f m3/h, P = %3.1f W"%(Vdot,P_comp))
        print ("p0: %2.2f bar, p1: =%3.1f bar, mdot: %1.4f (%2.1f)" % 
               (st[0,1]/1e5,st[1,1]/1e5, loes.x[-1],loes.x[-1]*3600))
    
        print("COP: %2.2f, COP(Carnot): %2.2f, eta(rational): %2.2f"
              % (cop, cop_carnot, eta_rat))
        print("P-Compressor/W: %3.2f, Q_dot-highT/W: %3.2f"
              % (p_compr, q_h_dot))
        
    if optimierung:
       
        opti_loes = opti.shgo(lambda x:  heat_pump_opti(x, eta, U, temp, 
                                             working_fluid, bou, WF), bou)
                                    
        print( "\nUsed model:\n", FLUIDMODEL, "\n\n", opti_loes)
        composition_opt = [opti_loes.x[-1],  1-opti_loes.x[-1]]
        st, qw, A_hg, cop, cop_carnot, eta_rat, p_compr, q_h_dot = \
            heat_pump_ht(opti_loes.x[:3], eta, U, opti_loes.x[3:5], 
                         temp, working_fluid, 
                         composition_opt, WF,
                         False)
        sortierung =[0, 1, 3, 4, 2, 0]  # order of the states along the cycle
        plt.axhline(temp[0])
        plt.axhline(temp[1])
        plt.plot(st[sortierung,4],st[sortierung,0], "-o")
    
        print("COP: %2.2f, COP(Carnot): %2.2f, eta(rational): %2.2f"
              % (cop, cop_carnot, eta_rat))
        print("P-Compressor/W: %3.2f, Q_dot-highT/W: %3.2f"
              % (p_compr, q_h_dot))
        
        
            
    if schleif:
        aa = []
        for a_factor in np.arange(1, 100, .75):
            areas = np.array([(1.9 * np.pi * 10e-3), (1 * np.pi * 10e-3)]) * a_factor
            try:
                loes = opti.root(heat_pump_ht, p0,
                                  args=(eta, U, areas, temp, working_fluid, 
                                  composition, WF),
                                  options={"factor": .25})
    
                st, qw, A_hg, cop, cop_carnot, eta_rat, p_compr, q_h_dot = \
                    heat_pump_ht(loes.x, eta, U, areas, temp, working_fluid, 
                    composition_opt, WF, False)
    
                print("COP: %2.2f, COP(Carnot): %2.2f, eta(rational): %1.3f"
                      % (cop, cop_carnot, eta_rat))
                print("P-Compressor/W: %3.2f, Q_dot-highT/W: %3.2f"
                      % (p_compr, q_h_dot))
                p0 = loes.x
            except:
                pass
            aa.append([a_factor, eta_rat, p_compr, cop, q_h_dot, *loes.x])
        aa = np.array(aa)
    
        f, ax = plt.subplots(2, 2, sharex=True)
        ax[0, 0].plot(aa[:, 0], aa[:, 1], "v")
        ax[0, 0].set_title("$\eta_{rat}$")
        ax[1, 1].plot(aa[:, 0], aa[:, 2]/1000, ".-")
        ax[1,  1].plot(aa[:, 0], -aa[:, 4] / 1000, "-")
        ax[1, 1].set_title("P, $\dot Q_h$ / kW")
        ax[1, 0].plot(aa[:, 0], aa[:, -3]/1e5, "v-")
        ax[1, 0].plot(aa[:, 0], aa[:, -2]/1e6, "-")
        ax[1, 0].set_title("$p_l, p_h/10$ / bar")
        ax[0, 1].plot(aa[:, 0], aa[:, -1]*1000, "o")
        ax[0, 1].set_title("$\dot m$ / g/s")
        f.suptitle(fluid_a+", T(l):%3i K, T(h):%3i K" % (temp[0], temp[1]))
        for axx in ax.flat:
            axx.set(xlabel='tubes')  # , ylabel='y-label')
    else:
        print(loes)
