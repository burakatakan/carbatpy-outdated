# -*- coding: utf-8 -*-
"""
Functions to obtain fluid properties from REFPROP or CoolProp using the low 
level interface. (Set _props accordingly)
Mainly used for thermodynamic cycle calculations.
So far for several combinations of (input) parameters in the single-phase or 
saturated regime. For the ouput it can be mostly chosen, whether transport 
properties are als returned, besides the thermodynamiic properties. 
Functions are also vectorized by hand (Name:"..._v). Physical exergies 
can be calculated for a state of known h and p.

For REFPROP two usages of hp, Tp etc. are possible:
Passing a string with the fluid name and compositon, properties etc.
(This sometimes leads to trouble, when having many function calls >1E5 and
it is slower). selected instance (RP) will be used throughout.

Or:First calling setPRFluid with the fluid name string. This generates an 
instance of REFPROP for this fluid (mixture), which has to be passed to the 
calls of hp, Tp etc. together with an empty string as fluid name.

Standard units: SI (J, kg, Pa, K,...)

"""


from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import numpy as np
import CoolProp.CoolProp as CP
import os
from time import time

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
_props = "REFPROP"  # or "CoolProp"
_fl_properties_names = ("Temperature", "Pressure", "spec. Enthalpy",
                        "spec. Volume","spec. Entropy", "quality", 
                        "spec. heat capacity", 
                        "viscosity", "thermal conductivity", 
                        "Prandtl number","k.viscosity")
# order for coolprop,alle_0:[_temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, prandtl, phase]"
if _props == "REFPROP":
    try:
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        # be careful pressure is in Pa!
        _units = RP.GETENUMdll(0, "MASS BASE SI").iEnum
    except:
        print("print refprop not installed!")
        
else:
    _units = 21


__Tenv__ = 283.15  # Temp. of the environment in K
__penv__ = 1.013e5  # Pressure of the environment in Pa

def setRPFluid(fluid):
    """
    A new instnce of Refpropdll for the given fluid. It can then be called using fluid =""

    Parameters
    ----------
    fluid : string
        fluid (mixture) name, as described in REFPROP.

    Returns
    -------
    SP : REFPROP Instance
        for further usage.

    """
    SP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) #secondary fluid
    SP.SETPATHdll(os.environ['RPPREFIX'])
    SP.SETFLUIDSdll(fluid)
    return SP
    
    
def mdot_area_function(m_dot, diameter):
    """
    Calculate the mass flux through a tube and the orthogonal area 

    Parameters
    ----------
    m_dot : float
        mass flow rate (kg/s).
    diameter : float
        of the tube (m).

    Returns
    -------
    m_dot_area : float
        in kg/(s m2).
    area : float
        in m2.

    """
    area = np.pi * (diameter / 2)**2
    m_dot_area = m_dot / area
    return m_dot_area, area


def hp_exergy(h, p, fluid="", T_env=__Tenv__, p_env=__penv__, props=_props,
              composition=[1.0], RP=RP):
    """
    calculate the specific exergy of a fluid from given enthalpy and pressure

    Parameters
    ----------
    h : float
        specific enthalpy (J/kg).
    p : float
        pressure (Pa).
    fluid : depends on model, as defined with props(REFPROP, CoolProp etc.)
        define the fluid(mixture).
    T_env : float, optional
        temperature (K) of the environment. The default is __Tenv__.
    p_env : float, optional
        pressure (Pa) of the environment. The default is __penv__.
    props : TYPE, optional
        define the fluid module to use. The default is _props.
    composition : depends on model, optional
        mole fractions of the different compounds, should add to 1. 
        The default is [1.0].
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, instead flud ="" works.

    Returns
    -------
    ex : float
        specific exergy (J/kg).

    """
    pr_test = False
    state_env = tp(T_env, p_env, fluid, option=1, props=props,
                   composition=composition)
    state = hp(h, p, fluid, props=props, composition=composition)
    dstate = state - state_env
    ex = dstate[2] - __Tenv__ * dstate[4]
    if pr_test:
        print(state, state_env, dstate, "\n")
    return ex


def hp(h, p, fluid="", composition=[1.0], option=1, units=_units, props=_props, RP=RP):
    """
    Properties needed for integration at given  values of pressure 
    and specific enthalpy.

    Parameters
    ----------
    h : float
        specific enthalpy in J/kg.
    p : float
        pressure in Pa.
    fluid :   an AbstractState in coolprop. or fluid name in Refprop
    option: integer
        if 0 also transport properties will be calculated 
        (not yet implemented for REFPROP)
        with option 1 : T,p,h,v,s,x is calculated
    units: integer
        units in Refprop (imposrtant: must be SI and mass based (kg, Pa etc.))      
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, instead flud ="" works.


    Returns
    -------
    alle : numpy array
        if option = 0 it includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

        if option =1: it is compatible to the high level output
        T,p,h,v,s,x.

    """
    if props == "REFPROP":

        if option == 0:
            o = RP.REFPROP2dll(
                fluid, "HP", "T;D;S;QMASS;CP;VIS;TCX;PRANDTL;KV", units, 0, h, p, composition)
            alle = [o.Output[0], p, h, 1/o.Output[1], *o.Output[2:9]]
            if alle[8] < 0:
                print("error in hp, probably saturated")
                print("alle in hp", alle, o)
        if option == 1:
            o = RP.REFPROP2dll(fluid, "HP", "T;D;S;QMASS",
                               units, 0, h, p, composition)
            
            alle = [o.Output[0], p, h, 1./o.Output[1], *o.Output[2:4]]

    elif props == "CoolProp":
        fluid.update(CP.HmassP_INPUTS, h, p)
        reihe = [CP.iT, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                 CP.iconductivity, CP.iCpmass, CP.iviscosity]
        props = [fluid.keyed_output(k) for k in reihe]
        _temp, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]

        if option == 0:
            print("warning, check variables order")
            alle = [_temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, prandtl, phase]
        elif option == 1:
            alle = [_temp, p,  h, 1/rho, s, x]
    else: print(f"Warning: props not set or invalid, actual value: {props}")
    return np.array(alle)


def uv(u, v, fluid="", composition=[1.0], option=1, units=_units, props=_props, RP=RP):
    """
    Properties needed for integration at given  values of internal energy 
    and specific volume.

    Parameters
    ----------
    u : float
        internal energy in J/kg.
    v : float
        specific volume in m3(kg.
    fluid :   an AbstractState in coolprop. or fluid name in Refprop
    option: integer
        if 0 also transport properties will be calculated 
        (not yet implemented for REFPROP)
        with option 1 : T,p,h,v,s,x is calculated
    units: integer
        units in Refprop (important: must be SI and mass based (kg, Pa etc.))
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, instead flud ="" works.

    Returns
    -------
    alle : numpy array
        if option = 0 it includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

        if option =1: it is compatible to the high level output
        T,p,h,v,s,x.

    """
    if props == "REFPROP":

        if option == 0:
            o = RP.REFPROP2dll(
                fluid, "ED", "T;p;S;q;CP;VIS;TCX;PRANDTL;KV;H", units, 0, u, 1/v, composition)
            alle = [*o.Output[0:2], o.Output[9], v,o.Output[2], *o.Output[2:9]]
            if alle[8] < 0:
                print("error in uv, probably saturated")
        if option == 1:
            o = RP.REFPROP2dll(fluid, "ED", "T;p;S;QMASS;H",
                               units, 0, u, 1/v, composition)
            alle = [*o.Output[0:2],  o.Output[4], v,*o.Output[2:4]]

    elif props == "CoolProp":
        print (" uv not implemented yet for coolprop")
        # fluid.update(CP.HmassP_INPUTS, h, p)
        # reihe = [CP.iT, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
        #          CP.iconductivity, CP.iCpmass, CP.iviscosity]
        # props = [fluid.keyed_output(k) for k in reihe]
        # _temp, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]

        # if option == 0:
        #     alle = [_temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, prandtl, phase]
        # elif option == 1:
        #     alle = [_temp, p,  h, 1/rho, s, x]
    return np.array(alle)


def sp(s, p, fluid="", composition=[1.0], option=1, units=_units, props=_props, RP=RP):
    """
    Properties needed for integration at given s and p, single phase.

    Parameters
    ----------
    s : float
        specific entropy in J/kg.
    p : float
        pressure in Pa.

    fluid :   an AbstractState in coolprop. or fluid name in Refprop
    option: integer
        if 0 also transport properties will be calculated 
        (not yet implemented for REFPROP)
        with option 1 : T,p,h,v,s,x is calculated
    units: integer
        units in Refprop (imposrtant: must be SI and mass based (kg, Pa etc.))
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, instead flud ="" works.


    Returns
    -------
    alle : numpy array
        if option = 0 it includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

        if option =1: it is compatible to the high level output
        T,p,h,v,s,x.

    """
    if props == "REFPROP":
        o = RP.REFPROP2dll(fluid, "SP", "T;D;H;QMASS", units, 0, s, p, composition)
        if option == 0:
            alle = []
        if option == 1:
            alle = [o.Output[0], p, o.Output[2],
                    1/o.Output[1], s, o.Output[-1]]

    elif props == "CoolProp":
        fluid.update(CP.PSmass_INPUTS, p, s)  # Pruefen
        reihe = [CP.iT, CP.iQ, CP.iHmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                 CP.iconductivity, CP.iCpmass, CP.iviscosity]
        props = [fluid.keyed_output(k) for k in reihe]
        _temp, x, h, rho, prandtl, phase, lambda_s, cp, mu = props[:]

        if option == 0:
            print("warning, check variables")
            alle = [_temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, prandtl, phase]
        elif option == 1:
            alle = [_temp, p,  h, 1/rho, s, x]
    return np.array(alle)


name_properties = [
    ["temperature", "p", "x", "h",  "s", "rho", "mu", "cp", "lambda_s",
     "phase", "prandtl"],
    ["temperature", "p",  "h", "v", "s", "x"]
]


def hp_v(h, p, fluid="", composition=[1.0], option=1, units=_units, props=_props, RP=RP):
    """ 
    Vectorization of the single phase properties function, see: hp
    
    """
    _n = len(h)
    if option == 1:
        alle = np.zeros((6, _n))
    else:
        alle = np.zeros((11, _n))
    for _i in range(_n):
        if np.isscalar(p):
            alle[:, _i] = hp(h[_i], p, fluid, composition, option,
                             units, props)
        else:
            alle[:, _i] = hp(h[_i], p[_i], fluid, composition, option,
                             units, props)
    return alle


def tp(temp, p,  fluid="", composition=[1.0], option=1, units=_units,
       props=_props, RP=RP):
    """
    Properties needed for integration at given temperature and pressure, 
    single phase. 
    option 0 (many properties) not yet implemented for REFPROP!

    Parameters
    ----------

    temp : float
         temperature in K.
    p : float
        pressure in Pa.

    fluid :   an AbstractState in coolprop.
    RP : REFPROP Instance as set with setRPfluid
    if set, the fluid name does not have to be passed again, instead flud ="" works.


    Returns
    -------
    alle : numpy array
        if option = 0 it includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

        if option =1: it is compatible to the high level output
        T,p,h,v,s,x.
        
        if option =2: molecular mass MM is included
        T,p,h,v,s,x, MM.

    """
    if props == "REFPROP":
        o = RP.REFPROP2dll(fluid, "TP", "H;D;S;QMASS", units,
                           0, temp, p, composition)
        if option == 0:
            alle = []

        elif option == 1:
            alle = [temp, p, o.Output[0], 1 / o.Output[1], *o.Output[2:4]]
        elif option == 2: ## for moelcular mass ouput:
            no_compounds = len(composition)
            mm_string = ";MM"*no_compounds
            o = RP.REFPROP2dll(fluid, "TP", "H;D;S;QMASS"+ mm_string, units,
                               0, temp, p, composition)
            mol_mass_mix = (o.Output[4:4+no_compounds]*np.array(composition)).sum()
            
            alle = [temp, p, o.Output[0], 1 / o.Output[1], 
                    *o.Output[2:4+no_compounds], mol_mass_mix]
    elif props == "CoolProp":
        fluid.update(CP.PT_INPUTS, p, temp)
        reihe = [CP.iHmass, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                 CP.iconductivity, CP.iCpmass, CP.iviscosity]
        props = [fluid.keyed_output(k) for k in reihe]
        h, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]

        if option == 0:
            print("warning, check variables")
            alle = np.array([temp, p,  h, 1/ rho, s,x,cp, mu,  \
                             lambda_s, prandtl, phase])
        elif option == 1:
            return [temp, p,  h, 1/rho, s,  x]

    return np.array(alle)
#  below must be checked!


def p_prop_sat(p,  fluid="", composition=[1.0], option=1, units=_units,
               props=_props, RP=RP):
    """
    Saturation state properties at given p for a certain fluid (mixture).

    Parameters
    ----------
    p : float
        saturation pressure in Pa.

    fluid :   an AbstractState in coolprop or a fluid in REFPROP.
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, 
        instead flud ="" works.


    Returns
    -------
    alle :  numpy array (2,4)
            includes:  properies in saturated state at given pressure p
            of  vapor(0,:),and liquid (1,:) 
            [1] T, p, h, v,s, q
            [0]  T p h v s q cp, viscosity, thermal conductivity, Pr, kinematic viscosity
            all in SI units.

    """
    vap_liq = []
    for qq in [1, 0]:
        if props == "REFPROP":

            if option == 0:
                o = RP.REFPROP2dll(
                    fluid, "PQ", "T;H;D;S;CP;VIS;TCX;PRANDTL;KV", units, 0, p, qq, composition)
                alle = [o.Output[0], p, o.Output[1], 1 /
                        o.Output[2], o.Output[3], qq, *o.Output[4:9]]

            elif option == 1:
                o = RP.REFPROP2dll(fluid, "PQ", "T;H;D;S",
                                   units, 0, p, qq, composition)
                alle = [o.Output[0], p, o.Output[1],
                        1 / o.Output[2], o.Output[3], qq]

        elif props == "CoolProp":
            fluid.update(CP.PQ, p, qq)
            reihe = [CP.iT, CP.iHmass, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                     CP.iconductivity, CP.iCpmass, CP.iviscosity]
            props = [fluid.keyed_output(k) for k in reihe]
            temp, h, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]

            if option == 0:
                alle = [temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, 
                        prandtl, phase]
                print("warning, check variables")
            elif option == 1:
                alle = [temp, p,  h, 1/rho, s,  x]

        vap_liq.append(np.array(alle))
        
    return np.array(vap_liq)


def T_prop_sat(temp,  fluid="", composition=[1.0], option=1, units=_units,
               props=_props, RP=RP):
    """
    Saturation state properties at given temperature for a certain fluid 
    (mixture).

    Parameters
    ----------
    temp : float
        Temperature in K.

    fluid :   depends
        an AbstractState in coolprop or a fluid in REFPROP.
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, 
        instead flud ="" works.


    Returns
    -------
    alle :  numpy array (2,4)
        includes:  properies in saturated state at given pressure p
        of liquid (0,:) and vapor(1,:),
        [1] T, p, h, v,s, q
        [0]  T p h v s q cp, viscosity, thermal conductivity, Pr, 
        kinematic viscosity
        all in SI units.

    """
    
    vap_liq = []
    alle=0
    for qq in [1, 0]:
        if props == "REFPROP":

            if option == 0:
                o = RP.REFPROP2dll(
                    fluid, "TQ", "P;H;D;S;CP;VIS;TCX;PRANDTL;KV", units, 0, 
                    temp, qq, composition)
                alle = [temp, o.Output[0], o.Output[1], 1 /
                        o.Output[2], o.Output[3], qq, *o.Output[4:9]]

            elif option == 1:
                o = RP.REFPROP2dll(fluid, "TQ", "P;H;D;S",
                                   units, 0, temp, qq, composition)
                alle = [temp, o.Output[0],  o.Output[1],
                        1 / o.Output[2], o.Output[3], qq]

        elif props == "CoolProp":
            fluid.update(CP.QT,  qq, temp)
            reihe = [CP.iP, CP.iHmass, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                     CP.iconductivity, CP.iCpmass, CP.iviscosity]
            props = [fluid.keyed_output(k) for k in reihe]
            p, h, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]

            if option == 0: # wrong order must be checked
                alle = [temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, 
                        prandtl, phase]
                print("warning, check variables")

            elif option == 1:
                alle = [temp, p,  h, 1/rho, s,  x]

        vap_liq.append(np.array(alle))
        
    return np.array(vap_liq)


def prop_pq(p, q, fluid="", composition=[1.0], option=1, units=_units,
            props=_props, RP=RP):
    """
    Saturation state properties at given pressure  and quality for a certain 
    fluid (mixture).

    Parameters
    ----------
    p : float
        pressure in Pa.
    q: float (0<q<1)
        quality

    fluid :   an AbstractState in coolprop or a fluid in REFPROP.
    composition: list of floats
        mole fraction of each fluid.
    option: integer
        [1] only T p h v s and q

    units: integer
        select units (SI of property model, for REFPROP: 21)
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, 
        instead flud ="" works.


    Returns
    -------
    alle : numpy array (2,4)
        includes:  properies in saturated state at given pressure p
        of liquid (0,:) and vapor(1,:),
        T, p, h, v,s, q
        all in SI units.

    """

    if props == "REFPROP":
        o = RP.REFPROP2dll(fluid, "PQ", "T;H;D;S", units, 0, p, q, composition)
        if option == 0:
            alle = []

        elif option == 1:
            alle = [o.Output[0], p, o.Output[1],
                    1 / o.Output[2], o.Output[3], q]

    elif props == "CoolProp":
        fluid.update(CP.PQ_INPUTS, p, q)
        reihe = [CP.iT, CP.iHmass, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
                 CP.iconductivity, CP.iCpmass, CP.iviscosity]
        props = [fluid.keyed_output(k) for k in reihe]
        temp, h, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]

        if option == 0:
            alle = np.array([temp, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, 
                             prandtl, phase])
        elif option == 1:
            alle = [temp, p,  h, 1/rho, s,  x]
    
    return np.array(alle)

def prop_Tq(T, q, fluid="", composition=[1.0], option=1, units=_units,
            props=_props, RP=RP):
    """
    Saturation state properties at given pressure  and quality for a certain 
    fluid (mixture).

    Parameters
    ----------
    Temp in K.
    q: float (0<q<1)
        quality

    fluid :   an AbstractState in coolprop or a fluid in REFPROP.
    composition: list of floats
        mole fraction of each fluid.
    option: integer
        [1] only T p h v s and q

    units: integer
        select units (SI of property model, for REFPROP: 21)
    RP : REFPROP Instance as set with setRPfluid
        if set, the fluid name does not have to be passed again, instead flud ="" works.


    Returns
    -------
    alle : numpy array (2,4)
        includes:  properies in saturated state at given pressure p
        of liquid (0,:) and vapor(1,:),
        T, p, h, v,s, q
        all in SI units.

    """

    if props == "REFPROP":
        o = RP.REFPROP2dll(fluid, "TQ", "P;H;D;S", units, 0, T, q, composition)
        if option == 0:
            alle = []

        elif option == 1:
            alle = [T, o.Output[0], o.Output[1],
                    1 / o.Output[2], o.Output[3], q]

    elif props == "CoolProp":
        fluid.update(CP.TQ_INPUTS, T, q)
        reihe = [CP.iP, CP.iHmass, CP.iQ, CP.iSmass, CP.iDmass, 
                 CP.iPrandtl, CP.iPhase,
                 CP.iconductivity, CP.iCpmass, CP.iviscosity]
        props = [fluid.keyed_output(k) for k in reihe]
        p, h, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]

        if option == 0:
            alle = np.array([T, p,  h, 1/ rho, s,x,cp, mu,  lambda_s, 
                             prandtl, phase])
        elif option == 1:
            alle = [T, p,  h, 1/rho, s,  x]
    print(props)
    return np.array(alle)


if __name__ == "__main__":

    # _vielPrint__ = False

    _props = "REFPROP"
    x_0 = 1.
    p_0 = 1e5     # Anfangsdruck Pa

    temp_sur = 283.15
    temp_0_s = 373.15

    p_sur = 1.013e5
    if _props == "CoolProp":
        fluid_a = "n-Propane"  # Working fluid
        fluid_s = "Water"
        p_c = CP.PropsSI('Pcrit', fluid_a)
        temp_0 = CP.PropsSI('T', 'P', p_0, 'Q', x_0, fluid_a)

        h_0 = CP.PropsSI('H', 'P', p_0, 'T', temp_0 + 1, fluid_a)
        working_fluid = CP.AbstractState("BICUBIC&HEOS", fluid_a)
        # mm = ht_properties_sat(1e6, working_fluid)
        secondary_fluid = CP.AbstractState("TTSE&HEOS", fluid_s)
        h_0_s = tp(temp_0_s, p_sur,  secondary_fluid, props=_props)[2]
        h_end = CP.PropsSI('H', 'P', p_sur, 'T', temp_0_s, fluid_a)
    elif _props == "REFPROP":
        # Fluid --------------------------------
        fluid_s = "Propane * Pentane"
        wf = setRPFluid(fluid_s)
        fluid_s =""
        comp = [.4, 0.6]
        #secondary_fluid = CP.AbstractState("TTSE&HEOS", fluid_s)
        # interesting, when using "BICUBIC&HEOS" the exergy of the ambient state is 0.15!
        t0 = time()
        state_data = tp(temp_0_s, p_sur, fluid_s, composition=comp,RP=wf)
        print(state_data, time()-t0)
        print(hp(state_data[2], p_sur, fluid_s, composition=comp, RP=wf))
        
        # new Fluid:
        fluid_s = "Water"
        wf = setRPFluid(fluid_s)  # new instance !
        fluid_s =""
        t0 = time()
        alles = p_prop_sat(p_0, fluid_s, option=0, RP=wf)
        print("Water with transport", p_prop_sat(
            p_0, fluid_s, option=0, RP=wf), "\n", time() - t0)
        print("Water with error(sat)", hp(
            alles[0][2], p_sur, fluid_s, option=0, RP=wf))
        print("Water single phase without error", hp(
            alles[0][2] + 1e3, p_sur, fluid_s, option=0, RP=wf),"\n\n")
    # h_0_s = tp(temp_0_s, p_sur,  secondary_fluid, props=_props)[2]
    # h_end = CP.PropsSI('H', 'P', p_sur, 'T', temp_0_s, fluid_a)
    # ex1 = hp_exergy(h_0_s, p_sur, secondary_fluid, props=_props)
    # ex2 = hp_exergy(h_0, p_sur, working_fluid, props=_props)
    # print( "Exergies (J/kg):", ex1, ex2)
    print(p_prop_sat(p_0, fluid_s, composition=comp, props=_props, RP=wf), "\n")
    
    # Example for an alternative call:
    print("\n Alternative: ", p_prop_sat(p_0, "water", composition=comp, props=_props))
