# -*- coding: utf-8 -*-
"""
Property evaluation from actual Refprop dll for the compressor model
of Dennis Roskosch

Created on Fri Dec  9 21:26:04 2022

@author: atakan
"""

import numpy as np
import carbatpy.fluid_properties_rp as fp
_props = "REFPROP" 
fluid_ = "Propane * Butane"
comp_ = [1.0, 0.]  # [0.5, 0.5]

def z_uv(u_, v, fluid=fluid_ , comp = comp_, props=_props):
    u =u_ * 1000.
    T, p, h,v,s,q = fp.uv(u, v, fluid, comp,props=props)   
    return np.array([T, p/1000., v, u/1000., h/1000., s])

def z_ps(p_, s, fluid=fluid_, comp = comp_, props=_props):
    p = p_ * 1000.
    T, p, h,v,s,q = fp.sp(s, p,fluid, comp,props=props)   
    u = h - p * v
    
    return np.array([T, p/1000., v, u/1000., h/1000., s])

def z_Tp(T, p_, fluid=fluid_, comp = comp_, props=_props):
    p=p_ * 1000.
    T, p, h,v,s,q = fp.tp(T,p, fluid, comp,props=props) 
    u = h - p * v
   
    return np.array([T, p/1000., v, u/1000., h/1000., s])


def z_Tx(T, x, fluid=fluid_, comp = comp_, props=_props):
    T, p, h,v,s,q  = fp.prop_Tq(T, x, fluid, comp, props=props)
    u = h - p * v
    return np.array([T, p/1000., v, u/1000., h/1000., s])


def z_px(p, x, fluid=fluid_, comp = comp_, props=_props):
    p = p * 1000.
    T, p, h,v,s,q  = fp.prop_pq(p, x, fluid, comp, props=props)
    u = h - p * v
    return np.array([T, p/1000., v, u/1000., h/1000., s])

def z_mm(T, p_, fluid=fluid_, comp = comp_, props=_props):
    p=p_ * 1000.
    alle= fp.tp(T,p, fluid=fluid, composition=comp,props=props, option=2) 
    T, p, h,v,s,q = alle[:6] 
    molecular_mass_mix = alle[-1]
    u = h - p * v
   
    return np.array([T, p/1000., v, u/1000., h/1000., s, molecular_mass_mix])

if __name__ == "__main__":
    state = z_Tx(300, 0, fluid_, comp_)
    state2 =z_uv(state[3],state[2], fluid_, comp_)
    state3 =z_mm(300,100, fluid_, comp_)
    state3 =z_ps(state[1],state[5], fluid_, comp_)
    print(state,"\n", state2, "\n", state3)