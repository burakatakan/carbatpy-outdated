# -*- coding: utf-8 -*-
"""
Created on Sun May 21 10:52:09 2023

@author: atakan
"""



import CoolProp.CoolProp as CP

import fluid_properties_rp as fprop
import numpy as np
import components.compressor_simple as comp
import components.throttle_simple as throt
import components.heat_exchanger_thermo as hext


working_fluid = "Propane * Butane * Pentane"
x1 =0.2
x2 = 0.4
x3 =1-x1-x2
composition_wf =[x1,x2,x3]
p_levels_wf = [1.1e5, 5.5e5]
T_pass =300
T_sat_low_p = fprop.prop_pq(p_levels_wf[0], 1, 
                      working_fluid, composition_wf)
T_sat_high_p = fprop.prop_pq(p_levels_wf[1], 1, 
                      working_fluid, composition_wf)
pass_low_p = fprop.tp(T_pass, p_levels_wf[0],
                      working_fluid, composition_wf)
pass_high_p = fprop.tp(T_pass, p_levels_wf[1],
                      working_fluid, composition_wf)

print(f"T-levels for this mixture {T_sat_low_p } K and {T_sat_high_p } K.")

eta_c = .65
state_in = fprop.tp(*pass_low_p[:2], working_fluid, composition_wf)
state_c_out = comp.compressor(state_in, p_levels_wf[1], eta_c, 
                              working_fluid, composition=composition_wf)
w_compressor =state_c_out[2]-state_in[2]

print(f"\nState after compresor: {state_c_out}, {state_c_out[0]-273.15} C")
dh_condenser = state_c_out[2] - pass_high_p[2]


Ts = [T_sat_high_p[0], T_pass]
dT_sup = state_c_out[0]-T_sat_high_p[0]
if dT_sup< 5: print(f"Warning superheating too low: {dT_sup} K")
ps = [state_c_out[1],  12e5]
fl2 = "Water"
flx = [working_fluid, fl2]
qs = [2, -2]
compositions=[composition_wf,[1.]]

hp_condenser = hext.static_heat_exchanger(
    flx, Ts, ps, dT_hex=1, h_enter=[state_c_out[2], -1e-9], dT_superh=dT_sup, 
    qs=qs,heating =False, dH_min=dh_condenser,
    compositions=compositions)
hp_condenser.pinchpoint()
hp_condenser.hex_plot()
q_condenser =hp_condenser.dh[0]

throttle_exit = throt.throttle(hp_condenser.t_all[0][:,0], p_levels_wf[0],
                               working_fluid, composition=composition_wf)
flb =[working_fluid, "Methanol"]
dh_ev = state_in[2]-throttle_exit[2]
qs =[throttle_exit[-1], -2]
ps[0] = p_levels_wf[0]
print("evaporator")
evaporator = hext.static_heat_exchanger(
    flb, Ts, ps, dT_hex=2, h_enter=[throttle_exit[2], -1e-9], dT_superh=dT_sup, 
    qs=qs,heating =True, dH_min=dh_ev,
    compositions=compositions)
evaporator.pinchpoint()
evaporator.hex_plot(hp_condenser)
# #storage_fluid_ht = 