# -*- coding: utf-8 -*-
"""
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


class pinch_point:
    def __init__(self, fluids, temps, ps,  points=30, dT_hex=0.5,
                 dT_superh=10, props="REFPROP",
                 compositions=[[1.0], [1.0]], calc_type="const", name="HEX_0",
                 units=21):

        calculate = True
        self.temps = temps
        self.ps = ps
        self.points = points
        self.dT_hex = dT_hex
        self.dT_superh = dT_superh
        self.compositions = compositions
        self.props = props
        self.units = units
        self.calc_type = calc_type
        self.name = name
        self.fluids =fluids
        if calculate:
            import sys
            from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary as modwf

            from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary as modsf
            self.RP = [flp.setRPFluid(self.fluids[0], modwf, 'RPPREFIX'),
                       flp.setRPFluid(self.fluids[1], modsf, 'RPPREFIXs')]
            state_in = np.zeros((2, 2))
            sat_v = flp.prop_pq(ps[0], 1.,
                               fluids[0], compositions[0],
                               props=props, units=units, RP=self.RP[0])
            
            sat_l = flp.prop_pq(ps[0], 0.,
                               fluids[0], compositions[0],
                               props=props, units=units, RP=self.RP[0])
            
            ev_out = flp.tp(sat_v[0]+dT_superh, ps[0],
                               fluids[0], compositions[0],
                               props=props, units=units, RP=self.RP[0])
            
            sf_in = flp.tp(sat_v[0]+dT_superh+ dT_hex, ps[1],
                               fluids[1], compositions[1],
                               props=props, units=units, RP=self.RP[1])
            sf_out = flp.tp(sat_l[0]+ dT_hex, ps[1],
                               fluids[1], compositions[1],
                               props=props, units=units, RP=self.RP[1])
            
            print(sat_v, sat_l, ev_out, "\nsf:",sf_in, sf_out)
            
            dh0 = ev_out[2] - sat_l[2]
            dh1 =sf_in[2] - sf_out[2]
            m_ratio =dh1 / dh0
            h_0 = np.linspace(sat_l[2], ev_out[2])
            h_1 = np.linspace(sf_out[2], sf_in[2])
            t0 = flp.hp_v(h_0,ps[0],
                               fluids[0], compositions[0],
                               props=props, units=units, RP=self.RP[0])
            t1 = flp.hp_v(h_1,ps[1],
                               fluids[1], compositions[1],
                               props=props, units=units, RP=self.RP[1])
            self.enthalpies =[h_0, h_1]
            self.m_ratio = m_ratio
            self.dh = [dh0,dh1]
            self.t_all =[t0,t1]
            self.points=[sat_l,sat_v, ev_out, sf_in, sf_out]
                                  
            
            


if __name__ == "__main__":
    fl1 = "Propane*Pentane*Butane"
    fl2 = "Methanol"
    compositions = [[.25, .3,.45], [1.]]
    Ts = [290., 350.]
    # fl = fl_names
    flx = [fl1, fl2]
    ps = [1.2e5,  12e5]
    hp0 = pinch_point(flx, Ts, ps, dT_superh=15, compositions=compositions)
    plt.figure()
    plt.plot(hp0.t_all[0][2,:]- hp0.t_all[0][2,0], hp0.t_all[0][0,:])
    plt.plot((hp0.t_all[1][2,:] -hp0.t_all[1][2,0])/hp0.m_ratio, hp0.t_all[1][0,:],"o")
