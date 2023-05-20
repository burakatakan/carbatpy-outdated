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


class static_heat_exchanger:
    def __init__(self, fluids, temps, ps, qs=[0, -2], points=30, dT_hex=0.5,
                 dT_superh=10, heating=True, dH_min=0, props="REFPROP",
                 compositions=[[1.0], [1.0]], calc_type="const",
                 name="evaporator",
                 units=21):  # qs def.: <0=> liquid entering, qs >1 => vapor entering

        self.temps = temps
        self.ps = ps
        self.qs = qs
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

        from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary as modwf

        from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary as modsf
        self.RP = [flp.setRPFluid(self.fluids[0], modwf, 'RPPREFIX'),
                   flp.setRPFluid(self.fluids[1], modsf, 'RPPREFIXs')]
        sat_p_initial = self.ps[0]
        sat_p = flp.prop_Tq(self.temps[0], 1.,
                            self.fluids[0], self.compositions[0],
                            props=self.props, units=self.units,
                            RP=self.RP[0])
        
        if self.heating:
            self.ps[0] = sat_p[1]
            print(
                f"Saturation pressure varied from {sat_p_initial:.2e} to {self.ps[0] :.2e} Pa!")
            sat_v = flp.prop_pq(self.ps[0], 1.,
                                self.fluids[0], self.compositions[0],
                                props=self.props, units=self.units,
                                RP=self.RP[0])
            if self.qs[0] < 0 or self.qs[0] > 1:
                raise Exception(f"working fluid quality is wrong!{self.qs}!")

            sat_l = flp.prop_pq(self.ps[0], self.qs[0],
                                self.fluids[0], self.compositions[0],
                                props=self.props, units=self.units,
                                RP=self.RP[0])

            ev_out = flp.tp(sat_v[0]+self.dT_superh, self.ps[0],
                            self.fluids[0], self.compositions[0],
                            props=self.props, units=self.units,
                            RP=self.RP[0])

            sf_in = flp.tp(sat_v[0]+self.dT_superh + self.dT_hex, self.ps[1],
                           self.fluids[1], self.compositions[1],
                           props=self.props, units=self.units,
                           RP=self.RP[1])
            sf_out = flp.tp(sat_l[0] + self.dT_hex, self.ps[1],
                            self.fluids[1], self.compositions[1],
                            props=self.props, units=self.units, RP=self.RP[1])

            if verbose:
                print(sat_v, sat_l, ev_out, "\nsf:", sf_in, sf_out)

            dh0 = ev_out[2] - sat_l[2]
            dh1 = sf_in[2] - sf_out[2]
            m_ratio = dh1 / dh0
            h_0 = np.linspace(sat_l[2], ev_out[2], self.points)
            h_1 = np.linspace(sf_out[2], sf_in[2], self.points)
            t0 = flp.hp_v(h_0, self.ps[0],
                          self.fluids[0], self.compositions[0],
                          props=self.props, units=self.units, RP=self.RP[0])
            t1 = flp.hp_v(h_1, self.ps[1],
                          self.fluids[1], self.compositions[1],
                          props=self.props, units=self.units, RP=self.RP[1])
            self.enthalpies = [h_0, h_1]
            self.m_ratio = m_ratio
            self.dh = [dh0, dh1]
            self.t_all = [t0, t1]
            self.points = [sat_l, sat_v, ev_out, sf_in, sf_out]
            h0_shifted = hp0.t_all[0][2, :] - \
                hp0.t_all[0][2, 0], hp0.t_all[0][0, :]
            h1_shifted = (hp0.t_all[1][2, :] - hp0.t_all[1]
                          [2, 0])/hp0.m_ratio, hp0.t_all[1][0, :]
            dT_min = np.min(h1_shifted[1]-h0_shifted[1])
            if dT_min < self.dT_hex*.999:

                print(f"Problem{dT_min}")
        else:
            con_in = flp.tp(self.temps[0]+self.dT_superh+self.dT_hex, self.ps[0],
                            self.fluids[0], self.compositions[0],
                            props=self.props, units=self.units,
                            RP=self.RP[0])
            con_dew = flp.prop_pq(self.ps[0], 1.,
                                 self.fluids[0], self.compositions[0],
                                 props=self.props, units=self.units,
                                 RP=self.RP[0])
            dew_T = con_dew[0] # here it must be prooved whether the pinch point T-difference is ok, if not: pressure change
            sf_out = flp.tp(self.temps[0]+self.dT_superh, self.ps[1],
                            self.fluids[1], self.compositions[1],
                            props=self.props, units=self.units, RP=self.RP[1])
            sf_dew = flp.tp(dew_T-self.dT_hex, self.ps[1],
                            self.fluids[1], self.compositions[1],
                            props=self.props, units=self.units, RP=self.RP[1])
            m_ratio =(sf_out[2]-sf_dew[2])/(con_in[2]-con_dew[2])
            
            if self.dH_min ==0:
                   
                sf_in = flp.tp(self.temps[1], self.ps[1],
                               self.fluids[1], self.compositions[1],
                               props=self.props, units=self.units, RP=self.RP[1])
                dh1 = sf_out[2] -sf_in[2]
                dh0 = dh1 / m_ratio
                con_out = flp.hp(con_in[2]-dh0, self.ps[0],
                                self.fluids[0], self.compositions[0],
                                props=self.props, units=self.units,
                                RP=self.RP[0])
            else:
                dh_remain = np.abs(self.dH_min)-(con_in[2]-con_dew[2])
                con_out = flp.hp(con_in[2]-self.dH_min, self.ps[0],
                                self.fluids[0], self.compositions[0],
                                props=self.props, units=self.units,
                                RP=self.RP[0])
                sf_in = flp.tp(con_out[0]-self.dT_hex, self.ps[1],
                               self.fluids[1], self.compositions[1],
                               props=self.props, units=self.units, RP=self.RP[1])
                m_ratio = (sf_dew[2]-sf_in[2])/(dh_remain)
                sf_out = flp.hp(sf_in[2]+self.dH_min*m_ratio, self.ps[1],
                               self.fluids[1], self.compositions[1],
                               props=self.props, units=self.units, RP=self.RP[1])
                dh0 =self.dH_min
                dh1 =dh0 / m_ratio # BA hier pruegfen!!
                
                                            
            h_0 = np.linspace(con_out[2], con_in[2], self.points)
            h_1 = np.linspace(sf_in[2], sf_out[2], self.points)
            t0 = flp.hp_v(h_0, self.ps[0],
                          self.fluids[0], self.compositions[0],
                          props=self.props, units=self.units, RP=self.RP[0])
            t1 = flp.hp_v(h_1, self.ps[1],
                          self.fluids[1], self.compositions[1],
                          props=self.props, units=self.units, RP=self.RP[1])
            self.enthalpies = [h_0, h_1]
            self.m_ratio = m_ratio
            self.dh = [dh0, dh1]
            self.t_all = [t0, t1]
            self.points = [con_in, con_out, con_dew,sf_dew, sf_in, sf_out]
            h0_shifted = hp0.t_all[0][2, :] - \
                hp0.t_all[0][2, 0], hp0.t_all[0][0, :]
            h1_shifted = (hp0.t_all[1][2, :] - hp0.t_all[1]
                          [2, 0])/hp0.m_ratio, hp0.t_all[1][0, :]
            dT_min = np.min(h0_shifted[1]-h1_shifted[1])
            if dT_min < self.dT_hex*.99:

                print(f"Problem:{dT_min}")

    def pp_root(self, var_in):
        if self.heating:
            if self.qs[0] >= 0 and self.qs[0] <= 1:
                pass


if __name__ == "__main__":
    fl1 = "Propane*Pentane*Butane"
    fl2 = "Methanol"
    compositions = [[.3, .4, .3], [1.]]
    Ts = [290., 250.]
    # fl = fl_names
    flx = [fl1, fl2]
    ps = [1.92e5,  12e5]
    qs = [0.05, -2]

    hp0 = static_heat_exchanger(
        flx, Ts, ps, dT_superh=15, qs=qs, compositions=compositions)
    hp0.pinchpoint()

    plt.figure()
    shift1 = hp0.t_all[0][2, -1] - hp0.t_all[0][2, 0]
    
    plt.plot(hp0.t_all[0][2, :] - hp0.t_all[0][2, 0], hp0.t_all[0][0, :])
    plt.plot((hp0.t_all[1][2, :] - hp0.t_all[1][2, 0]) /
             hp0.m_ratio, hp0.t_all[1][0, :], "o")
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
