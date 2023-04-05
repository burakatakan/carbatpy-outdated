# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:30:19 2022

definition of class to calculate heat transfer in two-phase region
@author: welp
"""

import CoolProp.CoolProp as CP
import scipy
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
#MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root


class BoundaryConditions:
    def __init__(self):
        self.Ra = 1.5e-6
        self.q_dot = 9000
        self.inlet_flow = 'true'
        self.RB = 'qzu'
        self.Rz = 7 * self.Ra
        self.g = 9.81
        self.lam_w = 320
        self.s = 2e-3


class PointND(BoundaryConditions):
    def __init__(self, T: float, q: float, fluid: str, name: str):
        BoundaryConditions.__init__(self)
        self.fluid = fluid
        self.T = T
        self.q = q
        fluid = 'REFPROP::' + fluid
        self.name = name
        self.rho = CP.PropsSI('D', 'T', T, 'Q', q, fluid)
        self.h = CP.PropsSI('H', 'T', T, 'Q', q, fluid)
        self.p = CP.PropsSI('P', 'T', T, 'Q', q, fluid)
        self.pcrit = CP.PropsSI('Pcrit', fluid)
        self.M = CP.PropsSI('molarmass', fluid)
        self.sigma = CP.PropsSI('surface_tension', 'T', T, 'Q', q, fluid)

        # saturation properties
        # v : vapour, l : liquid
        self.rho_v = CP.PropsSI('D', 'T', T, 'Q', 1, fluid)
        self.rho_l = CP.PropsSI('D', 'T', T, 'Q', 0, fluid)
        self.h_v = CP.PropsSI('H', 'T', T, 'Q', 1, fluid)
        self.h_l = CP.PropsSI('H', 'T', T, 'Q', 0, fluid)
        self.h_evap = self.h_v - self.h_l
        self.dyn_vis_v = CP.PropsSI('V', 'T', T, 'Q', 1, fluid)
        self.dyn_vis_l = CP.PropsSI('V', 'T', T, 'Q', 0, fluid)
        self.kin_vis_v = self.dyn_vis_v / self.rho_v
        self.kin_vis_l = self.dyn_vis_l / self.rho_l
        self.lam_l = CP.PropsSI('L', 'T', T, 'Q', 0, fluid)
        self.lam_v = CP.PropsSI('L', 'T', T, 'Q', 1, fluid)
        self.pr_l = CP.PropsSI('Prandtl', 'T', T, 'Q', 0, fluid)
        self.pr_v = CP.PropsSI('Prandtl', 'T', T, 'Q', 1, fluid)
        self.c_v = CP.PropsSI('SPEED_OF_SOUND', 'T', T, 'Q', 1, fluid)
        self.cp_v = CP.PropsSI('CP0MASS', 'T', T, 'Q', 1, fluid)
        self.cp_l = CP.PropsSI('CP0MASS', 'T', T, 'Q', 0, fluid)

    def reynolds_v(self, m_dot: float, d: float) -> float:
        """
        Reynolds number in vapour state
        :param m_dot: mass flow rate per area [kg/s/m**2]
        :param d: diameter [m]
        :return: Reynolds number in vapour state
        """
        Re_v = m_dot * self.q * d / self.dyn_vis_v  # H3.1 2.2 (11)
        return Re_v

    def reynolds_l(self, m_dot: float, d: float) -> float:
        """
        Reynolds number in liquid state
        :param m_dot: mass flow rate per area [kg/s/m**2]
        :param d: diameter [m]
        :return: Reynolds number in liquid state
        """
        re_l = m_dot * (1 - self.q) * d / self.dyn_vis_l  # H3.1 2.1 (6)
        return re_l

    def u_rohr(self, m_dot: float) -> float:
        """
        velocity in tube
        :param m_dot: mass flow rate per area [kg/s/m**2]
        :return: velocity in tube [m/s]
        """
        u = m_dot / self.rho
        return u

    def froude(self, m_dot: float, d: float) -> float:
        """
        Froude number
        :param m_dot: mass flow rate per area [kg/s/m**2]
        :param d: diameter [m]
        :return: Froude number [-]
        """
        fr = m_dot ** 2 * self.q ** 2 / (self.rho_l * self.rho_v * g * d)
        return fr

    def froude_l(self, m_dot: float, d: float) -> float:
        """
        Froude number in liquid state
        :param m_dot: overall mass flow rate per area [kg/s/m**2]
        :param d: diameter [m]
        :return: Froude number for liquid state [-]
        """
        fr_l = m_dot ** 2 * (1 - self.q) ** 2 / (self.rho_l ** 2 * g * d)
        return fr_l

    def froude_v(self, m_dot: float, d: float) -> float:
        """
        Froude number in vapour state
        :param m_dot: overall mass flow rate per area [kg/s/m^2]
        :param d: diameter [m]
        :return: Froude number for vapour state [-]
        """
        fr_v = m_dot ** 2 * self.q ** 2 / (self.rho_v ** 2 * g * d)
        return fr_v

    def alpha_1p_l(self, m_dot: float, d: float, z_i: float) -> float:
        """
        calculates heat transfer under assumption that whole flow is liquid
        :param m_dot: overall mass flow rate per area [kg/s/m^2]
        :param d: diameter [m]
        :param z_i: local variable in tube [m]
        :return: heat transfer coefficient [W/m^2*K]
        """
        if z_i == 0:
            z_i = 1e-5
        if d / z_i >= 1:
            dz = 1
        else:
            dz = d / z_i

        re_l0 = m_dot * d / self.dyn_vis_l
        xi = (1.82 * np.log10(re_l0) - 1.64) ** (-2)
        nu_inf = xi / 8 * (re_l0 - 1000) * self.pr_l / (1 + 12.7 * (xi / 8) ** 0.5 * (self.pr_l ** (2 / 3) - 1))

        # laminar case:
        nu_x_lam = 0
        if re_l0 < 5 * 10 ** 4:

            if self.RB == 'Tw':
                if self.inlet_flow == 'true':
                    nu_x_lam = 0.332 * self.pr_l ** (1 / 3) * (re_l0 * dz) ** (1 / 2)
                if self.inlet_flow == 'false':
                    nu_x_lam = (3.66 ** 3 + 1.077 ** 3 * re_l0 * self.pr_l * dz) ** (1 / 3)

            if self.RB == 'qzu':
                if self.inlet_flow == 'true':
                    nu_x_lam = 0.455 * self.pr_l ** (1 / 3) * (re_l0 * dz) ** 0.5
                if self.inlet_flow == 'false':
                    nu_x_lam = (4.36 ** 3 + 1.302 ** 3 * re_l0 * self.pr_l * dz) ** (1 / 3)

        # turbulent case:
        nu_x_turb = 0
        if re_l0 >= 2300:
            if self.inlet_flow == 'true':
                if d / z_i >= 1:
                    nu_x_turb = 4 / 3 * nu_inf
                else:
                    nu_x_turb = nu_inf * (1 + 1 / 3 * dz ** (2 / 3))
            if self.inlet_flow == 'false':
                nu_x_turb = nu_inf

        if nu_x_turb > nu_x_lam:
            nu_x = nu_x_turb
        else:
            nu_x = nu_x_lam
        alpha = nu_x * self.lam_l / d
        return alpha

    def alpha_1p_v(self, m_dot: float, d: float, z_i: float) -> float:
        """
        calculates heat transfer under assumption that whole flow is gaseous
        :param m_dot: overall mass flow rate per area [kg/s/m^2]
        :param d: diameter [m]
        :param z_i: z_i: local variable in tube [m]
        :return: heat transfer coefficient [W/m^2*K]
        """
        if z_i == 0:
            z_i = 1e-5
        if d / z_i >= 1:
            dz = 1
        else:
            dz = d / z_i

        re_v0 = m_dot * d / self.dyn_vis_v
        xi = (1.82 * np.log10(re_v0) - 1.64) ** (-2)
        nu_inf = xi / 8 * (re_v0 - 1000) * self.pr_v \
                 / (1 + 12.7 * (xi / 8) ** 0.5 * (self.pr_v ** (2 / 3) - 1))

        # laminar case:
        nu_x_lam = 0
        if re_v0 < 5 * 10 ** 4:

            if RB == 'Tw':
                if inlet_flow == 'true':
                    nu_x_lam = 0.332 * self.pr_v ** (1 / 3) * (re_v0 * dz) ** (1 / 2)
                if inlet_flow == 'false':
                    nu_x_lam = (3.66 ** 3 + 1.077 ** 3 * re_v0 * self.pr_v * dz) ** (1 / 3)

            if RB == 'qzu':
                if inlet_flow == 'true':
                    nu_x_lam = 0.455 * self.pr_v ** (1 / 3) * (re_v0 * dz) ** 0.5
                if inlet_flow == 'false':
                    nu_x_lam = (4.36 ** 3 + 1.302 ** 3 * re_v0 * self.pr_v * dz) ** (1 / 3)

        # turbulent case:
        nu_x_turb = 0
        if re_v0 >= 2300:
            if self.inlet_flow == 'true':
                if d / z_i >= 1:
                    nu_x_turb = 4 / 3 * nu_inf
                else:
                    nu_x_turb = nu_inf * (1 + 1 / 3 * dz ** (2 / 3))
            if self.inlet_flow == 'false':
                nu_x_turb = nu_inf

        if nu_x_turb > nu_x_lam:
            nu_x = nu_x_turb
        else:
            nu_x = nu_x_lam
        alpha = nu_x * self.lam_v / d
        return alpha

    def alpha_g_cal(self, m_dot: float, d: float, z_i: float, dhg: float) -> float:
        """

        :param m_dot:
        :param d:
        :param z_i:
        :param dhg:
        :return:
        """
        if z_i == 0:
            z_i = 1e-5
        if dhg / z_i >= 1:
            dz = 1
        else:
            dz = dhg / z_i

        ReG = m_dot * self.q * dhg / self.dyn_vis_v / eps
        xi = (1.82 * np.log10(ReG) - 1.64) ** (-2)
        Nu_inf = xi / 8 * (ReG - 1000) * self.pr_v \
                 / (1 + 12.7 * (xi / 8) ** 0.5 * (self.pr_v ** (2 / 3) - 1))

        # laminar case:
        Nu_x_lam = 0
        if ReG < 5 * 10 ** 4:
            if RB == 'Tw':
                if inlet_flow == 'true':
                    Nu_x_lam = 0.332 * self.pr_v ** (1 / 3) * (ReG * dz) ** (1 / 2)
                if inlet_flow == 'false':
                    Nu_x_lam = (3.66 ** 3 + 1.077 ** 3 * ReG * self.pr_v * dz) ** (1 / 3)

            if RB == 'qzu':
                if inlet_flow == 'true':
                    Nu_x_lam = 0.455 * self.pr_v ** (1 / 3) * (ReG * dz) ** 0.5
                if inlet_flow == 'false':
                    Nu_x_lam = (4.36 ** 3 + 1.302 ** 3 * ReG * self.pr_v * dz) ** (1 / 3)

        # turbulent case:
        Nu_x_turb = 0
        if ReG >= 2300:
            if inlet_flow == 'true':
                if d / z_i >= 1:
                    Nu_x_turb = 4 / 3 * Nu_inf
                else:
                    Nu_x_turb = Nu_inf * (1 + 1 / 3 * (dz) ** (2 / 3))
            if inlet_flow == 'false':
                Nu_x_turb = Nu_inf

        if Nu_x_turb > Nu_x_lam:
            Nu_x = Nu_x_turb
        else:
            Nu_x = Nu_x_lam
        alpha = Nu_x * self.lam_v / d
        print(f"alpha_Gas: {alpha}")
        return alpha

    # TODO: create testcases, then merge three similar functions above, differ only in Re
    def alpha_sieden(self, m_dot, d, z_i):
        """
        calculates heat transfer coefficient in two-phase region, uses nucleate boiling and convective boiling
        :param m_dot: mass flow rate density, kg/(s*m**2)
        :param d: diameter, m
        :param z_i: axial position in tube, m
        :return: alpha, W/(m*K)
        """

        # step 1: does nucleate boiling need to be considered?
        r_cr = 0.3e-6  # [m]

        alpha_LO_start = self.alpha_1p_l(m_dot, d,
                                         0)  # lokale einphasige Wärmeübergangskoeffizient an der Stelle z=0 für Flüssigkeit
        print(f"alpha_LO_start: {alpha_LO_start}")

        q_dot_onb = 2 * self.sigma * self.T * alpha_LO_start / (r_cr * self.rho_v * self.h_evap)
        # print("q_onb: ", q_dot_onb)

        # q_dot_m = m_dot1 * h_v / (np.pi * d * l)
        # print("q_dot_m: ", q_dot_m)

        zustand, phi, eps = self.flow_character(m_dot, d)

        if self.q_dot < q_dot_onb:
            alpha_b = 0
        else:
            M_H2 = CP.PropsSI('molarmass', 'Hydrogen')
            p_stern = self.p / self.pcrit
            n = 0.9 - 0.36 * p_stern ** 0.13
            CF = 0.789 * (self.M / M_H2) ** 0.11
            if (zustand == 'Wellenströmung' or zustand == 'Schichtenströmung') and (self.RB == 'Tw'):
                CF = CF * 0.86 * self.lam_w * self.s
            if self.RB == 'qzu' and (self.lam_w * self.s) < 0.7:
                kappa = 0.675 + 0.325 * np.tanh(3.711 * (self.lam_w * self.s - 3.24e-2))  # H3.4 2.2 (33a)
                n = kappa * n  # H3.4 2.2 (34a)
                if zustand == 'Wellenströmung' or zustand == 'Schichtenströmung':  # P: kleines psi
                    p = 0.46 + 0.4 * np.tanh(3.387 * (self.lam_w * self.s - 8.62e-3))  # H3.4 2.2 (33b)
                if zustand == 'Pfropfenströmung':
                    p = 0.671 + 0.329 * np.tanh(3.691 * (self.lam_w * self.s - 8.42e-3))  # H3.4 2.2 (33c)
                if zustand == 'Ringströmung':
                    p = 0.755 + 0.245 * np.tanh(3.702 * (self.lam_w * self.s - 1.25e-2))  # H3.4 2.2 (33d)
                CF = CF * p  # H3.4 2.2 (35)

            table2 = np.array(pd.read_excel("fluid_prop_VDI_Tab2.xlsx"))

            for i in range(len(table2)):
                q0 = 0
                alpha_0 = 0
                if table2[i, 0] == self.fluid:
                    q0 = table2[i, 5]
                    alpha_0 = table2[i, 6]
                    break
            if self.fluid == "R12":
                q0 = 20000
                alpha_0 = 3290
            if q0 == 0 and alpha_0 == 0:
                print('Fluid ', self.fluid, ' not found.')
                del q0, alpha_0

            ### subpart for p_stern = 0.1 ###
            p_stern_tem = 0.1
            p_tem = p_stern_tem * self.pcrit
            T_tem = CP.PropsSI('T', 'P', p_tem, 'Q', self.q, self.fluid)
            point_tem = PointND(T_tem, self.q, self.fluid, 'Hilfspunkt für p_stern_tem = 0.1')

            q_dot_kr = 0.144 * point_tem.h_evap * ((point_tem.rho_l \
                                                    - point_tem.rho_v) * point_tem.rho_v) ** 0.5 \
                       * ((self.g * point_tem.sigma) / point_tem.rho_l) ** 0.25 \
                       * point_tem.pr_l ** (-0.245)  # TODO: which pradtl number has to be used here?

            del (point_tem)

            ### end subpart ###

            if p_stern >= 0.1:
                q_cr_PB = q_dot_kr * 3.2 * p_stern ** 0.45 * (1 - p_stern) ** 1.2
            if p_stern < 0.1:
                q_cr_PB = q_dot_kr * 1.2 * (p_stern ** 0.17 + p_stern ** 0.8)

            alpha_b = alpha_0 * CF * (self.q_dot / q0) ** n \
                      * (2.692 * p_stern ** 0.43 + 1.6 * p_stern ** 6.5 / (1 - p_stern ** 4.4)) \
                      * (1e-2 / d) ** 0.5 * (self.Ra / 1e-6) ** 0.133 * (m_dot / 100) ** 0.25 \
                      * (1 - p_stern ** 0.1 * (self.q_dot / q_cr_PB) ** 0.3 * self.q)  # H3.4 2.2 (31)

            alpha_b_hzw = alpha_0 * CF * (self.q_dot / q0) ** n \
                          * (2.816 * p_stern ** 0.45 + (3.4 + 1.7 / (1 - p_stern ** 7)) * p_stern ** 3.7) \
                          * (1e-2 / d) ** 0.4 * (self.Ra / 1e-6) * 0.133  # H3.4 2.1 (21)

            if alpha_b_hzw < alpha_b:
                alpha_b = alpha_b_hzw

        ### step 2: convective boiling
        alpha_LO = self.alpha_1p_l(m_dot, d, z_i)
        # print(f"alpha_LO: {alpha_LO}")
        alpha_GO = self.alpha_1p_v(m_dot, d, z_i)
        # print(f"alpha_GO: {alpha_GO}")
        # print("alpha_GO :", alpha_GO)
        alpha_quotient = ((1 - self.q) ** 0.01 * ((1 - self.q) + 1.2 * self.q ** 0.4 \
                                                  * (self.rho_l / self.rho_v) ** 0.37) ** (-2.2) + self.q ** 0.01 \
                          * (alpha_GO / alpha_LO * (1 + 8 * (1 - self.q) ** 0.7 \
                                                    * (self.rho_l / self.rho_v) ** 0.67)) ** (-2)) ** (-0.5)
        # print(f"alpha_quotient: {alpha_quotient}")
        alpha_k = alpha_LO * alpha_quotient

        if zustand == 'Schichtenströmung':

            PHI = 0.5 * phi / (2 * np.pi)
            dm = d + s
            dhg = d * (phi - np.sin(180 / np.pi * phi)) / (phi + 2 * np.sin(90 / np.pi * phi))

            alpha_g = self.alpha_g_cal(m_dot, d, z_i, dhg)
            a = alpha_g / alpha_k

            M = alpha_k * d * np.pi ** 2 * dm / (lam_w * 4 * s)

            if M <= 1:
                alpha_k = alpha_k * (1 - PHI) + alpha_g * PHI

            else:
                f1 = a * M ** 0.5 * PHI * (1 - PHI) * 1 / (np.tanh(M ** 0.5 * (1 - PHI)))
                f2 = a * M ** 0.5 * PHI * (1 - PHI) * 1 / np.tanh((a * M) ** 0.5 * PHI)
                PSI = 1 + PHI * (1 - PHI) * (1 - a) ** 2 / a * (1 - ((1 - (1 - a) * PHI) / (f1 + f2)))
                alpha_k = alpha_k * (1 - (1 - a) * PHI / PSI)

        print(f"x: {self.q}")
        print("alpha_b: ", alpha_b)
        print("alpha_k: ", alpha_k)

        alpha = (alpha_b ** 3 + alpha_k ** 3) ** (1 / 3)
        print("alpha: ", alpha)

        return alpha, alpha_b, alpha_k

    # alpha_sieden = np.vectorize(alpha_sieden)

    def flow_values(self, m_dot, d):
        """
        calculates values to determine flow character

        Parameters
        ----------
        m_dot : TYPE float.
            DESCRIPTION.  mass flow rate per area [kg/s/m**2]
        d : TYPE float.
            DESCRIPTION. diamater [m]

        Returns
        -------
        X : TYPE float.
            DESCRIPTION. Martinelli Parameter [-]
        ReFr05 : TYPE float.
            DESCRIPTION. sqrt(Reynoldsnumber (liquid) * Froudenumber (vapour)) 
        Fr05 : TYPE float.
            DESCRIPTION. sqrt(Froudenumber (vapour, average))
        FrEu05 : TYPE float.
            DESCRIPTION. sqrt(Froudenumber * Eulernumber) (liquid)
        WeFrL : TYPE float.
            DESCRIPTION. Webernumber/Froudenumber (liquid)

        """
        if self.q == 0:
            X = 1e5
        else:
            # print(f"In Schleife: q: {self.q}, rho_v: {self.rho_v}, rho_l: {self.rho_l} dyn_vis_l: {self.dyn_vis_l}, dyn_vis_v: {self.dyn_vis_v}")
            X = ((1 - self.q) / self.q) ** 0.875 * \
                (self.rho_v / self.rho_l) ** 0.5 * (
                        self.dyn_vis_l / self.dyn_vis_v) ** 0.125  # Martinelli-Parameter, H3.1 2.1 (1)
        ReFr05 = (m_dot ** 3 * self.q ** 2 * (1 - self.q) / (self.rho_v \
                                                             * (self.rho_l - self.rho_v) * self.dyn_vis_l * self.g)) ** 0.5  # H3.1 2.1 (2)
        Fr05 = (m_dot ** 2 * self.q ** 2 / (self.g * d * self.rho_l * self.rho_v)) ** 0.5  # H3.1 2.1 (3)
        xiL = 0.3164 / self.reynolds_l(m_dot, d) ** 0.25  # H3.1 2.1 (7)
        FrEu05 = (xiL * m_dot ** 2 * (1 - self.q) ** 2 / (2 * d * self.rho_l \
                                                          * (self.rho_l - self.rho_v) * self.g)) ** 0.5  # H3.1 2.1 (4)
        WeFrL = self.g * d ** 2 * self.rho_l / self.sigma  # H3.1 2.1 (5)

        return X, ReFr05, Fr05, FrEu05, WeFrL

    def epsilon(self, m_dot):
        """
        steam volume fraction after Rouhani

        Parameters
        ----------
        m_dot : TYPE float.
            DESCRIPTION. mass flow rate per area [kg/s/m**2]

        Returns
        -------
        eps : TYPE float.
            DESCRIPTION. steam volume fraction [-]

        """

        eps = self.q / self.rho_v * ((1 + 0.12 * (1 - self.q)) \
                                     * (self.q / self.rho_v + (1 - self.q) / self.rho_l) \
                                     + 1.18 * (1 - self.q) * (self.g * self.sigma * (self.rho_l - self.rho_v)) ** 0.25 \
                                     / (m_dot * self.rho_l ** 0.5)) ** (-1)  # H3.1 2.2 (26)
        return eps

    def hlstart(self, eps, phi):
        """
        initial value for heigth of liquid level

        Parameters
        ----------
        eps : TYPE float.
            DESCRIPTION. vapour volume fraction [-]
        phi : TYPE float.
            DESCRIPTION. vapour angle in tube [rad]

        Returns
        -------
        hl0 : TYPE float.
            DESCRIPTION. initial value for relative liquid level heigth

        """
        hl0 = 15 * np.pi * (1 - eps) / (8 * (3 * np.sin(phi / 2) + 4 * np.sin(phi / 4)))  # H3.1 2.2 (28)
        return hl0

    def hl_values(self, hl):
        """
        geometric values for determination of flow character (sketch p. 901, H3.1 Abb4.)

        Parameters
        ----------
        hl : TYPE float.
            DESCRIPTION. relative height of liquid level [-]

        Returns
        -------
        Ui : TYPE float.
            DESCRIPTION. relative width of liquid level [-]
        UL : TYPE float.
            DESCRIPTION. relative circumference of liquid [-]
        UG : TYPE float.
            DESCRIPTION. relative circumference of vapour [-]
        fL : TYPE float.
            DESCRIPTION. relative cross-sectional area of liquid [-]
        fG : TYPE float.
            DESCRIPTION. relative cross-sectional area of vapour [-]

        """
        Ui = 2 * (hl * (1 - hl)) ** 0.5  # H3.1 2.2 (15)

        if hl <= 0.5:
            psi = 2 * np.arcsin(Ui)  # H3.1 2.2 (16)
            UL = psi / 2  # H3.1 2.2 (17)
            UG = np.pi - UL  # H3.1 2.2 (18)
            fL = (psi - np.sin(psi)) / 8  # H3.1 2.2 (19)   removed factor in sinus
            fG = np.pi / 4 - fL  # H3.1 2.2 (20)
        if hl > 0.5:
            phi = 2 * np.arcsin(Ui)  # H3.1 2.2 (21)
            UG = phi / 2  # H3.1 2.2 (22)
            UL = np.pi - UG  # H3.1 2.2 (23)
            fG = (phi - np.sin(phi)) / 8  # H3.1 2.2 (24)  removed factor in sinus
            fL = np.pi / 4 - fG  # H3.1 2.2 (25)

        return Ui, UL, UG, fL, fG

    def schichtenstroemung(self, fL, fG):
        k = 226.3 ** 2 / np.pi ** 3 * fL * fG ** 2  # H3.1 2.2 (29)
        return k

    def wellenstroemung(self, fG, hl, WeFrL):
        k = 16 * fG ** 3 / (np.pi ** 2 * (1 - (2 * hl - 1) ** 2) ** 0.5) \
            * (np.pi ** 2 / (25 * hl ** 2) * WeFrL ** (-1) + 1)  # H3.1 2.2 (30)
        # print(f"fG: {fG}, hl: {hl}, WeFrL: {WeFrL}")
        # print(f"Wellenströmung k : {k}")
        return k

    def blasenstroemung(self, fG, fL, Ui):
        k = 128 * fG * fL ** 2 / (np.pi ** 2 * Ui)  # H3.1 2.2 (31)
        return k

    def nebelstroemung(self, fL, fG, WeFrL):
        xi_ph = (1.138 + 2 * np.log10(np.pi / (1.5 * fL))) ** (-2)  # H3.1 2.2 (35)
        k = 7680 * fG ** 2 / (np.pi ** 2 * xi_ph) * WeFrL ** (-1)  # H3.1 2.2 (34)
        return k

    def flow_character(self, m_dot, d):
        [X, ReFr05, Fr05, FrEu05, WeFrL] = self.flow_values(m_dot, d)
        # print(f"X: {X}")
        # get flow condition
        global eps
        eps = self.epsilon(m_dot)

        # print(f"eps: {eps}")

        def find_phi(phi, eps):
            rest = 2 * np.pi * eps + np.sin(phi) - phi  # H3.1 2.2 (27)    removed factor in sinus
            return rest

        sol = scipy.optimize.root(find_phi, np.pi, args=(eps))
        phi = sol.x

        # print(f"phi: {phi}")

        def find_hl(hl, self, m_dot, d, X):
            if hl > 1:
                print('whole tube filled with liquid')
                hl = 1 - 1e-5
            if hl < 0:
                print('invalid relative liquid height')
                hl = 1e-5
            # else:
            Ui = 2 * (hl * (1 - hl)) ** 0.5

            if hl <= 0.5:
                psi = 2 * np.arcsin(Ui)
                UL = psi / 2
                UG = np.pi - UL
                fL = (psi - np.sin(psi)) / 8  # removed factor in sinus
                fG = np.pi / 4 - fL
            if hl > 0.5:
                phi = 2 * np.arcsin(Ui)
                UG = phi / 2
                UL = np.pi - UG
                fG = (phi - np.sin(phi)) / 8  # removed factor in sinus
                fL = np.pi / 4 - fG
            # for horizontal tubes FrEu = infinite
            # xi = 0.3164 / (self.Reynolds_v(m_dot,d))**0.25          # H3.1 2.2 (10)
            # FrEu_G = xi * (m_dot*self.q)**2 / ( 2 * d * g * self.rho_v * (self.rho_l - self.rho_v))  # H3.1 2.2 (9)
            FrEu_G = 1e10
            rest = (((UG + Ui) / np.pi) ** 0.25 * np.pi ** 2 / (64 * fG ** 2) \
                    * ((UG + Ui) / fG + Ui / fL) - 1 / FrEu_G) * (np.pi / UL) ** 0.25 \
                   * (64 * fL ** 3 / (np.pi ** 2 * UL)) - X ** 2  # H3.1 2.2 (8)
            # print(f"hl: {hl}, rest: {rest}")
            return rest

        hl0 = self.hlstart(eps, phi)
        # print(f"hl0: {hl0}")

        sol2 = scipy.optimize.root(find_hl, hl0, args=(self, m_dot, d, X), method='lm')
        # print(f"hl: {sol2.x}")
        if sol2.success != True:
            print("ERROR! Solver hl not converged")
        hl = sol2.x

        Ui, UL, UG, fL, fG = self.hl_values(hl)

        # Überprüfung des Strömungszustandes
        if ReFr05 ** 2 <= self.schichtenstroemung(fL, fG):
            zustand = 'Schichtenströmung'
        elif Fr05 ** 2 <= self.wellenstroemung(fG, hl, WeFrL):
            zustand = 'Wellenströmung'
        elif FrEu05 ** 2 >= self.blasenstroemung(fG, fL, Ui):
            zustand = 'Blasenströmung'
        elif X >= 0.36 and self.reynolds_l(m_dot, d) >= 1187 and self.reynolds_v(m_dot, d) >= 1187:
            if Fr05 ** 2 > self.wellenstroemung(fG, hl, WeFrL):
                zustand = 'Pfropfenströmung'
        elif X >= 0.51 and self.reynolds_l(m_dot, d) <= 1187 and self.reynolds_v(m_dot, d) >= 1187:
            if Fr05 ** 2 > self.wellenstroemung(fG, hl, WeFrL):
                zustand = 'Pfropfenströmung'
        elif X < 0.51 and Fr05 ** 2 >= self.nebelstroemung(fL, fG, WeFrL):
            zustand = 'Nebelströmung'
        elif X < 0.51 and Fr05 ** 2 < self.nebelstroemung(fL, fG, WeFrL) and Fr05 ** 2 > self.wellenstroemung(fG, hl,
                                                                                                              WeFrL):
            zustand = 'Ringströmung'
        else:
            print('Kein Zustand gefunden!')
            print(f"Zustand mit {self.q} und fL {fL} und hl {hl}")
        return zustand, phi, eps

    def pressure_loss(self, m_dot, d, z_i):
        # acceleration pressure loss
        deltap_b = m_dot ** 2 * (1 / self.rho_v - 1 / self.rho_l)

        # fraction pressure loss
        # dispersed or continuous?

        a = self.q * self.rho_l / ((1 - self.q) * self.rho_v)
        Fr = self.froude(m_dot, d)
        if a <= 12 * Fr ** 0.5 / (1 + Fr ** 0.5 / 7):
            # dispersed
            Re_zp = m_dot * d / (self.dyn_vis * (1 - self.q(1 - self.dyn_vis_v / self.dyn_vis_l)))
            if (1 / a) <= 0.4:
                K2 = 1 + 0.09 * (1 / a)
            elif (1 / a) > 0.4:
                K2 = (1 - (2.97 / (1 / a) ** (2 / 3) + 1) / (6 * (1.83 / \
                                                                  (1 / a) ** (2 / 3) + 1) * (
                                                                     3.43 / (1 / a) ** (2 / 3) + 1))) ** (-1)

            def find_xi(xi, d, Re_zp):
                rest = -2 * np.log10(Rz / d / 3.7) + 2.51 / (Re_zp * xi ** 0.5) - 1 / xi ** 0.5
                return rest

            xi = scipy.optimize.root(find_xi, 0.1, args=(Re_zp))
            dp_dl = xi * m_dot ** 2 / (2 * self.rho_v * d) * (1 + self.q * (self.rho_l / self.rho_v - 1)) * \
                    (1 - self.q * (self.rho_l / self.rho_v - 1) * (K2 - 1))

        if a > 12 * Fr ** 0.5 / (1 + Fr ** 0.5 / 7):

            # continuous
            Re_G = self.Reynolds(m_dot, d)
            ReL = self.reynolds_l(m_dot, d)
            FrL = self.froude_l(m_dot, d)

            if Re_G > 2300:

                def find_xi2(xi, Re_G):
                    rest = 2 * np.log10(Re_G * xi ** 0.5) - 0.8 - xi ** (-0.5)
                    return rest

                xi = scipy.optimize.root(find_xi2, 0.1, args=(Re_G))
                E = 1.857 + 0.815 * np.log10((m_dot * self.q / (self.rho_v * self.c_v)) ** 2 \
                                             * (1 + 4575 * self.rho_v ** 2 / self.rho_l ** 2))
                psi = (1 - self.q) / self.q * (ReL * FrL) ** (-1 / 6) \
                      * (self.rho_l / self.rho_v) ** (-0.9) * (self.dyn_vis_l / self.dyn_vis_v) ** (-0.5)
                eps2 = psi * 9.1

                if (Rz / d) < 5e-4:
                    # smooth tube
                    eps1 = 1.71 * psi ** 0.2 * ((1 - self.q) / self.q) ** 0.15 \
                           * (self.rho_v / self.rho_l) ** 0.5 * (self.dyn_vis_v / self.dyn_vis_l) ** 0.1
                elif (Rz / d) >= 5e-4:
                    eps1 = 1.71 * psi ** 0.2 * ((1 - self.q) / self.q) ** 0.15 \
                           * (self.rho_v / self.rho_l) ** 0.5 * (self.dyn_vis_v / self.dyn_vis_l) ** 0.1 * \
                           (5e4 / (Rz / d)) ** 0.1

                eps = (eps1 ** (-3) + eps2 ** (-3)) ** (-1 / 3)
                gammaF = 1 - (1 + (1 - self.q) * self.rho_v / (self.q * eps * self.rho_l)) ** (-1.19)
                gammaE = (1 + 6.67 / ((1 - self.q) / self.q) ** 0.45 \
                          * (1 + 3 * self.q ** 4) * (self.dyn_vis_l / self.dyn_vis_v - 1) ** 0.25) ** (-1)

                Phi = (1 / (1 - (1 - E) * gammaF - E * gammaE)) ** 2

                dp_dl = xi.x * m_dot ** 2 / (2 * self.rho_v * d) * Phi

        deltap_r = dp_dl * z_i
        deltap = deltap_r + deltap_b
        return deltap

    def length(self, m_dot, d, q_dot, deltax):
        """
        calculates length for average q_dot

        Parameters
        ----------
        m_dot : TYPE float.
            DESCRIPTION. mass flow per area [kg/s/m**2]
        d : TYPE float.
            DESCRIPTION. diameter [m]
        q_dot : TYPE float.
            DESCRIPTION. average heat flow density [kg/m**2/s]
        deltax : TYPE float.
            DESCRIPTION. total difference in vapour quality []

        Returns
        -------
        l : TYPE float.
            DESCRIPTION. length [m]

        """
        print(f"deltax: {deltax}, hv: {self.h_evap}, q_dot: {q_dot}")
        l = m_dot * np.pi * 0.25 * d ** 2 * deltax * self.h_evap / (q_dot * np.pi * d)
        return l


if __name__ == '__main__':
    d_i = 12e-3
    T = 372
    m_dot = 0.01 / (0.25 * d_i ** 2 * np.pi)
    fluid_name = 'Butane'
    x_var = np.linspace(0.1, 1 - 1e-5, 10)
    alpha_complete = []
    for vapourquality in x_var:
        p1 = PointND(T, vapourquality, fluid_name, 'Testpunkt Isobutan')
        alpha = p1.alpha_sieden(m_dot, d_i, 0.1)
        # print(p1.flow_character(300, 0.03))
        # print("alpha: ",alpha[1])
        alpha_complete.append(alpha)
        print("Strömungsform:", p1.flow_character(m_dot, d_i)[0], "\n")

    print(f"alpha_complete: {alpha_complete}")
    plt.plot(x_var, alpha_complete[0], 'r.', label="Gesamt")
    plt.plot(x_var, alpha_complete[1], 'g.', label="blasen")
    plt.plot(x_var, alpha_complete[2], 'b.', label="konvektiv")
    plt.legend()
    plt.show()

    # delta_p =p1.druckverlust(2.1,0.015,0.1)

    # plt.figure(2)
    # plt.plot(np.linspace(0.1,1-1e-5,10),h)
    # plt.title("Höhe Flüssigkeitsspiegel")

    if (p1.inlet_flow != 'true') and (p1.inlet_flow != 'false'):
        print('Unbekannter hydrodynamischer Status. Check Einlauf.')
        sys.exit()
    if p1.RB == 'qzu':
        if p1.q_dot == None or p1.q_dot == 0:
            print("Bitte Wärmestromdichte vorgeben")
            sys.exit

    if (p1.RB != 'Tw') and (p1.RB != 'qzu'):
        print('Unbekannter hydrodynamischer Status. Check Einlauf.')
        sys.exit()
