# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 20:40:12 2023

@author: Max
"""

"==========================================   IMPORT   ==========================================="

from SALib.sample import saltelli

import numpy as np
import datetime
import os

from scipy.integrate import solve_ivp
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from ht.conv_internal import Nu_conv_internal as nu_corr
import fluids.friction as frict
import fluids.core as flc

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum
SI = RP.GETENUMdll(0, "SI").iEnum


def diffeq_enthalpy_ivp(ort, y, input_values, wf='Isobutane', sf='Water'):
    '''
    Differential equations of enthalpy for working and secondary fluid in an double-pipe heat exchanger.
    '''
    p_wf_0, p_sf_0, di, ds, dh, da, area_wf, area_sf, lam_rohr, m_wf, m_sf, z, T_sf_0 = input_values
    h_wf, h_sf = y

    # Stoffeigenschaften der Fluide: Dichte, Wärmekapazität, kinematische Viskosität, Thermische Leitfähigkeit, Prandtl-Zahl, Temperatur
    Zustand_wf_x = RP.REFPROPdll(wf, "PH", "D;CP;VIS;TCX;PRANDTL;T", MASS_BASE_SI, 0, 0, 10e5, h_wf, [0]).Output[0:6]
    Zustand_sf_x = RP.REFPROPdll(sf, "PH", "D;CP;VIS;TCX;PRANDTL;T", MASS_BASE_SI, 0, 0, 1e5, h_sf, [0]).Output[0:6]

    # Temperaturen der Fluide in: K
    T_wf = Zustand_wf_x[5]
    T_sf = Zustand_sf_x[5]

    # Wärmeübergangskoeffizient
    alpha_wf = 1000
    alpha_sf = 1000

    # Übertragener Wärmestrom
    R = (alpha_wf * np.pi * di) ** -1 + np.log((di + 2 * ds) / di) / (2 * np.pi * lam_rohr) + (
                alpha_sf * np.pi * (di + 2 * ds)) ** -1
    dQdx_R = (T_wf - T_sf) / R

    # DGL
    dhwf = 1 / (m_wf) * dQdx_R
    dhsf = 1 / (m_sf) * dQdx_R

    return [dhwf, dhsf]


def run_model(parameter_values, sf="Water", wf="Isobutane"):
    global L
    lam_rohr = 50
    ds = 1e-3
    z = [0]

    p_wf_0, p_sf_0, dT_max, m_wf, m_sf, c_wf, c_sf = parameter_values

    T_wf_sat = RP.REFPROPdll(wf, "PQ", "T", MASS_BASE_SI, 0, 0, p_wf_0, 1, z).Output[0]
    T_sf_0 = T_wf_sat - dT_max

    roh_wf = RP.REFPROPdll(wf, "PQ", "D", MASS_BASE_SI, 0, 0, p_wf_0, 1, z).Output[0]
    roh_sf = RP.REFPROPdll(sf, "PT", "D", MASS_BASE_SI, 0, 0, p_wf_0, T_sf_0, [0]).Output[0]

    area_wf = m_wf / (c_wf * roh_wf)
    area_sf = m_sf / (c_sf * roh_sf)

    di = np.sqrt(4 * area_wf / np.pi)
    da = di + 2 * ds
    dA = np.sqrt(4 * area_sf / np.pi + da ** 2)
    dh = dA - da

    input_values = [p_wf_0, p_sf_0, di, ds, dh, da, area_wf, area_sf, lam_rohr, m_wf, m_sf, z, T_sf_0]

    h_wf_sat = RP.REFPROPdll(wf, "PQ", "H", MASS_BASE_SI, 0, 0, p_wf_0, 1, z).Output[0]
    h_sf_0 = RP.REFPROPdll(sf, "PT", "H", MASS_BASE_SI, 0, 0, p_sf_0, T_sf_0, [0]).Output[0]
    y_bc = [h_wf_sat, h_sf_0]

    orte = np.linspace(0, L, 100)  # Definition der Stützstellen
    int_method = 'BDF'

    res = solve_ivp(lambda ort, y: diffeq_enthalpy_ivp(ort, y, input_values, wf=wf, sf=sf), (0, L),
                    y_bc, t_eval=orte, method=int_method)

    # res = solve_ivp(lambda ort, y: diffeq_enthalpy_ivp(ort,y, input_values, wf=wf, sf=sf), (0, L),
    #                  y_bc, t_eval = orte, method=int_method, events = (event_saturatedstate, event_dT_max))

    return res


def event_saturatedstate(x, h):
    global sf, deltaT
    T_sat = RP.REFPROPdll(sf, "PQ", "T", MASS_BASE_SI, 0, 0, 1e5, 0, [0]).Output[0]
    T_abbruch = T_sat - deltaT
    h_abbruch = RP.REFPROPdll(sf, "PT", "H", MASS_BASE_SI, 0, 0, 1e5, T_abbruch, [0]).Output[0]
    # print('event_saturatedstate:' + str(h[1]-h_abbruch))
    return h[1] - h_abbruch


event_saturatedstate.terminal = True


def event_dT_max(x, h):
    global sf, deltaT_max
    T_wf = RP.REFPROPdll(wf, "PH", "T", MASS_BASE_SI, 0, 0, 15e5, h[0], [0]).Output[0]
    T_sf = RP.REFPROPdll(sf, "PH", "T", MASS_BASE_SI, 0, 0, 1e5, h[1], [0]).Output[0]
    dT = T_wf - T_sf
    # print('event_dT_max:' + str(dT-deltaT_max))
    return dT - deltaT_max


event_dT_max.terminal = True

if __name__ == "__main__":
    sf = "Water"
    wf = "Isobutane"

    int_method = 'BDF'
    L = 2

    deltaT_max = 50  # Eventbedingung, größtest lokales Delta T
    deltaT = 10  # Eventbedingung, delta T zur Sattdampftemperatur des Wassers

    p_wf_0 = 15e5
    p_sf_0 = 1e5
    dT_max = 10
    m_wf = 10e-3
    m_sf = 40e-3
    c_wf = 1
    c_sf = 1

    parameter_values = p_wf_0, p_sf_0, dT_max, m_wf, m_sf, c_wf, c_sf

    print("\n______________________________________________\nrun sensitivity analyis....\n")

    for i in range(10000):
        solution = run_model(parameter_values)
        print('Run: ' + str(i) + ', Prozent: ' + str(round(i / 10000 * 100, 3)) + ' %')




