# -*- coding: utf-8 -*-
"""
Druckabfall und Wärmeübergang im Doppelrohr-Waermeuebertrager im
Zweiphasengebiet.
Druckabfall zunächst nur im Innenrohr.
Sekundärfluid strömt im Gegenstrom und ist bisher einphasig!
Löser nun: Randwertlöser
Ueberarbeitung 19.12.2020: low-level-Interface aus Coolprop
          mit tabellierten Daten wird genutzt

Benutzung der Bibliotheken ht und fluid für DRuckabfall und alpha
Benutzung von Coolprop für Stoffeigenschaften

Created on Wed Dec  9 17:37:20 2020

@author: atakan
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.optimize import differential_evolution
import CoolProp.CoolProp as CP
from fluids.two_phase import two_phase_dP
import fluids.friction as frict
import fluids.core as flc
from ht.condensation import Shah as alpha_shah
from ht.conv_internal import Nu_conv_internal as nu_conv_in
# frfactV = np.vectorize(frict.friction_factor)
# nu_conv_inV = np. vectorize(nu_conv_in)
# one_Phase_dp_V = np.vectorize(frict.one_phase_dP)

__vielPrint__ = False

m_dot = 0.0085  # Massenstrom kg/s
x_0 = 1.
p_0 = 20e5     # Anfangsdruck Pa
dtemp_ueberhitz = 20
diameter = 6.8e-3  # Innenrohrdurchmesser
l_max = 10.0  # Rohrlaenge
fluid_a = "n-Propane"  # Working fluid
temp_sur = 283
p_sur = 1.0135e5
area = np.pi * (diameter / 2)**2
m_dot_area = m_dot / area  # Massenstromdichte
p_c = CP.PropsSI('Pcrit', fluid_a)
temp_0 = CP.PropsSI('T', 'P', p_0, 'Q', x_0, fluid_a) + \
          dtemp_ueberhitz  # initial temeperature/state

h_0 = CP.PropsSI('H', 'P', p_0, 'T', temp_0, fluid_a)
working_fluid = CP.AbstractState("BICUBIC&HEOS", fluid_a)
# Sekundärfluid --------------------------------
fluid_s = "Water"
secondary_fluid = CP.AbstractState("BICUBIC&HEOS", fluid_s)
p_s = 4e5
m_dot_s = m_dot *1.9 # 0.05125
diameter_s = 16e-3
area_s = np.pi * (diameter_s / 2)**2 - area
m_dot_area_s = m_dot_s / area_s  # Massenstromdichte
m_dot_area = m_dot / area
temp_0_s = temp_sur
h_0_s = CP.PropsSI('H', 'P', p_s, 'T', temp_0_s, fluid_s)
h_end = CP.PropsSI('H', 'P', p_0, 'T', temp_0_s, fluid_a)


n_space = 40  # initial discretization in space/no. of points

if __vielPrint__:
    print("Fluid %10s, T: %3.1f K, p: %6.1f, h: %5.2f" %
          (fluid_a, temp_0, p_0, h_0))


def mdot_area_function(m_dot, diameter):
    area = np.pi * (diameter / 2)**2
    m_dot_area = m_dot / area
    return m_dot_area, area


def properties_sat(p, fluid):
    """
    Properties needed for integration at given p and h.

    Parameters
    ----------
    p : float
        pressure in Pa.

    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        includes: tranport properies in saturated state at given pressure p
        densities of liquid and vapor,
        viscosities of liquid and vapor, cp of liquid, conductivity of liquid

        all in SI units.

    """

    # liquid:
    fluid.update(CP.PQ_INPUTS, p, 0)
    reihe = [CP.iDmass, CP.iCpmass, CP.iviscosity, CP.iconductivity]
    props = [fluid.keyed_output(k) for k in reihe]
    rhol_0, cpl, mul_0, lambda_l = props[:]

    # vapor:
    fluid.update(CP.PQ_INPUTS, p, 1)
    reihe = [CP.iDmass, CP.iviscosity]
    props = [fluid.keyed_output(k) for k in reihe]
    rhov_0,  muv_0, = props[:]

    alle = [rhol_0, rhov_0, mul_0, muv_0, cpl, lambda_l]

    return alle


def properties_satV(p, h, fluid): # unbenutzte vektorisierung
    _n = len(p)
    alle = np.zeros((14, _n))
    for _i in range(_n):
        alle[:, _i] = properties_sat(p[_i], h[_i], fluid)
    return alle


def properties(p, h, fluid):
    """
    Properties needed for integration at given p and h, single phase.

    Parameters
    ----------
    p : float
        pressure in Pa.
    h : float
        specific enthalpy in J/kg.
    fluid :   an AbstractState in coolprop.

    Returns
    -------
    alle : numpy array
        includes: T, p , quality, specific enthalpy,entropy
        densities,
        viscosities, cp , conductivity,phase prandtl-number
        at the defined state (p,h)
        all in SI units.

    """
    fluid.update(CP.HmassP_INPUTS, h, p)
    reihe = [CP.iT, CP.iQ, CP.iSmass, CP.iDmass, CP.iPrandtl, CP.iPhase,
             CP.iconductivity, CP.iCpmass, CP.iviscosity]
    props = [fluid.keyed_output(k) for k in reihe]
    _temp, x, s, rho, prandtl, phase, lambda_s, cp, mu = props[:]
    alle = [_temp, p, x, h,  s, rho, mu,
            cp, lambda_s, phase, prandtl]

    return alle


def properties_V(p, h, fluid):
    _n = len(p)
    alle = np.zeros((11, _n))
    for _i in range(_n):
        alle[:, _i] = properties(p[_i], h[_i], fluid)
    return alle

def bilanzen(ort, ph, diameters, m_dots, fluids, p_c_w):
    """
    ODE balances for pressure drop and enthalpy change in tube flow.

    The other fluid is at constant temperature temp_sec

    Parameters
    ----------
    ort : float
        position (in m).
    ph : np.array of 4 floats
        ph[0] = presure in Pa
        ph[1] =specific enthalpy in J/kg.
        ph[2] = presure in Pa of secondary/outer fluid
        ph[3] =specific enthalpy in J/kg. of secondary/outer fluid
    fluids : list of 2 strings
        two CoolProp fluids [0]: inner fluid,[1] secondary/outer fluid.
    m_dots : array of 2 floats
        mass flow rate in kg/s.
    diameters: array of two floats
            inner/outer dimeter[0],[1].
    p_c_w: critical pressure of working fluid


    Returns
    -------
    diffs : numpy-array of 2 floats
        differntial change in pressure and enthalpy
        dp/dort in Pa/m
        dh/dort in kJ/kg/m.

    """

    # print("ph:", ort, m_dot, ph, diameter)
    temp, p, _x, h,  s, rho, mu, cp, lambda_s, phase, prandtl =\
        properties(ph[0], ph[1], fluids[0])
    temp_s, p_s, x_s, h_s, s_s, rho_s, mu_s, cp_s, lambda_s,  phase_s, \
        prandtl_s = properties(ph[2], ph[3], fluids[1])
    m_dot_area, area = mdot_area_function(m_dots[0], diameters[0])
    area_s = np.pi * (diameters[1] / 2)**2 - area
    m_dot_area_s = m_dots[1] / area_s
    d_h = diameters[1] - diameters[0]
    rhol, rhov, mul, muv, cpl, lambda_l =\
        properties_sat(ph[0], fluids[0])

    if _x >= 0:  # 2 phase flow
        alpha = alpha_shah(m_dots[0], _x, diameters[0], rhol,
                           mul, lambda_l, cpl, p, p_c_w)
        dp = -two_phase_dP(m_dots[0], _x, rhol, diameters[0], rhog=rhov,
                           mul=mul, mug=muv, P=p, Pc=p_c_w)
        if __vielPrint__:
            print("alpha: %5.2f W /(m2 K), Ort: %2.4f m" % (alpha, ort))
    elif _x < 0:  # must be changed to single phase corelations
        rey = flc.Reynolds(m_dot_area/rho, diameters[0], rho, mu=mul)
        friction_fac = frict.friction_factor(rey)

        if ort == 0:
            nusselt = nu_conv_in(
                rey, prandtl, Di=diameters[0],  fd=friction_fac)
        else:
            nusselt = nu_conv_in(rey, prandtl, Di=diameters[0],
                                 x=ort, fd=friction_fac)

        alpha = nusselt * lambda_l / diameters[0]
        dp = -frict.one_phase_dP(m_dots[0], rho, mul, diameters[0])
        # print("alpha: %5.2f W /(m2 K), Ort: %2.4f m ### x=0, Re: %5.1f"
        #         % (alpha, ort, rey))
    else:
        print("else:", _x)
        alpha = 0
        dp = 0

    # Secondary fluid ------------------------------
    rey_s = flc.Reynolds(m_dot_area_s/rho_s, d_h, rho_s, mu=mu_s)
    friction_fac_s = frict.friction_factor(rey_s)

    if ort == 0:
        nusselt_s = nu_conv_in(
            rey_s, prandtl_s, Di=d_h,  fd=friction_fac_s)
    else:
        nusselt_s = nu_conv_in(rey_s, prandtl_s, Di=d_h,
                               x=ort, fd=friction_fac_s)
    alpha_s = nusselt_s * lambda_s / d_h * 0.86 * (diameters[0] /
                                                   diameters[1])**(-.16)
    alpha_t = 1 / (1 / alpha_s + 1 / alpha)

    dh = alpha_t * diameters[0] * np.pi / m_dots[0] * (temp_s - temp)
    dh_s = alpha_t * diameters[0] * np.pi / m_dots[1] * (temp_s - temp)
    # all_alpha.append(np.array([alpha_t, alpha, alpha_s]))

    # if True :
    # print("dh: %4.4f , T:%3.3f K" % (dh, temp))
    #print("-# %5f, %5f, %5f, %5f , %5f" %( rey_s, temp, alpha, alpha_s, alpha_t))
    diffs = np.array([dp, dh, 1e-7, dh_s])
    return diffs


def bilanzenV(ort, ph, diameters, m_dots, fluids, p_c_w):  # vectorization ...
    _n = len(ort)
    diffs = np.zeros((4, _n))
    for _i in range(_n):
        diffs[:, _i] = bilanzen(ort[_i], ph[:, _i], diameters,
                                m_dots, fluids, p_c_w)
    return diffs


def bcs(right, left):
    global p_0, h_0, p_s, h_0_s
    return np.array((right[0]-p_0, right[1]-h_0, left[2]-p_s, left[3]-h_0_s))


def h_exch_entropy(op, m_dot, fluids, p_c_w=p_c):
    m_dot_s, diameter, diameter_s = op
    diameters = [diameter, diameter_s]
    m_dots = [m_dot, m_dot_s]
    x_all = np.linspace(0, l_max, n_space)
    anfwert = np.zeros((4, n_space))
    para = np.polyfit([0, l_max], [h_0, h_end], 1)

    anfwert[0, :] = p_0
    anfwert[2, :] = p_s
    anfwert[3, :] = h_0_s
    #anfwert += 1
    anfwert[1, :] = np.polyval(para, x_all)
    print("A: program running, please wait ...")
    ph_verlauf = solve_bvp(fun=
            lambda x, y: bilanzenV(x, y, diameters, m_dots, fluids, p_c_w),
            bc=bcs, x=x_all, y=anfwert,  max_nodes=500, tol=0.0051, verbose=0)
    # ph_verlauf = solve_bvp(bilanzen, bcs, x_all, anfwert)
    n_val = len(ph_verlauf.x)
    alle = properties(ph_verlauf.y[0, 0], ph_verlauf.y[1, 0], working_fluid)
    alle_s = properties(ph_verlauf.y[2, 0], ph_verlauf.y[3, 0], secondary_fluid)
    h_w_sur = CP.PropsSI('H', 'P', p_sur, 'T', temp_sur, fluid_a)
    s_w_sur = CP.PropsSI('S', 'P', p_sur, 'T', temp_sur, fluid_a)
    exergy_start = m_dot * ((alle[3] - h_w_sur) -
                        temp_sur * (alle[4] - s_w_sur))
    h_s_sur = CP.PropsSI('H', 'P', p_sur, 'T', temp_sur, fluid_s)
    s_s_sur = CP.PropsSI('S', 'P', p_sur, 'T', temp_sur, fluid_s)
    exergy_end = m_dot_s * ((alle_s[3] - h_s_sur) -
                        temp_sur * (alle_s[4] - s_s_sur))
    print(exergy_start, exergy_end,exergy_start - exergy_end,alle[0], alle_s[0])
    #print( "%3.3f %3.3f %3.3f " % (op[0],op[1],op[2]))
    return  - exergy_end/m_dot_s

fluids = [working_fluid, secondary_fluid]
bc_opt = ((.004, .02), (3.5e-3, 14e-3), (13e-3, 30e-3))
nm = 4
nid = 4
res_eval=np.zeros((nm,nid))
msec = np.linspace(bc_opt[0][0], bc_opt[0][1], nm)
ds =  np.linspace(bc_opt[1][0], bc_opt[1][1], nid)
plt.figure()
for ii, mds in enumerate(msec):
    for jj, idia in enumerate(ds):
        res_eval[ii,jj] = h_exch_entropy([mds, idia, 20e-3], m_dot, fluids)
        print(ii, jj)
    plt.plot(ds, -res_eval[ii, :], label="%2.2f g/s" % (mds * 1000))
plt.legend()
plt.xlabel("d / m")
plt.ylabel("d_sp_exergy_sec. fl. / (kJ/kg)")
plt.title("L: %2i m,  %8s_%8s_m%3i g/s" %
          (l_max, fluid_a, fluid_s, m_dot*1000))
plt.savefig("sp_Length_%2i_F%8s_%8s_m%3i.png" %
            (l_max, fluid_a, fluid_s, m_dot*1000))

try:
    optimal = differential_evolution(h_exch_entropy, bc_opt,
                                     args=(m_dot, fluids, p_c))
    print(optimal)
except:
    print("Optimierung mit Fehler")

m_dot_s, diameter, diameter_s = [0.08999934, 0.00454419, 0.02790231]
x_all = np.linspace(0, l_max, n_space)
anfwert = np.zeros((4, n_space))
para = np.polyfit([0, l_max], [h_0, h_end], 1)

anfwert[0, :] = p_0
anfwert[2, :] = p_s
anfwert[3, :] = h_0_s
#anfwert += 1
anfwert[1, :] = np.polyval(para, x_all)
print("program running, please wait ...")
diameters =[diameter, diameter_s]
m_dots = [m_dot, m_dot_s]

ph_verlauf = solve_bvp(fun=
            lambda x,y: bilanzenV(x,y, diameters, m_dots, fluids, p_c),
            bc=bcs, x=x_all, y=anfwert,  max_nodes=500, tol=0.0051, verbose=0)
# ph_verlauf = solve_bvp(bilanzen, bcs, x_all, anfwert)
n_val = len(ph_verlauf.x)


# zum Auswerten und Plotten: -----------------------------------------
alle = properties_V(ph_verlauf.y[0, :], ph_verlauf.y[1, :], working_fluid)
alle_s = properties_V(ph_verlauf.y[2, :], ph_verlauf.y[3, :], secondary_fluid)
q_dot = m_dot * (alle[3][-1] - alle[3][0])
# Evaluating all changes:
entropy_change_w = m_dot * (alle[4][-1] - alle[4][0])
entropy_change_s = m_dot_s * (alle_s[4][0] - alle_s[4][-1])
entropy_production = entropy_change_w + entropy_change_s
dex_w = m_dot * ((alle[3][-1] - alle[3][0]) -
                 temp_sur * (alle[4][-1] - alle[4][0]))
dex_s = -m_dot_s * ((alle_s[3][-1] - alle_s[3][0]) -
                    temp_sur * (alle_s[4][-1] - alle_s[4][0]))
exergetic_efficiency = dex_s / dex_w
ex_loss = entropy_production * temp_sur

print("Q_dot: %3.3f kW, m_dot/Area: %3.3f kg/(s m2)" % (q_dot/1e3, m_dot_area))
print("Ex_eff: %3.3f , entropy-prod-rate: %3.3f W/K" %
      (exergetic_efficiency, entropy_production))
print("Ex_loss: %3.3f  W" % (ex_loss))
#  Evaluating only the useful exergy of the secondary fluid:
# Pressure levels are neglected/set to initial values for exergy calculations here
h_w_sur = CP.PropsSI('H', 'P', p_0, 'T', temp_sur, fluid_a)
s_w_sur = CP.PropsSI('S', 'P', p_0, 'T', temp_sur, fluid_a)
exergy_start = m_dot * ((alle[3][0] - h_w_sur) -
                        temp_sur * (alle[4][0] - s_w_sur))
h_s_sur = CP.PropsSI('H', 'P', p_s, 'T', temp_sur, fluid_s)
s_s_sur = CP.PropsSI('S', 'P', p_s, 'T', temp_sur, fluid_s)
exergy_end = m_dot_s * ((alle_s[3][0] - h_s_sur) -
                        temp_sur * (alle_s[4][0] - s_s_sur))

print("Ex_start: %3.3f  W, Ex_end_sec: %3.3f  W, Ex_eff: %3.3f " %
      (exergy_start, exergy_end, exergy_end/exergy_start))
h_w_sur = CP.PropsSI('H', 'P', p_sur, 'T', temp_sur, fluid_a)
s_w_sur = CP.PropsSI('S', 'P', p_sur, 'T', temp_sur, fluid_a)
exergy_start = m_dot * ((alle[3][0] - h_w_sur) -
                        temp_sur * (alle[4][0] - s_w_sur))
h_s_sur = CP.PropsSI('H', 'P', p_sur, 'T', temp_sur, fluid_s)
s_s_sur = CP.PropsSI('S', 'P', p_sur, 'T', temp_sur, fluid_s)
exergy_end = m_dot_s * ((alle_s[3][0] - h_s_sur) -
                        temp_sur * (alle_s[4][0] - s_s_sur))

print("\nEx_start: %3.3f  W, Ex_end_sec: %3.3f  W, Ex_eff: %3.3f relative to p_sur " %
      (exergy_start, exergy_end, exergy_end/exergy_start))
# plotting ------------------------------------------
y_achse = ["T / K", "p / (Pa)", "x / [-]", "h / (J / kg)"]
zahl_plots = 4

fig, ax1 = plt.subplots(zahl_plots, 1, sharex=True)

for ii in range(4):
    plt.subplot(zahl_plots, 1, ii+1)
    plt.plot(ph_verlauf.x, alle[ii])
    plt.plot(ph_verlauf.x, alle_s[ii], "r.")
    if ii == 3:
        plt.xlabel("x / m")
    plt.ylabel(y_achse[ii])

print("p(Ende):  %6.3f Pa, T_ende %3.2f K, T_sec_out %3.2f K" %
      (alle[1][-1], alle[0][-1], alle_s[0][0]))

# Vergleich kinetische Energie mit Totalenthalpie ----------------
ekin_plot = False
if ekin_plot:
    plt.figure()
    plt.plot(ph_verlauf.x, 0.5 * (m_dot_area / (alle[-3]))**2, "k-")  # E_kin
    plt.plot(ph_verlauf.x, alle[3]-(alle[3])[0], "b.")  # Totalenthalpie
    plt.xlabel("x / m")
    plt.ylabel(" (E/m) / (m2 / s2)")
