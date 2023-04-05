'''
calculate benchmarks for lab setup
author: Alexandra Welp
1: compressor inlet
2: compressor outlet
3: condenser outlet
4: evaporator inlet
'''

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import ht
import fluid_properties_rp as rp

fluid = "REFPROP::IsoButane[0.5]&Pentane[0.5]"

p_o = 103000        # lowest possible pressure
#T_o = CP.PropsSI("T", "P", p_o, "Q", 1, fluid)     # zeotropic fluid mixtures temperature for x = 1
T_4c = CP.PropsSI("T", "P", p_o, "Q", 1, fluid)
T_4b = CP.PropsSI("T", "P", p_o, "Q", 0, fluid)

p_c = 13e5
#T_c = CP.PropsSI("T", "P", p_c, "Q", 1, fluid)
T_2c = CP.PropsSI("T", "P", p_c, "Q", 0, fluid)
T_2b = CP.PropsSI("T", "P", p_c, "Q", 1, fluid)

delta_T_4c_1 = 20          # delta T superheated low pressure
eta_compressor = 0.65

T_1 = (T_4c + delta_T_4c_1)
h_1 = CP.PropsSI("H", "T", T_1, "P", p_o, fluid)
s_1 = CP.PropsSI("S", "T", T_1, "P", p_o, fluid)
h_2s = CP.PropsSI("H", "S", s_1, "P", p_c, fluid)

h_2 = (h_2s - h_1) / eta_compressor + h_1
T_2 = CP.PropsSI("T", "P", p_c, "H", h_2, fluid)
s_2 = CP.PropsSI("S", "T", T_2, "P", p_c, fluid)

w_compressor = h_2 - h_1
T_2 = CP.PropsSI("T", "H", h_2, "P", p_c, fluid)
delta_T_2_2b = T_2 - T_2b

h_2b = CP.PropsSI("H", "P", p_c, "Q", 1, fluid)
s_2b = CP.PropsSI("S", "P", p_c, "Q", 1, fluid)
h_2c = CP.PropsSI("H", "P", p_c, "Q", 0, fluid)
s_2c = CP.PropsSI("S", "P", p_c, "Q", 0, fluid)

delta_h_sh = h_2 - h_2b
delta_h_ws = h_2b - h_2c


complete_subcooling = 'off'
if complete_subcooling == 'on':
    delta_T_4_4b = 10          # deltaT subcooled low pressure
    T_4 = T_4b - delta_T_4_4b
    h_4 = CP.PropsSI("H", "T", T_4, "P", p_o, fluid)
    h_3 = h_4                   # throttle isenthalp
    T_3 = CP.PropsSI("T", "P", p_c, "H", h_3, fluid)

limited_subcooling = 'on'
if limited_subcooling == 'on':
    delta_Tu_T3 = 20
    Tu = 25 + 273.15
    T_3 = Tu + delta_Tu_T3
    h_3 = CP.PropsSI("H", "T", T_3, "P", p_c, fluid)
    h_4 = h_3
    T_4 = CP.PropsSI("T", "H", h_4, "P", p_o, fluid)
    if T_4 >= T_4b:
        T_4b = T_4

h_4b = CP.PropsSI("H", "P", p_o, "T", T_4b, fluid)
s_4b = CP.PropsSI("S", "P", p_o, "T", T_4b, fluid)
h_4c = CP.PropsSI("H", "P", p_o, "Q", 1, fluid)
s_4c = CP.PropsSI("S", "P", p_o, "Q", 1, fluid)
s_3 = CP.PropsSI("S", "T", T_3, "P", p_c, fluid)
s_4 = CP.PropsSI("S", "P", p_o, "H", h_4, fluid)

delta_h_sc = h_2c - h_3

delta_h_o = h_1 - h_4




print(f"delta_h superheated hp: {delta_h_sh} J/kg")
print(f"delta_h wet steam hp: {delta_h_ws} J/kg")
print(f"delta_h subcooled: {delta_h_sc} J/kg")


# Calculation counter flow heat exchanger
m_var = np.linspace(0.01, 0.025, 10)
m_dot = 0.010
rho_l = CP.PropsSI("D", "Q", 0, "P", p_c, fluid)
rho_v = CP.PropsSI("D", "Q", 1, "P", p_c, fluid)
mu_l = CP.PropsSI("VISCOSITY", "Q", 0, "P", p_c, fluid)
mu_v = CP.PropsSI("VISCOSITY", "Q", 1, "P", p_c, fluid)
k_l = CP.PropsSI("CONDUCTIVITY", "Q", 0, "P", p_c, fluid)
cp_l = CP.PropsSI("C", "Q", 0, "P", p_c, fluid)
cp_v = CP.PropsSI("C", "Q", 1, "P", p_c, fluid)
p_crit = CP.PropsSI("PCRIT", fluid)
d_i = 12e-3
s = 1.5e-3
d_a = 18e-3
d_ai = d_i + 2 * s

def T_log(T1ein, T1aus, T2ein, T2aus):
    deltaTein = T1ein - T2aus
    deltaTaus = T1aus - T2ein
    if deltaTein == deltaTaus:
        delta_T_log = deltaTein
    else:
        delta_T_log = (deltaTein - deltaTaus) / np.log(deltaTein / deltaTaus)
    return abs(delta_T_log)

# secondary fluid

T_pinch = 10
T_sc_in = T_3 - T_pinch
T_sc_out = T_2c - T_pinch
T_sh_in = T_2b - T_pinch
T_sh_out = T_2 - T_pinch
T_ws_out = T_2b - 5
delta_T_ws = 25
T_ws_in = T_ws_out - delta_T_ws
T_evap_out = T_4 + T_pinch
T_evap_in = T_1 + T_pinch
cp_w = 4181.9
cp_methanol = CP.PropsSI("C", "T", 0.5 * (T_evap_in + T_evap_out), "P", 1e5, "REFPROP::methanol")



# size of heat storages
rho_w = 972.42      # density water
timestorage = 3 * 3600
rho_thermooil = 900
cp_thermooil = 2300
Q_storage_sc = m_dot * delta_h_sc * timestorage
Q_storage_ws = m_dot * delta_h_ws * timestorage
Q_storage_sh = m_dot * delta_h_sh * timestorage
Q_storage_cold = m_dot * delta_h_o * timestorage

delta_T_evap = T_evap_in - T_evap_out
delta_T_sc = T_sc_out - T_sc_in

delta_T_sh = T_sh_out - T_sh_in

print(CP.PhaseSI( "T", 281, "D", 801, "Methanol"))
rho_methanol = CP.PropsSI("D" ,"T", T_1, "P", 1e5, "Methanol")

m_methanol = Q_storage_cold / (cp_methanol * delta_T_evap)
m_sc = Q_storage_sc / (cp_w * delta_T_sc)
m_ws = Q_storage_ws / (cp_w * delta_T_ws)
m_sh = Q_storage_sh / (cp_thermooil * delta_T_sh)

V_methanol = m_methanol / rho_methanol
V_sc = m_sc / rho_w
V_ws = m_ws / rho_w
V_sh = m_sh / rho_thermooil

print(f"Q_storage_methanol: {Q_storage_cold/1000} kJ = {Q_storage_cold/3600/1000} kWh, {V_methanol} m^3")
print(f"Q_storage_sc: {Q_storage_sc/1000} kJ = {Q_storage_sc/3600/1000} kWh, {V_sc} m^3")
print(f"Q_storage_ws: {Q_storage_ws/1000} kJ = {Q_storage_ws/3600/1000} kWh, {V_ws} m^3")
print(f"Q_storage_sh: {Q_storage_sh/1000} kJ,= {Q_storage_sh/3600/1000} kWh, {V_sh} m^3")


# calculate heat transfer coefficients

def alpha_1P_i(p, T, fluid, m_dot, d_i):
    vis = CP.PropsSI("VISCOSITY", "P", p, "T", T, fluid)
    rho = CP.PropsSI("D", "P", p, "T", T, fluid)
    Pr = CP.PropsSI("PRANDTL", "P", p, "T", T, fluid)
    lam = CP.PropsSI("L", "P", p, "T", T, fluid)
    Re = 4 * m_dot / (vis * np.pi * d_i)
    Nu = ht.conv_internal.Nu_conv_internal(Re, Pr, Di=d_i)
    alpha = lam * Nu / d_i
    return alpha


def alpha_1P_annulus(p, T, fluid, m_dot, d_ai, d_a):
    if fluid == "thermooil":
        # data from Therminol 66 data sheet for 130 °C
        vis = 2.05e-3
        rho = 935
        Pr = 35.82
        lam = 0.111
    else:
        vis = CP.PropsSI("VISCOSITY", "P", p, "T", T, fluid)
        rho = CP.PropsSI("D", "P", p, "T", T, fluid)
        Pr = CP.PropsSI("PRANDTL", "P", p, "T", T, fluid)
        lam = CP.PropsSI("L", "P", p, "T", T, fluid)
    Re = 4 * m_dot * (d_a - d_ai) / (vis * np.pi * (d_a ** 2 - d_ai ** 2))
    if Re > 2300:
        Nu = 0.023 * Re ** 0.8 * Pr ** 0.4
    else:
        Nu = 3.66 + 1.2 * (d_ai / d_a) ** -0.8
    alpha = lam * Nu / (d_a - d_ai)
    return alpha


m_sc_dot = m_dot * delta_h_sc / (cp_w * delta_T_sc)
m_sh_dot = m_dot * delta_h_sh / (cp_thermooil * delta_T_sh)
m_ws_dot = m_dot * delta_h_ws / (cp_w * delta_T_ws)
m_evap = m_dot * delta_h_o / (cp_methanol * delta_T_evap)

u_sc = m_sc_dot / (rho_w * np.pi * 0.25 * (d_a ** 2 - d_ai ** 2))
u_ws = m_ws_dot / (rho_w * np.pi * 0.25 * (d_a ** 2 - d_ai ** 2))
u_wf = m_dot / (rho_v * np.pi * 0.25 * d_i ** 2)
u_sh = m_sh_dot / (935 * np.pi * 0.25 * (d_a ** 2 - d_ai ** 2))
u_evap = m_evap / (rho_methanol * np.pi * 0.25 * (d_a ** 2 - d_ai ** 2))

alpha_sc_i = alpha_1P_i(p_c, 0.5 * (T_2c + T_3), fluid, m_dot, d_i)
alpha_sh_i = alpha_1P_i(p_c, 0.5 * (T_2b + T_2), fluid, m_dot, d_i)
alpha_sc_o = alpha_1P_annulus(1e5, 0.5 * (T_sc_in + T_sc_out), 'water', m_sc_dot, d_ai, d_a)
alpha_sh_o = alpha_1P_annulus(1e5, T_2, 'thermooil', m_sh_dot, d_ai, d_a)
alpha_ws_o = alpha_1P_annulus(1e5, 0.5 * (T_ws_in + T_ws_out), 'water', m_ws_dot, d_ai, d_a)
alpha_ws_i = ht.condensation.Shah(m_dot, 0.5, d_i, rho_l, mu_l, k_l, cp_l, p_c, p_crit)

rho_l_evap = CP.PropsSI("D", "P", p_o, "Q", 0, fluid)
rho_v_evap = CP.PropsSI("D", "P", p_o, "Q", 1, fluid)
mu_l_evap = CP.PropsSI("VISCOSITY", "P", p_o, "Q", 0, fluid)
mu_v_evap = CP.PropsSI("VISCOSITY", "P", p_o, "Q", 1, fluid)
k_l_evap = CP.PropsSI("L", "P", p_o, "Q", 0, fluid)
cp_l_evap = CP.PropsSI("C", "P", p_o, "Q", 0, fluid)
h_evap = CP.PropsSI("H", "P", p_o, "Q", 1, fluid) - CP.PropsSI("H", "P", p_o, "Q", 0, fluid)
sigma_evap = CP.PropsSI("I", "P", p_o, "Q", 0, fluid)
T_excess = T_pinch + T_4b
dPsat = CP.PropsSI("P", "T", T_excess, "Q", 0, fluid) - p_o

alpha_evap_i = ht.boiling_flow.Chen_Bennett(m_dot, 0.5, d_i, rho_l_evap, rho_v_evap, mu_l_evap, mu_v_evap, k_l_evap
                                            , cp_l_evap, h_evap, sigma_evap, dPsat, T_excess)
alpha_evap_o = alpha_1P_annulus(1e5, 0.5 * (T_evap_in + T_evap_out), 'methanol', m_evap, d_ai, d_a)


def cal_U(alpha_i, alpha_o, d_a, d_ai, lam_tube):
    U = (1 / alpha_i + np.log(d_a / d_ai) * d_ai / (2 * lam_tube) + d_ai / (alpha_o * d_a)) ** (-1)
    return U


lam_tube = 300
U_sc = cal_U(alpha_sc_i, alpha_sc_o, d_a, d_ai, lam_tube)
U_sh = cal_U(alpha_sh_i, alpha_sh_o, d_a, d_ai, lam_tube)
U_ws = cal_U(alpha_ws_i, alpha_ws_o, d_a, d_ai, lam_tube)
U_evap = cal_U(alpha_evap_i, alpha_evap_o, d_a, d_ai, lam_tube)

delta_T_log_sc = T_log(T_2c, T_3, T_sc_in, T_sc_out)
delta_T_log_sh = T_log(T_2, T_2b, T_sh_in, T_sh_out)
delta_T_log_ws = T_log(T_2b, T_2c, T_ws_in, T_ws_out)
delta_T_log_evap = T_log(T_4, T_1, T_evap_in, T_evap_out)

A_sc = m_dot * delta_h_sc / delta_T_log_sc / U_sc
A_sh = m_dot * delta_h_sh / delta_T_log_sh / U_sh
A_ws = m_dot * delta_h_ws / delta_T_log_ws / U_ws
A_evap = m_dot * delta_h_o / delta_T_log_evap / U_evap
l_sc = A_sc / np.pi / d_i
l_sh = A_sh / np.pi / d_i
l_ws = A_ws / np.pi / d_i
l_evap = A_evap / np.pi /d_i
print(f"length for subcooled: {l_sc} m")
print(f"length for superheated: {l_sh} m")
print(f"length for wetsteam: {l_ws} m")
print(f"length for evaporator: {l_evap} m")

point_label = ["1", "2", "2b", "2c", "3", "4", "4b", "4c"]
x = [s_1, s_2, s_2b, s_2c, s_3, s_4, s_4b, s_4c]

y = [T_1, T_2, T_2b, T_2c, T_3, T_4, T_4b, T_4c]
plt.figure(1)
for i in range(len(x)):
    plt.plot(x[i], y[i], '*', markersize=15)
    plt.annotate(point_label[i], (x[i]+25, y[i]), fontsize=12)
plt.xlabel("s in J/kg/K")
plt.ylabel("T in K")
plt.legend()

# Berechnung des Nassdampfbereichs #
s_i = []
s_j = []
t_step = np.linspace(200, 442, 50)
for t_i in t_step:
    s_i1 = CP.PropsSI('S', 'T', t_i, 'Q', 0, fluid)
    s_i2 = CP.PropsSI('S', 'T', t_i, 'Q', 1, fluid)
    s_i.append(s_i1)
    s_j.append(s_i2)

plt.plot(s_i, t_step, 'k-')
plt.plot(s_j, t_step, 'k-', label="wet steam region")
#plt.xlabel('s in J/kg/K')
#plt.ylabel('T in K')
plt.title('T-s-Diagramm für ' + fluid)

# Berechnung Isobare #
s_step = np.linspace(200, 2700, 100)
for px in [p_o, p_c]:
    t_isobar = []
    for si in s_step:
        t_iso = CP.PropsSI('T', 'S', si, 'P', px, fluid)
        t_isobar.append(t_iso)

    plt.plot(s_step, t_isobar, 'b:', label="isobare")
plt.legend()

# h_dot-T diagram
x2 = np.array([h_1, h_2, h_2b, h_2c, h_3, h_4, h_4b, h_4c])
x2 = m_dot * x2
plt.figure(2)
for i in range(len(x2)):
    plt.plot(x2[i], y[i], '*', markersize=15)
    plt.annotate(point_label[i], (x2[i]+25, y[i]), fontsize=12)
plt.xlabel("h_dot in J/s")
plt.ylabel("T in K")
plt.legend()

# Berechnung des Nassdampfbereichs #
h_i = []
h_j = []
t_step = np.linspace(250, 442, 50)
for t_i in t_step:
    h_i1 = CP.PropsSI('H', 'T', t_i, 'Q', 0, fluid)
    h_i2 = CP.PropsSI('H', 'T', t_i, 'Q', 1, fluid)
    h_i.append(h_i1)
    h_j.append(h_i2)

plt.plot(np.array(h_i) * m_dot, t_step, 'k-')
plt.plot(np.array(h_j) * m_dot, t_step, 'k-', label="wet steam region")
#plt.xlabel('h in J/kg')
#plt.ylabel('T in K')
plt.title('T-h-Diagramm für ' + fluid)

# Berechnung Isobare #
h_step = np.linspace(0, 700000, 100)
for px in [p_o, p_c]:
    t_isobar = []
    for hi in h_step:
        t_iso = CP.PropsSI('T', 'H', hi, 'P', px, fluid)
        t_isobar.append(t_iso)

    plt.plot(h_step * m_dot, t_isobar, 'b:', label="isobare")
plt.legend()

# adding secondary fluids to plot figure 2
x_sec_evap = np.linspace(h_4 * m_dot, (h_4 + delta_h_o) * m_dot, 100)
y_sec_evap = np.linspace(T_evap_out, T_evap_in, 100)
plt.plot(x_sec_evap, y_sec_evap, 'b')

x_sec_sc = np.linspace(h_3 * m_dot, (h_3 + delta_h_sc) * m_dot, 100)
y_sec_sc = np.linspace(T_sc_in, T_sc_out, 100)
plt.plot(x_sec_sc, y_sec_sc, 'r')

x_sec_ws = np.linspace(h_2c * m_dot, (h_2c + delta_h_ws) * m_dot, 100)
y_sec_ws = np.linspace(T_ws_in, T_ws_out, 100)
plt.plot(x_sec_ws, y_sec_ws, 'r')

x_sec_sh = np.linspace(h_2b * m_dot, (h_2b + delta_h_sh) * m_dot, 100)
y_sec_sh = np.linspace(T_sh_in, T_sh_out, 100)
plt.plot(x_sec_sh, y_sec_sh, 'r')



plt.show()

