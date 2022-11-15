# -*- coding: utf-8 -*-
"""
Reversibler Kondensator mit Druckabfall: wie muss der Verlauf von p, q und h sein?
Created on Tue Oct 11 10:58:53 2022

@author: atakan
"""

import fluid_properties_rp as cb
import numpy as np
import matplotlib.pyplot as plt
x_0 = 1.
p_0 = 20e5     # Anfangsdruck Pa
_props = "REFPROP"

temp_sur = 283.15
temp_0_s = 373.15
fluid_s = "Propane * Pentane"
comp = [1., 0.0]
state_props = []
p_all = np.linspace(p_0, 1e5, 2)
for n_p, p in enumerate(p_all):
    state_props.append(cb.p_prop_sat(p, fluid_s, composition=comp,
                                     props=_props))
p0, h0 = state_props[0][0, 1:3]
p1, h1 = state_props[-1][1, 1:3]

qT = []
no_points = 150
hs = np.linspace(h0, h1, no_points)
ps = np.linspace(p0, p1, no_points)
state_old = np.zeros((6))
for n_p, p in enumerate(ps):
    state = cb.hp(hs[n_p], ps[n_p], fluid_s, composition=comp,
                  props=_props)
    if n_p > 0:
        T_m = (state + state_old) / 2
        q = state - state_old
        ds = q[4]
        q = q[2]
        qT.append(np.array([T_m[0], q, ds, q/T_m[0]]))
    state_old = state
qT = np.array(qT)
fig, ax = plt.subplots(1, 3)
hs = hs*1e-3
ax[0].plot(hs[1:], qT[:, 0])
ax[0].set_xlabel("h")
ax[0].set_ylabel("T_mean")
ax[1].plot(hs[1:], qT[:, 2] - qT[:, 3])
ax[1].set_xlabel("h")
ax[1].set_ylabel("Delta s_irr")
ax[2].plot(hs[1:], qT[:, 2])
ax[2].plot(hs[1:], qT[:, 3], "o")
ax[2].set_xlabel("h")
ax[2].set_ylabel("delta s and dq/T_mean")

""" naechste Frage: wie muss der Massenstrom vom Sekundärfluid und der Massenstrom
bzw. der Durchmesser in axialer Richtung verändert werden, damit das passt?"""
