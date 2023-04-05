"""
calling semi-empirical Roskosch compressor to check accordance with polynoms
author: Alexandra Welp
"""

import compressor_roskosch_orig_rp as comros
import matplotlib.pyplot as plt
import numpy as np
import bitzer_refprop_interface_ba_copy_welp as bitzer
from fl_props_compressor import z_uv, z_ps, z_Tp, z_Tx, z_mm
import fluid_properties_rp as fp

_props = "REFPROP"
fluid = "Propane"
comp = [1.0]
name_compressor = "2CESP-3P.xlsx"
Te_min, Te_max, Tc_min, Tc_max, d, h, f = bitzer.get_data(name_compressor)
points = 10
Te = np.linspace(Te_min, Te_max, points)+273.15
Tc = np.linspace(Tc_min, Tc_max, points)+273.15
V_h = 0.25 * np.pi * d ** 2 * h
T_in = 20 + 273.15
IS = 360  # number of differential steps for one cycle
IS0 = IS
pZyk = np.zeros(2, float)
z_it = np.zeros([IS, 16])
pZ = np.zeros(7, float)
###############################     parameter set specific for compressor   ###################################

pV = [d, h, 3.5, .04, .06071, 48.916, 50., 50. / 2., 2.]  # parameter see above
fo = open("Daten.txtx", "w")
#################################################################################################################
count = 0
for te in Te:
    eta_s = []

    eta_v = []
    pe = fp.T_prop_sat(te, fluid, comp)[1][1]/1000  # fl.zs_kg(['T','q'],[0.,0.],['p'],fluid)[0]
    for tc in Tc:
        pa = fp.T_prop_sat(tc, fluid, comp)[1][1]/1000  # fl.zs_kg(['T','q'],[35.,0.],['p'],fluid)[0]
        #fo = open("Daten.txtx", "w")
        #print("Drücke %2.2f kPa %2.2f kPa" % (pe, pa))
        #print(f"{T_in} {pe} {pa} {fluid} {comp} {pV} {pZ} {z_it} {IS} {pZyk} {IS0}")
        o1 = comros.getETA(T_in, pe, pa, fluid, comp, pV, pZ, z_it, IS, pZyk, IS0)
        eta_s.append(o1[0])
        eta_v.append(o1[1])
        count += 1
        print(count)
        print(o1)

    fo.write(str(eta_s, eta_v))


    plt.figure(1)
    plt.plot(Tc, eta_s)
    plt.figure(2)
    plt.plot(Tc, eta_v)

plt.show()
fo.close()
