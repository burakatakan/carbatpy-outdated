"""
testcase for compressor_roskosch model

Verification values obtained before modification of global variables etc
but AFTER changing to fluid mixtures and refprop, validation of these parts still to be performed

author: Alexandra Welp
"""

from carbatpy import compressor_roskosch_orig_rp as comros
import numpy as np
import pytest


def test_overall_script():
    fluid = 'Propane * Butane'
    comp = [1.0, 0.]
    pe = comros.z_Tx(263, 0, fluid, comp)[1]  # fl.zs_kg(['T','q'],[0.,0.],['p'],fluid)[0]
    pa = comros.z_Tx(355, 0, fluid, comp)[1]  # fl.zs_kg(['T','q'],[35.,0.],['p'],fluid)[0]
    dt_all = np.linspace(9.5, 20.5, 3)
    out = []
    for dt in dt_all:
        o1 = comros.getETA(dt + 273.15, pe, pa, fluid, comp)
        # o1.append((np.max(z_it[:,11]) - np.min(z_it[:,11]) * pV[7]))  # Massenstrom
        out.append(o1)
        print(dt, o1)
    actual_out = [[0.69124745, 0.40855443],[0.69800568, 0.41455612],[0.70437112, 0.42013682]]
    assert round(pe, 2) == 343.57
    assert round(pa, 2) == 3241.82
    assert (round(out[0][0],8) == actual_out[0][0])
    assert (round(out[1][1],8) == actual_out[1][1])
