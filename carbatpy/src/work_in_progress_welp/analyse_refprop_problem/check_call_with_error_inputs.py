import numpy as np
import os

from scipy.integrate import solve_ivp
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum
SI = RP.GETENUMdll(0, "SI").iEnum

p = 10e5
h_wf = 663819.6890146584
wf = "Isobutane"
Zustand_wf_x = RP.REFPROPdll(wf, "PH", "D;CP;VIS;TCX;PRANDTL;T", MASS_BASE_SI, 0, 0, 10e5, h_wf, [0])

print(Zustand_wf_x)