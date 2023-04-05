import numpy as np
import os
import matplotlib.pyplot as plt
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum
SI = RP.GETENUMdll(0, "SI").iEnum

x = np.linspace(0,1,11)
fluid = "PROPANE;ISOBUTANE"
h_v_store = []
for i in x:
    comp = [1-i,i]
    h_0 = RP.REFPROPdll(fluid, "PQ","H",MASS_BASE_SI,0,0,2e6,0,comp).Output[0]
    h_1 = RP.REFPROPdll(fluid,"PQ","H",MASS_BASE_SI,0,0,2e6,1,comp).Output[0]
    h_v = h_1 - h_0
    h_v_store.append(h_v)

plt.plot(x,h_v_store)
plt.show()
