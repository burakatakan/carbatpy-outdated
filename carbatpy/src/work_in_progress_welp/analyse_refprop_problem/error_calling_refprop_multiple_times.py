from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
import numpy as np
from scipy.integrate import solve_ivp

os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum
SI = RP.GETENUMdll(0, "SI").iEnum

call = int(1e4)
checkpoints = np.linspace(0, call, 101)

def get_state(p, q, fluid, z=[1.0]):
    rho = RP.REFPROPdll(fluid, "PQ", "D", MASS_BASE_SI, 0, 0, p, q, z).Output[0]
    return rho

def try_fun(x, T, p, q, fluid, z=[1.0]):
    global counter_refprop
    res = RP.REFPROPdll(fluid, "PQ", "D;CP;VIS;TCX;PRANDTL;T;H", MASS_BASE_SI, 0, 0, p, q, z).Output[0:7]
    h = res[6]
    rho = res[0]
    counter_refprop +=1
    if rho == -9999990.:
        print(f"counter_refprop: {counter_refprop}")
        raise ValueError("Unplausible state")
    return [-0.5 * rho * T[0], h]

counter_refprop = 0
count = 0
for i in range(call):
    count +=1
    sol = solve_ivp(lambda x,T: try_fun(x, T, 1e6, 0.5, "water", [1.0]), [0,100], [4,200000])
    if (count%(call/100)) == 0:
        print(f"{count*100/call} %")

print(f"ivp part finished successfully after calling refprop {counter_refprop} times")

#count = 0
#for i in range(call):
#    rho_i = get_state(1e6, 0.5, "water")
#    count += 1
#    if rho_i == -9999990.:
#        print(f"counter_refprop: {counter_refprop}")
#        raise ValueError("Unplausible state")
#    if (count % (call / 100)) == 0:
#        print(f"{count * 100 / call} %")

#print(f"Script finished successfully after {call} calls.")