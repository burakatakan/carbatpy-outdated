"""
created on 13.04.2023 to check if new version of fluid_properties now calculates q correctly
@author: welp
"""
import carbatpy.fluid_properties_rp as fprop

fluid = "Propane;Isobutane"
comp = [0.1, 0.9]
p = 1e6
h = 525e3

q = fprop.hp(h, p, fluid, comp)[5]

if q <= 0:
    raise ValueError("vapour quality was not calculated correctly")
else:
    print(f"q = {q}")
