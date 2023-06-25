# -*- coding: utf-8 -*-
"""
Using pymoo for the optimization of an  Carnot cycle with irreversible
heat transfer (twice dT). Objective functions are power and efficiency
The total heat exchanger area (times overall heat transfer coefficient)
shall be divided between the two heat exchangers.

Created on Wed Jan  4 13:36:11 2023

@author: atakan
"""

import numpy as np

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from scipy.optimize import root
from fluid_screening import dewpoints, v_names
import fluid_props as fprop
import pandas as pd
import seaborn as sns


#  Variables:
fixed_points = {
    "p_low": 1e5,
    "T_hh": 370,
    "T_hl": 290,
    "T_lh": 290,
    "T_ll": 276.0,
    "eta_s": 0.65
    
}
FLUID = "Ethane * Propane * Isobutane * Butane * Pentane"
FLS = "Methanol"  # "Water"  #
comp = [.10, 0.4, .19, 0.2, .01]
flm = fprop.FluidModel(FLUID)
myFluid = fprop.Fluid(flm, comp)
min_max = np.array([-1, -1, -1, -1, -1,1, -1,  1, -1, -1])

#  Function to approximate (Carnot with irreversible heat transfer):


class MyProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(n_var=4,
                         n_obj=10,
                         n_ieq_constr=3,
                         xl=np.array([0.00, 0.01, 0.01, 0.00]),
                         xu=np.array([0.1, 0.9, 0.9, 0.9]))

    def _evaluate(self, x, out, *args, **kwargs):
        f_out = dewpoints(x, myFluid, fixed_points)

        # print(f_out)

        # the total T difference between Th and Tl should not be reverted
        g1 = (x[0] + x[1] + x[2] + x[3]) - 1
        g2 = -f_out[-3] + 1e5
        g3 = -25e5 + f_out[-2]

        out["F"] = f_out * min_max  # minimization: minus
        out["G"] = [g1, g2, g3]


# new_in = [.3, .3,]

problem = MyProblem()

algorithm = NSGA2(pop_size=100,
                  eliminate_duplicates=True)

res = minimize(problem,
               algorithm,
               ("n_gen", 200),
               verbose=True,
               save_history=True,
               seed=1)

plot = Scatter()
plot.add(-res.F, edgecolor="red", facecolor="none")
plot.save("fluidOpti.png")
plot.show()
# for ii in range(2):
#     print(f"Mittelwerte-Obj.F {-res.F.mean(axis=0)[ii]:.2f}")
#     print(f"Standardabweichung-Obj.F {res.F.std(axis=0)[ii]:.2f}")
# for ii, n in enumerate(v_names):
#     print(
#         f"{n} - mean : {res.X.mean(axis=0)[ii]:.2f}, std-dev: {res.X.std(axis=0)[ii]:.1e}")


pareto = pd.DataFrame( min_max * res.F, columns=v_names)
paretox = pd. DataFrame(res.X, columns=FLUID.split(" * ")[:-1])
paretofx = pd. DataFrame(-res.X.sum(axis=1)+1,
                         columns=[FLUID.split(" * ")[-1]])
pareto_combine = pd.concat([pareto, paretox, paretofx], axis=1)
T_hh = fixed_points["T_hh"]
fname = f"fluidOpti_{T_hh}_00"
pareto_combine.to_csv(fname+"csv")
# pareto_combine.to_csv(fname+"csv")

# Plotting --------------------------
sns.set_theme(style="whitegrid")

# Load the example planets dataset
planets = sns.load_dataset("planets")

cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
g = sns.relplot(
    data=pareto_combine,
    x="delt_ex_destruct_compr", y="delt_ex_destruct",
    hue="p_ratio", size="T_liq(p_low)",
    palette=cmap, sizes=(10, 200),
)
# g.set(xscale="log", yscale="log")
g.ax.xaxis.grid(True, "minor", linewidth=.25)
g.ax.yaxis.grid(True, "minor", linewidth=.25)
g.despine(left=True, bottom=True)

g.savefig(fname+".png")
