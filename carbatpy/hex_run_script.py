# -*- coding: utf-8 -*-
"""
Example file for counter flow heat exchamnger evaluation
Reading from an input yaml file
calculating
plotting and writimg the results to a second file

(within carbatpy)
Created on Sat May  6 10:02:19 2023

@author: atakan 
"""

import heat_exchanger as h_ex
import pandas as pd



# Name of the yaml-file
y_file_name= "st_hex_parameters_file"
T0 =298.15

# read the yaml file as input
neu = h_ex.st_heat_exchanger_input.read_yaml(y_file_name+".yaml")

# new heat exchanger with the parameters from the file
hex2 = h_ex.counterflow_hex(*neu.all_out())

# solve the boundary value problem for the heat exchanger
res = hex2.he_bvp_solve()
print(f"Solution found: {res.success},  {res.message}")

#plotting and evalution:
f1, f2, ds, dq = hex2.he_state(res, 6, y_file_name)  # evaluate results (and plot)

# also adding the input to the Excel file as a new sheet
zzz = pd.DataFrame(dict( (key, value) for (key, value) in neu.__dict__.items()))
with pd.ExcelWriter(y_file_name+".xlsx", mode ="a") as writer: 
    zzz.to_excel(writer, sheet_name="input")


ex_in = hex2.exergy_entering()
print("Entropy production rate: %2.2e W/K, exergy loss rate %3.3f W, dq %3.2f"
      % (ds, ds * T0, dq))
print("Exergy flow rate, entering: %3.3f W" % (ex_in))
