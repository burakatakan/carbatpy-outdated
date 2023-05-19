# -*- coding: utf-8 -*-
"""
Example script for a (simple) counter flow heat exchamnger evaluation
Reading from an input yaml file
calculating, plotting and writimg the results to new files and storing them with 
date and time in a results directory. The main result are now stored as an
excel file (see heat_exchanger.he_state())

(within carbatpy)
Created on Sat May  6 10:02:19 2023

@author: atakan 
"""

import heat_exchanger as h_ex
# import xl_read_heatexchanger_new as xlread
import pandas as pd
import os
import shutil
from datetime import datetime
import csv

directory_path = os.path.dirname(os.path.abspath(__file__))
resultsDir = "results"
logFile = "HeatExchangerCalculations_log.csv"

# Name of the excel-file
y_file_name0 = "heat_exchanger_input1"
y_file_name = os.path.join(directory_path, y_file_name0)
T0 = 298.15

# read the yaml file as input
neu = h_ex.st_heat_exchanger_input.read_hex_file(
    y_file_name+".xlsx", "HEXsimple")
# The file should be copied to a file with the same name and the date added as below
now = datetime.now()
dt_string = now.strftime("_%Y_%m_%d_%H_%M_%S")

resDir = os.path.join(directory_path, resultsDir)
exists = os.path.exists(resDir)
if not exists:
    os.makedirs(resDir)
res_file_name = y_file_name+dt_string

shutil.copy(y_file_name+".xlsx", res_file_name+".xlsx")

# new heat exchanger with the parameters from the file
hex2 = h_ex.counterflow_hex(*neu.all_out())

# solve the boundary value problem for the heat exchanger
res = hex2.he_bvp_solve()
print(f"Solution found: {res.success},  {res.message}")

# plotting and evalution:
# evaluate results (and plot)
# y_file_name+="calc"
f1, f2, ds, dq = hex2.he_state(res, 6, res_file_name)


# also adding the input to the Excel file as a new sheet
zzz = pd.DataFrame(dict((key, value) for (key, value) in neu.__dict__.items()))
with pd.ExcelWriter(res_file_name+".xlsx", mode="a") as writer:
    zzz.to_excel(writer, sheet_name="input")


ex_in = hex2.exergy_entering()
print("Entropy production rate: %2.2e W/K, exergy loss rate %3.3f W, dq %3.2f"
      % (ds, ds * T0, dq))
print("Exergy flow rate, entering: %3.3f W" % (ex_in))

# copying the results to a results directory, including date and time


dircont = os.listdir(directory_path)
for file in dircont:
    allparts = file.split(".")
    if len(allparts) > 1:
        nbase, nend = allparts
        if nbase == y_file_name0+dt_string:
            fname = nbase+dt_string+"."+nend
            resDir0 = os.path.join(resDir, fname)
            shutil.move(file, resDir0)


# write in logfile:
try:
    inp_dict = neu.__dict__.copy()
    inp_dict["filename"] = res_file_name

    with open(os.path.join(resDir, logFile), 'a') as csvfile:
        writer = csv.DictWriter(csvfile, inp_dict)

        writer.writerow(inp_dict)
except IOError:
    print("I/O error")
