# TREND 5.0 Python Interface 
# Example Call
# @authors: David Celný , Sven Pohl
# Bochum, 16.05.2019
from fluid import *
import numpy as np
import time 

# For path seperation use D:\\...\\
path = 'C:\\Users\\atakan\\sciebo\\Progs\\TREND_50'

# Path to the TREND 5.0 DLL
dll_path = 'C:\\TREND_MIX\\TREND 5.0\\TREND_x64.dll'
dll_path = 'C:\\Users\\atakan\\sciebo\\Progs\\TREND_50\TREND_x64.dll'
_props ="dummy"
# Create object of fluid
fld = fluid('TP','H',['water'],[1],[1],1,path,'specific',dll_path)

# Calculate density with temperature and pressure input
T =373.15
p =1.0135
h, err = fld.TREND_EOS(T, p)
print('enth: ',h)
print('Error: ',err)


def hp(h, p, _fluid, composition=[1.0], option=1, 
           units ='specific', props=_props,
           eostype=[1],mixtype=1,path=path,dll_path=dll_path):
    p=p/1e5
    all_val=[]
    if option ==1:
        request = ["T", "D", "S", "PHASE"]
    elif option ==0:
        request = ["T", "D", "S", "PHASE","ETA", "TCX"]
    for what in request:
            hp_fld =fluid('HP',what,_fluid,composition,eostype,mixtype,path,
                          units,dll_path)
            value, err = hp_fld.TREND_EOS(h,p)
            
            if err.value ==0:
                all_val.append(value)
            else: print("Fehler", what, err, value)
            
    return np.array([all_val[0], p*1e5, h, 1/all_val[1], *all_val[2:]])

if __name__ == "__main__":
    t0 = time.time()
    b = hp(h, p * 1e5,["pentane","propane"], composition=[.50, .5], eostype=[1, 1], mixtype=1)
    print(b, time.time()-t0)
            