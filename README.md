# carbatpy
Modeling Carnot Batteries (Thermal Energy Storage)

## This is the old version, we migrated!
The new Version is here: https://git.uni-due.de/spp-2403/td-ude/carbatpy


This is a project aiming to model thermal energy storages using heat pumps for 
charging, organic Rankine cycles (ORC) for discharging and different kinds of 
storages.
For this, it is planned to use detailed fluid models (as implemented e.g. in 
REFPROP, CoolProp, or TREND ) and setting up systems which can either be steady 
state or (later) also unsteady.
Since this project just starts, do not expect too much.
It is aimed to have heat exchangers, machines and storages as compounds, which 
can be combined to different charging and dicharging configurations. For these, 
the energy balance, mass balance and further relations will be applied/solved.
Later on also thermo-economic calculations are planned.

For the beginning, the solution of the spatially resolved heat exchanger 
profiles, a  boundary value problem, and its irreversibility will be 
implemented. An optimization will follow. 


Burak Atakan, University of Duisburg-Essen, Germany

You can contact us at: batakan [a t ]uni-duisburg.de or atakan.thermodynamik.duisburg [ a t] gmail.com

Some documentation is found here: https://carbatpy.readthedocs.io/en/latest/index.html

