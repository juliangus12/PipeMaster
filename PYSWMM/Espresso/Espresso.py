'''
PySWMM Code Example A
Author: Bryant McDonnell
Version: 1
Date: Oct 21, 2022
'''
from pyswmm import Simulation

with Simulation(r'Example1.inp') as sim:
    # Launch a simulation!
    for ind, step in enumerate(sim):
        if ind % 100 == 0:
            print(round(sim.percent_complete*100))
