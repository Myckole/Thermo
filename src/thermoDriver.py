""" thermoDriver.py : Big Data Do-er

Outputs three data sets. Fugacity coefficient vs y1 with constant P and T, Fugacity 
coefficent vs P with constant y1 and T, and Fugacity coefficient vs T with constant
y1 and P. This only works for binary mixtures.

This script requires that 'scipy', 'numpy', 'pandas', and 'matplotlib' be installed within
the python environment you're running it in with the exception of 'matplotlib' if figure and
jpg == False.

thermoDriver.py and thermo.py must be in the same directory and any folder named data before
initial run should not exist.
"""

import thermo

"""
Parameters
----------
p,t,y1 : float
    The constant values of pressure, temperature, and y1 respectively
thermo.txt : boolean
    If true, saves three .csv files in data/ containing each point that could be graphed
thermo.jpg : boolean
    If true, saves three .jpg files in data/ containing graphs
thermo.figure : boolean
    If true, opens three instances of matplotlib figures (they're easier to look at than .jpg or .txt)
thermo.species1 : string
    The name of species 1
thermo.a, thermo.b : list
    Lists containing VdW's coefficients in the order of [y1, y2].
"""

p=40*10**5
t=473
y1 = 0.6

thermo.txt = False
thermo.jpg = False
thermo.figure = True
thermo.species1 = '{methane}'
thermo.a = [0.1303, 0.939]
thermo.b = [4.31e-5 , 9.05e-5]
thermo.graph(p,t,y1)
