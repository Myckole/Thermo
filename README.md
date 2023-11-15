# Thermo

Various scripts that can be used to find patterns/values of thermodynamic properties. Specifically developed for Rutgers Course 155:309 Thermodynamics 2 for Chemical Engineering.

## Author
* Michael Nguyen (https://github.com/Myckole)

## Dependencies
The following program(s) is/are necessary to have installed:
* Python 3

## Scripts

* `thermoDriver.py`: Outputs three data sets. Fugacity coefficient vs $y_1$ with constant $P$ and $T$, Fugacity 
coefficent vs $P$ with constant $y_1$ and $T$, and Fugacity coefficient vs $T$ with constant
$y_1$ and $P$. This only works for binary mixtures.

* `thermo.py`: Carries supplemental methods for `thermoDriver.py`.
* `GenFug.py`: A fugacity calculator in progress.