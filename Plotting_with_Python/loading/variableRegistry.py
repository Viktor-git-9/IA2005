"""
Define the ensemble of variables, as it is produced by the numerical model.
The names of the variables in this registry have to agree exactly with the
files that are to be read, minus the scale index at the end.
"""

from .variables import VariableSpec

VARIABLES = {
    "heterogeneity" : VariableSpec("heterogeneity", ("nx", "ny"), False),
    "ruptureTimes"  : VariableSpec("ruptureTimes", ("nx", "ny", "nt"), True),
    "slipHistories" : VariableSpec("slipHistories", ("nx", "ny", "nt"), True),
    "offPlaneStress": VariableSpec("offPlaneStress", ("nx", "ny", "nt"), True),
    "onPlaneStress" : VariableSpec("onPlaneStress", ("nx", "ny", "nt"), True),
    "slipVelocities": VariableSpec("slipVelocities", ("nx", "ny", "nt"), True),
    "moment"        : VariableSpec("moment", ("nt",), True), # comma turns it into tuple with 1 element.
    "momentRate"    : VariableSpec("momentRate", ("nt",), True),
    "magnitude"     : VariableSpec("magnitude", ("nt",), True)}
