# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from plotSnapshot import plotSnapshot
from loadBinData  import loadSnapshot, loadParameters

scaleIndex = 0
timeIndex  = 100

paramFile = "params4python.dat"
nx, ny, nt, offset = loadParameters(paramFile)
nx = int(nx)
ny = int(ny)
nt = int(nt)

baseNames = ["heterogeneity", "ruptureTimes", "slipHistories", "offPlaneStress"]
titles2by2 = ["Heterogeneity", "Rupture time", "Slip", "Offplane Stress"]
shapes = [(nx, ny), (nx, ny, nt), (nx, ny, nt), (nx, ny, nt)]

hetero, rupTimes, slip, offPlaneStress = loadSnapshot(scaleIndex, timeIndex, baseNames, shapes)
plotSnapshot(hetero, rupTimes, slip, offPlaneStress, titles2by2, globalTitle=f"Offset = {offset:.3f}")






