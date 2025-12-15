# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path # import robust path handling
from loading.runData import runData # import data object
from loading.variableRegistry import VARIABLES # import variable scheme
from loading.loaders  import loadParameters # import loading function for run parameters
from plotting.plotSnapshot import plotContour2x2, plotProfiles2x2 # import plotting functions

PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/offPlaneStress/").resolve() # defines root directory for project
RUNSDIR     = PROJECTROOT / "runs" # path to directory containing all runs

runs = {} # prepare dictionary to hold all runs
for runDir in RUNSDIR.iterdir(): # iterates over all runs in the RUNSDIR directory
    if runDir.is_dir(): # checks if an object is a run folder
        paramFile = runDir / "params4python.dat" # prepare path to the parameter file
        nx, ny, nt, offset = loadParameters(paramFile) # loads parameters for specific run
        nx = int(nx) # turn parameters into ints
        ny = int(ny)
        nt = int(nt)
        runs[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=nt) # now use data object to fill the runs directionary
        # runs stored as objects with name, location, variables, dimensions in x, y, t directions
     
runName = "offset=20" # name of selected run
scaleIndex = 1 # scale index
timeIndex  = 100 # time index

run = runs[runName] # access run from dictionary
hetero         = run.load("heterogeneity", scaleIndex) # load heterogeneity map
rupTimes       = run.load("ruptureTimes", scaleIndex, timeIndex) # load rupture times and other time dependent variables
slip           = run.load("slipHistories", scaleIndex, timeIndex)
offPlaneStress = run.load("offPlaneStress", scaleIndex, timeIndex)
        
titles2by2 = ["Heterogeneity", "Rupture time", "Slip", "Offplane Stress"] # set titles for plotting

plotContour2x2(hetero, rupTimes, slip, offPlaneStress, titles2by2, globalTitle=f"Offset = {offset:.3f}") # make 2x2 contour map
plotProfiles2x2(hetero, rupTimes, slip, offPlaneStress, "y", 32, titles2by2, globalTitle=f"Offset = {offset:.3f}") # make line plot






