# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path # import robust path handling
from loading.runData import runData # import data object
from loading.variableRegistry import VARIABLES # import variable scheme
from loading.loaders  import loadParameters # import loading function for run parameters
from plotting.plotContour import plotContours, plotContour2x2, plotContour2x2plus1 # import plotting functions
from plotting.plotProfiles import getProfiles, plotProfiles2x2, plotProfiles2x2plus1, plotProfiles

#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/offPlaneStress/heterogeneous/").resolve() # defines root directory for project
PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/offPlaneStress/homogeneous/nmax64").resolve() # defines root directory for project
RUNSDIR     = PROJECTROOT / "runs" # path to directory containing all runs

runs = {} # prepare dictionary to hold all runs
for runDir in RUNSDIR.iterdir(): # iterates over all runs in the RUNSDIR directory
    if runDir.is_dir(): # checks if an object is a run folder
        paramFile = runDir / "params4python.dat" # prepare path to the parameter file
        nx, ny, nt, offset = loadParameters(paramFile) # loads parameters for specific run
        nx = int(nx) # turn parameters into ints
        ny = int(ny)
        nt = int(nt)
        runs[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=nt, offset=offset) # now use data object to fill the runs directionary
        # runs stored as objects with name, location, variables, dimensions in x, y, t directions
     
runName = "offset=5" # name of selected run
scaleIndex = 0 # scale index
timeIndex  = 10 # time index

run = runs[runName] # access run from dictionary
hetero         = run.load("heterogeneity", scaleIndex) # load heterogeneity map
rupTimes       = run.load("ruptureTimes", scaleIndex, timeIndex) # load rupture times and other time dependent variables
slip           = run.load("slipHistories", scaleIndex, timeIndex)
offPlaneStress = run.load("offPlaneStress", scaleIndex, timeIndex)
onPlaneStress  = run.load("onPlaneStress", scaleIndex, timeIndex)
runData = [hetero, rupTimes, slip, offPlaneStress, onPlaneStress]
        
runDataTitles = ["Heterogeneity", "Rupture time", "Slip", "Offplane Stress", "Onplane Stress"] # set titles for plotting
cbarLabels    = ["Dc [mm]", "time [dt]", "Slip [mm]", "Stress [MPa]", "Stress [MPa]"]

### Plotting contour plots of run selected above ###
plotContours(runData, runDataTitles, cbarLabels, globalTitle=f"Offset = {run.offset:.3f}")

### Plot profiles of run selected above along x or y axis ###
profiles = []
for iArr in runData:
    profiles.append(getProfiles(iArr, "x", 32))
 
#plotProfiles2x2plus1(profiles, titles2by2, globalTitle=f"Offset = {offset:.3f}") # make line plot

### Stress field profiles at different offsets ###
#offPlaneStresses = []
#for iRun in runs:
#    offPlaneStresses.append(getProfiles(runs[iRun].load("offPlaneStress", scaleIndex, timeIndex), "x", 32))
 
#plotProfiles(offPlaneStresses, runs, "Offplane stresses")

### Stress field contours at different time steps
scaleIndex2 = 1
plotTimeSteps = [10, 20, 30, 40, 50, 60, 70, 80]
timeStepLabels = [f"{v}dt" for v in plotTimeSteps]
cbarLabels2 = ["Stress [MPa]"] * len(plotTimeSteps)
profileDir2 = "y"
temp_rupTimes = []
temp_slip = []
temp_offPlaneStress = []
temp_onPlaneStress  = []

for iTime in plotTimeSteps:
    temp_rupTimes.append(run.load("ruptureTimes", scaleIndex2, iTime))
    temp_slip.append(run.load("slipHistories", scaleIndex2, iTime))
    temp_offPlaneStress.append(run.load("offPlaneStress", scaleIndex2, iTime))
    temp_onPlaneStress.append(run.load("onPlaneStress", scaleIndex2, iTime))
    
plotContours(temp_onPlaneStress, timeStepLabels, cbarLabels2, globalTitle="Onplane Stresses")
plotContours(temp_offPlaneStress, timeStepLabels, cbarLabels2, globalTitle="Offplane Stresses")


### Stress field profiles at different times
temp_ProfrupTimes = []
temp_Profslip = []
temp_ProfoffPlaneStress = []
temp_ProfonPlaneStress  = []
for iTime in range(len(plotTimeSteps)):
    temp_ProfrupTimes.append(getProfiles(temp_rupTimes[iTime], profileDir2, 32))
    temp_Profslip.append(getProfiles(temp_slip[iTime], profileDir2, 32))
    temp_ProfoffPlaneStress.append(getProfiles(temp_offPlaneStress[iTime], profileDir2, 32))
    temp_ProfonPlaneStress.append(getProfiles(temp_onPlaneStress[iTime], profileDir2, 32))

axesLabels = ["x [element count]", "Stress [MPa]"]

plotProfiles(temp_ProfonPlaneStress, axesLabels, timeStepLabels, "Onplane stress at different times")
plotProfiles(temp_ProfoffPlaneStress, axesLabels, timeStepLabels, "Offplane stress at different times")











