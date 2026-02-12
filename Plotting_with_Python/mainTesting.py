import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path # import robust path handling
from loading.runData import runData # import data object
from loading.variableRegistry import VARIABLES # import variable scheme
from loading.loaders  import loadParameters # import loading function for run parameters
from plotting.plotContour import plotContours, plotContour2x2, plotContour2x2plus1 # import plotting functions
from plotting.plotProfiles import getProfiles, plotProfiles2x2, plotProfiles2x2plus1, plotProfiles
from plotting.plotMesh import plotMeshSlices

PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/homogeneous_dcmap/")
#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/537/")
RUNSDIR     = PROJECTROOT / "runs" # path to directory containing all runs

runs = {} # prepare dictionary to hold all runs
for runDir in RUNSDIR.iterdir(): # iterates over all runs in the RUNSDIR directory
    if runDir.is_dir(): # checks if an object is a run folder
        paramFile = runDir / "params4python.dat" # prepare path to the parameter file
        nx, ny, nt, offset, runTime, kernelTime = loadParameters(paramFile) # loads parameters for specific run
        nx = int(nx) # turn parameters into ints
        ny = int(ny)
        nt = int(nt)
        runs[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=nt, offset=offset, runTime=runTime) # now use data object to fill the runs directionary
        # runs stored as objects with name, location, variables, dimensions in x, y, t directions
#newOrder4runs = ["Dc_0", "Dc_0point5", "Dc_1", "Dc_2"]
#runs = {k: runs[k] for k in newOrder4runs}

### contour plots of specified run at specified time
scaleIndex = 1
runName = "renorm_on" # name of selected run
run = runs[runName]

#plotTimeSteps = [75, 76, 77, 78, 80, 90, 100, 200]
plotTimeSteps = [60, 70, 80, 90, 100, 120]
timeStepLabels = [f"{v}dt" for v in plotTimeSteps]
cbarLabels2 = ["Stress [MPa]"] * len(plotTimeSteps)
profileDir2 = "y"
temp_rupTimes = []
temp_slipHis = []
temp_slipVelo = []
temp_onPlaneStress  = []

for iTime in plotTimeSteps:
    temp_rupTimes.append(run.load("ruptureTimes", scaleIndex, iTime))
    temp_slipHis.append(run.load("slipHistories", scaleIndex, iTime))
    temp_slipVelo.append(run.load("slipVelocities", scaleIndex, iTime))
    temp_onPlaneStress.append(run.load("onPlaneStress", scaleIndex, iTime))

plotContours(temp_rupTimes, timeStepLabels, cbarLabels2, globalTitle="Rupture Times")
# plotContours(temp_slipHis, timeStepLabels, cbarLabels2, globalTitle="Acc. Slip")
plotContours(temp_slipVelo, timeStepLabels, cbarLabels2, globalTitle="Curr. Slip") 
plotContours(temp_onPlaneStress, timeStepLabels, cbarLabels2, globalTitle="Onplane Stress")

### plot Dc heterogeneity ###
heterogeneity = run.load("heterogeneity", scaleIndex, iTime)
plotContours([heterogeneity],  ["Dc"], globalTitle = "Dc heterogeneity")
