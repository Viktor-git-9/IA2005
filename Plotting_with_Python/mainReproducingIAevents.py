import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path # import robust path handling
from loading.runData import runData # import data object
from loading.variableRegistry import VARIABLES # import variable scheme
from loading.loaders  import loadParameters # import loading function for run parameters
from plotting.plotContour import plotContours, plotContour2x2, plotContour2x2plus1 # import plotting functions
from plotting.plotProfiles import getProfiles, plotProfiles2x2, plotProfiles2x2plus1, plotProfiles
from plotting.plotMesh import plotMeshSlices

def patch_extents(arrays):
    """
    Determine x- and y-extent of the non--1 patch in each 2D array.
    To be used on the rupture time patches.

    Parameters
    ----------
    arrays : list of 2D numpy arrays

    Returns
    -------
    extents : list of tuples
        Each tuple is (x_extent, y_extent)
    """
    extents = []

    for arr in arrays:
        y_idx, x_idx = np.where(arr != -1)

        if len(x_idx) == 0:
            # no patch found
            extents.append((0, 0))
            continue

        x_extent = x_idx.max() - x_idx.min() + 1
        y_extent = y_idx.max() - y_idx.min() + 1

        extents.append((x_extent, y_extent))

    return extents


#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/homogeneous_dcmap/")
PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/6_4_event116/")
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
scaleIndex = 0
runName = "myRun" # name of selected run
run = runs[runName]

#moment = run.load("moment", scaleIndex)
#momentRate = run.load("momentRate", scaleIndex)
#magnitude = run.load("magnitude", scaleIndex)

#plotTimeSteps = [10, 100, 215] # for event 7987
#plotTimeSteps = [10, 50, 100, 192] # for event 806
#plotTimeSteps = [300, 400, 500] #  for stage 1 nrn
#plotTimeSteps = [75, 100, 125] # for stage 1 rn
plotTimeSteps = [280]
#plotTimeSteps = [10, 500, 750, 1000]
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
    
    #temp_slipVelo.append(run.load("slipVelocities", scaleIndex, iTime)[96:160,96:160])
    #temp_slipVelo.append(run.load("slipVelocities", scaleIndex, iTime))
    
    #temp_onPlaneStress.append(run.load("onPlaneStress", scaleIndex, iTime)[96:160,96:160])
    temp_onPlaneStress.append(run.load("onPlaneStress", scaleIndex, iTime))
    
rupture_extents = patch_extents(temp_rupTimes)

plotContours(temp_rupTimes, timeStepLabels, ["Rupture times [dt]"], globalTitle="Rupture Times")
plotContours(temp_slipHis, timeStepLabels, ["Slip [m?]"], globalTitle="Acc. Slip")
#plotContours(temp_slipVelo, timeStepLabels, cbarLabels2, clims=[0, 0.5], globalTitle="Curr. Slip") 
#plotContours(temp_onPlaneStress, timeStepLabels, cbarLabels2, globalTitle="Onplane Stress for large event at density = 6", clims=[0, 6])

### plot Dc heterogeneity ###
#heterogeneity = run.load("heterogeneity", scaleIndex, iTime)
#plotContours([heterogeneity],  ["Dc"], globalTitle = "Dc heterogeneity", clims=[0, 5])

### plot moment, moment rate, moment magnitude ###
#plotProfiles([moment], ["time step", "seismic moment"], ["moment"])
#plotProfiles([momentRate], ["time step", "moment release rate"], ["momentRate"])
#plotProfiles([magnitude], ["time step", "magnitude"], ["magnitude"])
