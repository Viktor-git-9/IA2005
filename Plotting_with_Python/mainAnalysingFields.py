import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path # import robust path handling
from loading.runData import runData # import data object
from loading.variableRegistry import VARIABLES # import variable scheme
from loading.loaders  import loadParameters # import loading function for run parameters
from plotting.plotContour import plotContours, plotContour2x2, plotContour2x2plus1 # import plotting functions
from plotting.plotProfiles import getProfiles, plotProfiles2x2, plotProfiles2x2plus1, plotProfiles
from plotting.plotMesh import plotMeshSlices
import scipy.integrate as integrate

mapIndex = 50
savePath = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/IA2005/Plotting_with_Python/fieldAnalysis/"
PROJECTROOT  = Path(f"/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/3_bins_5_15_25_35_50/datafiles/data_{mapIndex}/{mapIndex}/fieldData")
RUNSDIR     = PROJECTROOT / "runs" # path to directory containing all runs

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


runsField = {} # prepare dictionary to hold all runs
runsLine  = {} # setting up a second directory for the lineData, because of the nt conflict
# this is very ugly and kind of defeats the purpose of the run dictionaries. Should be fixed
# maybe by including 2 different nt for fieldData and lineData respectively?
for runDir in RUNSDIR.iterdir(): # iterates over all runs in the RUNSDIR directory
    if runDir.is_dir(): # checks if an object is a run folder
        paramFile = runDir / "params4python.dat" # prepare path to the parameter file
        nx, ny, nt, offset, runTime, kernelTime = loadParameters(paramFile) # loads parameters for specific run
        nx = int(nx) # turn parameters into ints
        ny = int(ny)
        nt = int(nt)
        #runs[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=nt, offset=offset, runTime=runTime) # now use data object to fill the runs directionary
        runsField[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=2, offset=offset, runTime=runTime)
        runsLine[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=nt, offset=offset, runTime=runTime)
        # runs stored as objects with name, location, variables, dimensions in x, y, t directions
        
scaleIndex = 0
iTime = 1
eventSizes = []
eventMoments = []
averagedStress = np.zeros((nx, ny))
fieldSize = nx * ny

runCount = 0
for runDir in RUNSDIR.iterdir():
    if runDir.is_dir():
        runField = runsField[runDir.name]
        runLine  = runsLine[runDir.name]
        ruptureTimes  = runField.load("ruptureTimes", scaleIndex, iTime)
        stressField   = runField.load("onPlaneStress", scaleIndex, iTime)
        momentHistory = runLine.load("moment", scaleIndex)
        
        eventSizes.append(fieldSize - np.count_nonzero(ruptureTimes == -1))
        averagedStress = averagedStress + stressField
        eventMoments.append(np.max(momentHistory))
        
        runCount = runCount + 1
     
averagedStress = averagedStress/runCount
np.savetxt(savePath + f"eventMoments/eventMoments_{mapIndex}.txt", eventMoments, delimiter=',', newline="\n")
np.savetxt(savePath + f"eventSizes/eventSizes_{mapIndex}.txt", eventSizes, delimiter=',', newline="\n")

plotContours([averagedStress], cbarLabels=["Stress [MPa]"], titles=[f"Averaged over {runCount} events."], globalTitle="Initial Stress")
fig101 = plt.gcf()
fig101.savefig(savePath + f"/stressField_{mapIndex}.png")
plt.close(fig101)