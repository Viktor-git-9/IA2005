import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path # import robust path handling
from loading.runData import runData # import data object
from loading.variableRegistry import VARIABLES # import variable scheme
from loading.loaders  import loadParameters # import loading function for run parameters
from plotting.plotContour import plotContours, plotContour2x2, plotContour2x2plus1 # import plotting functions
from plotting.plotProfiles import getProfiles, plotProfiles, plot2Profiles
from plotting.plotMesh import plotMeshSlices

#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/homogeneous_dcmap/")
PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/806/")
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
stageChangeInd = None # index of move from stage 0 to stage 1
scaleIndex1 = 0
scaleIndex2 = 0
cutOffInd1 = 220
cutOffInd2 = 220
runName1 = "renorm_on" # names of selected runs
runName2 = "renorm_off"
correctionFactor = 1
plotTitle = "Event 806, Stage 0"

run_on = runs[runName1]
run_off = runs[runName2]

moment_on = run_on.load("moment", scaleIndex1)[0:cutOffInd1]
momentRate_on = run_on.load("momentRate", scaleIndex1)[0:cutOffInd1]
magnitude_on = run_on.load("magnitude", scaleIndex1)[0:cutOffInd1]

moment_off = run_off.load("moment", scaleIndex2)[0:cutOffInd2]
momentRate_off = run_off.load("momentRate", scaleIndex2)[0:cutOffInd2]
magnitude_off = run_off.load("magnitude", scaleIndex2)[0:cutOffInd2]

### Plot moment, moment rate, magnitude ###
plot2Profiles([moment_on, moment_off*correctionFactor], xlabel_bottom="time step (rn on)", xlabel_top="time step (rn off)",
                    ylabel="Moment [N m]", title=plotTitle, labels=["rn on", "rn off"], lineInd=stageChangeInd)

plot2Profiles([momentRate_on, momentRate_off*correctionFactor], xlabel_bottom="time step (rn on)", xlabel_top="time step (rn off)",
                    ylabel="Moment Rate [N m / s]", title=plotTitle, labels=["rn on", "rn off"], lineInd=stageChangeInd)

plot2Profiles([magnitude_on, magnitude_off], xlabel_bottom="time step (rn on)", xlabel_top="time step (rn off)",
                    ylabel="Magnitude", title=plotTitle, labels=["rn on", "rn off"], lineInd=stageChangeInd)


