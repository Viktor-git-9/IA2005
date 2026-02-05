# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path # import robust path handling
from loading.runData import runData # import data object
from loading.variableRegistry import VARIABLES # import variable scheme
from loading.loaders  import loadParameters # import loading function for run parameters
from plotting.plotContour import plotContours, plotContour2x2, plotContour2x2plus1 # import plotting functions
from plotting.plotProfiles import getProfiles, plotProfiles2x2, plotProfiles2x2plus1, plotProfiles
from plotting.plotMesh import plotMeshSlices


#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/offPlaneStress/heterogeneous/").resolve() # defines root directory for project
#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/offPlaneStress/homogeneous/nmax64").resolve() # defines root directory for project
PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/timing")
RUNSDIR     = PROJECTROOT / "runs" # path to directory containing all runs

runs = {} # prepare dictionary to hold all runs
for runDir in RUNSDIR.iterdir(): # iterates over all runs in the RUNSDIR directory
    if runDir.is_dir(): # checks if an object is a run folder
        paramFile = runDir / "params4python.dat" # prepare path to the parameter file
        nx, ny, nt, offset, runTime, kernelTime = loadParameters(paramFile) # loads parameters for specific run
        nx = int(nx) # turn parameters into ints
        ny = int(ny)
        nt = int(nt)
        runs[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=nt, offset=offset, runTime = runTime, kernelTime = kernelTime) # now use data object to fill the runs directionary
        # runs stored as objects with name, location, variables, dimensions in x, y, t directions
newOrder4runs = ["nmax_32", "nmax_64", "nmax_128", "nmax_256", "nmax_512"]
runs = {k: runs[k] for k in newOrder4runs}

scaleIndex = 0 # scale index for plots
runTimes = []
kernelTimes = []
nmax = []
kernelSizes = []
for iRun in runs.values():
    runTimes.append(iRun.runTime)
    kernelTimes.append(iRun.kernelTime)
    nmax.append(iRun.sizes["nx"])
    # calculate kernel sizes as: dimensions of kernel * FFT factor (=2*2) * bits per element * conversion constant
    kernelSizes.append(iRun.sizes["nx"]*iRun.sizes["ny"]*iRun.sizes["nt"]*16*10**(-9)*4)
    
elementCount = np.power(nmax,2)
scatterLabels = newOrder4runs

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot()
plt.scatter(elementCount, kernelSizes, s=100)
plt.xlabel('element count', fontsize=18)
plt.ylabel('Kernel memory size [GB]', fontsize=18)
plt.title('kernel size over # elements', fontsize=18)
plt.grid()
for i, txt in enumerate(scatterLabels):
    if i == 0:
        ax.annotate(txt, (elementCount[i], kernelSizes[i]-0.38))
    else:
        ax.annotate(txt, (elementCount[i], kernelSizes[i]-0.3))

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot()
plt.scatter(elementCount, runTimes, s=100)
plt.xlabel('element count', fontsize=18)
plt.ylabel('total runtime [s]', fontsize=18)
plt.title('total runtime over # elements', fontsize=18)
plt.grid()
for i, txt in enumerate(scatterLabels):
    if i == 0:
       ax.annotate(txt, (elementCount[i], runTimes[i]-70))
    else:
       ax.annotate(txt, (elementCount[i], runTimes[i]-50))

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot()
plt.scatter(elementCount, np.divide(kernelTimes, runTimes), s=100)
plt.xlabel('element count', fontsize=18)
plt.ylabel('time [s]', fontsize=18)
plt.title('run time / kernel compilation time over # elements', fontsize=18)
plt.grid()
for i, txt in enumerate(scatterLabels):
    ax.annotate(txt, (elementCount[i], np.divide(kernelTimes, runTimes)[i]-0.01))


### Currently unused ###  
def plotStackedBars(x, y1, y2, color_bottom="tab:blue", color_top="tab:orange"):
    """
    Plot two-colored bars.

    Parameters
    ----------
    x:  list of numbers
        X positions of the bars.
    y1: list of numbers
        total heights of bars
    y2: list of numbers
        height of lower bar sections
    color_bottom : str
        Color of the lower segment.
    color_top : str
        Color of the upper segment.
    """
    
    topHeights = [a - b for a, b  in zip(y1, y2)]
    
    fig, ax = plt.subplots()

    # plot stacked bars
    ax.scatter(x, y1)
    ax.bar(x, y2, color = color_bottom)
    ax.bar(x, topHeights, bottom = y2, color = color_top)
    # plt.scatter(x, y1)
    # plt.bar(x, y1, bottom=y2, color=color_top)

    plt.xlabel("x")
    plt.ylabel("height")
    plt.tight_layout()
    plt.show()

