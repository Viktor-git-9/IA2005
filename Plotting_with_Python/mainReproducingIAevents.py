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

def analyticSlipProfile(deltaSigma, mu, a, r):
    D = []
    for i in range(0, len(r)):
        if r[i] <= a:
            D.append(24/(7*np.pi) * deltaSigma/mu * np.sqrt(a**2 - r[i]**2))
        elif r[i] > a:
            D.append(0)
    return np.array(D)

def compareSlipProfile(slipProfile, deltaSigma, mu, a):
    maxInd = np.argmax(slipProfile)
    halfProfile = slipProfile[maxInd:-1]
    
    r = np.linspace(0, (len(halfProfile)-1)*4, len(halfProfile))
    analyticProfile = analyticSlipProfile(deltaSigma, mu, a, r)
    
    return [halfProfile, analyticProfile, r]

moment2mag  = lambda M: (np.log10(M) - 9.1) / 1.5
    


#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/homogeneous_dcmap/")
#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/reproducing_IA_events/6_4_event116/")
#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/6_4_1000_rupTimes")
#PROJECTROOT = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/single_asperities/6_1_single/fieldData")
#PROJECTROOT  = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/9Maps_1/corrected_moment_2/data/fieldData")
PROJECTROOT  = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/3_bins_5_15_25_35_50/datafiles/data_30/30/fieldData")
RUNSDIR     = PROJECTROOT / "runs" # path to directory containing all runs

runs = {} # prepare dictionary to hold all runs
for runDir in RUNSDIR.iterdir(): # iterates over all runs in the RUNSDIR directory
    if runDir.is_dir(): # checks if an object is a run folder
        paramFile = runDir / "params4python.dat" # prepare path to the parameter file
        nx, ny, nt, offset, runTime, kernelTime = loadParameters(paramFile) # loads parameters for specific run
        nx = int(nx) # turn parameters into ints
        ny = int(ny)
        nt = int(nt)
        #runs[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=nt, offset=offset, runTime=runTime) # now use data object to fill the runs directionary
        runs[runDir.name] = runData(runDir.name, runDir, VARIABLES, nx=nx, ny=ny, nt=2, offset=offset, runTime=runTime)
        # runs stored as objects with name, location, variables, dimensions in x, y, t directions
#newOrder4runs = ["Dc_0", "Dc_0point5", "Dc_1", "Dc_2"]
#runs = {k: runs[k] for k in newOrder4runs}

### contour plots of specified run at specified time
scaleIndex = 0
runName = "2" # name of selected run
run = runs[runName]

#moment = run.load("moment", scaleIndex)
#momentRate = run.load("momentRate", scaleIndex)
#magnitude = run.load("magnitude", scaleIndex)

#plotTimeSteps = [10, 100, 215] # for event 7987
#plotTimeSteps = [10, 50, 100, 192] # for event 806
#plotTimeSteps = [300, 400, 500] #  for stage 1 nrn
#plotTimeSteps = [75, 100, 125] # for stage 1 rn
#plotTimeSteps = [10, 15, 20, 200]
#plotTimeSteps = [0, 1]
plotTimeSteps = [1]
#plotTimeSteps = [10, 500, 750, 1000]
#timeStepLabels = [f"{v}dt" for v in plotTimeSteps]
#timeStepLabels = ["Simulation start", "Simulation stop"]
timeStepLabels = ["Simulation stop"]
cbarLabels_rupTimes = ["Rupture times [dt]"] * len(plotTimeSteps)
cbarLabels_slipHis = ["Accumulated slip [m]"] * len(plotTimeSteps)
cbarLabels_slipVelo = ["Slip velo. [m/s]"] * len(plotTimeSteps)
cbarLabels_stress = ["Stress [MPa]"] * len(plotTimeSteps)
profileDir2 = "y"
temp_rupTimes = []
temp_slipHis = []
temp_slipVelo = []
temp_onPlaneStress  = []

for iTime in plotTimeSteps:
    temp_rupTimes.append(run.load("ruptureTimes", scaleIndex, iTime))
    temp_slipHis.append(run.load("slipHistories", scaleIndex, iTime)*0.004) # add correction factor to computed slip profiles
    
    #temp_slipVelo.append(run.load("slipVelocities", scaleIndex, iTime)[96:160,96:160])
    #temp_slipVelo.append(run.load("slipVelocities", scaleIndex, iTime))
    
    #temp_onPlaneStress.append(run.load("onPlaneStress", scaleIndex, iTime)[96:160,96:160])
    temp_onPlaneStress.append(run.load("onPlaneStress", scaleIndex, iTime))
    
    
rupture_extents = patch_extents(temp_rupTimes)
print(rupture_extents)

plotContours(temp_rupTimes, timeStepLabels, cbarLabels_rupTimes, globalTitle = "Rupture times for a large event")
plotContours(temp_slipHis, timeStepLabels, cbarLabels_slipHis, globalTitle="Accumulated Slip for a large event")
#plotContours(temp_slipVelo, timeStepLabels, cbarLabels_slipVelo, clims=[0, 0.1], globalTitle="Slip velo. for large event " + runName) 
plotContours(temp_onPlaneStress, timeStepLabels, cbarLabels_stress, globalTitle=f"Onplane Stress for large event {runName}", clims=[0, 6])

### plot Dc heterogeneity ###
#heterogeneity = run.load("heterogeneity", scaleIndex, iTime)
#plotContours([heterogeneity],  ["Dc"], globalTitle = "Dc heterogeneity", clims=[0, 5])

### plot moment, moment rate, moment magnitude ###
#plotProfiles([moment], ["time step", "seismic moment"], ["moment"])
#plotProfiles([momentRate[0:300]], ["time step", "moment release rate [N m /s]"], ["Moment release rate"], globalTitle="Moment release rate")
#plotProfiles([magnitude], ["time step", "magnitude"], ["magnitude"])

#plotProfiles([temp_slipHis[0][:,156]], ["x [m]", "slip"], ["slip"]) # center at y = 156

deltaSigma = 3*1e6 # [MPa]
mu = 32.4*1e9
a = 192 # radius in m
A = np.pi * a**2
#myProfiles = compareSlipProfile(temp_slipHis[1][:,156], deltaSigma, mu, a)

#plotProfiles(myProfiles, ["x [els]", "slip [m]"], ["Simulation", "Prediction"], globalTitle="Predicted and simulated slip profiles" ,lineInd = 49)
#plotProfiles([temp_onPlaneStress[1][:,156]], ["x [els]", "stress [N m]"], ["Simulation", "Prediction"])

# momentSimulation = 2*np.pi*mu * integrate.simpson(myProfiles[0] * myProfiles[2], myProfiles[2])
# momentPrediction = 2*np.pi*mu * integrate.simpson(myProfiles[1] * myProfiles[2], myProfiles[2])
# momentDifference = 2*np.pi*mu * integrate.simpson(np.abs(myProfiles[0] - myProfiles[1]) * myProfiles[2], myProfiles[2])
# momentControl    = 16/7 * deltaSigma * a**3

# magTestSimulation = moment2mag(momentSimulation)
# magTestPrediction = moment2mag(momentPrediction)
# magTesttControl   = moment2mag(momentControl)

# fig, ax = plt.subplots()
# ax.plot(myProfiles[2]/a, myProfiles[1]/np.max(myProfiles[1]), label="Prediction", linewidth=3)
# ax.plot(myProfiles[2]/a, myProfiles[0]/np.max(myProfiles[1]), label="Simulation", linewidth=3)
# ax.set_xlim([0, 1.5])
# plt.grid()
# plt.legend(fontsize="large")
# plt.title("Self-similar slip profiles")
# plt.xlabel("x / asperity radius")
# plt.ylabel("slip / max. slip")


