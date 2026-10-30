#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 13:10:27 2026

@author: viktor
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ------------------------------------------------------------------
# Data extraction
# ------------------------------------------------------------------
def getMagnitudes(dataPath, mapCount):

    all_files = sorted(os.listdir(dataPath))
    dataList = []
    
    for ifile, file in enumerate(all_files):
            if file.endswith('agnitudes.txt'): # changed this to retrieve magnitude files
                dataList.append(np.loadtxt(dataPath + file, delimiter=" "))
                
    return dataList

# ------------------------------------------------------------------
# Magnitude analysis
# ------------------------------------------------------------------
# def runaway_probability(arr, rel_tol=0.02, runaway_magnitude=None):
#     a = np.asarray(arr, dtype=float).ravel()
#     a = a[~np.isnan(a)]
#     m = a.max() if runaway_magnitude is None else runaway_magnitude
#     close = np.abs(a - m) <= rel_tol * abs(m)      # entries within ±2 % of m
#     if runaway_magnitude is None and close.sum() - 1 < 2 or m < 3:
#         print("No runaway magnitude detected, returning 0.")
#         return 0.0                                 # fewer than 2 others near the max
    
#     if runaway_magnitude is None:
#         print(f"Runaway magnitude detected at {np.round(m,3)}")
#     else:
#         print(f"Runaway magnitude set to {np.round(m,3)}")
        
        
#     return close.mean()

def runaway_probability(arr, rel_tol=0.05, runaway_magnitude=None,
                        gap_tol=0.1, return_gap=False):
    a = np.asarray(arr, dtype=float).ravel()
    a = a[~np.isnan(a)]
    m = a.max() if runaway_magnitude is None else runaway_magnitude
    close = np.abs(a - m) <= rel_tol * abs(m)      # entries within ±2 % of m
    if runaway_magnitude is None and close.sum() - 1 < 2 or m < 3:
        print("No runaway magnitude detected, returning 0.")
        print("Gap size set to NaN")
        return (0.0, np.nan) if return_gap else 0.0

    if runaway_magnitude is None:
        print(f"Runaway magnitude detected at {np.round(m,3)}")
    else:
        print(f"Runaway magnitude set to {np.round(m,3)}")

    # --- gap detection ---
    below = a[~close & (a < m)]                    # entries below the runaway window
    if below.size > 0:
        gap = (m - below.max()) / abs(m)           # relative distance to next largest value
        has_gap = gap >= gap_tol
        print(f"Next largest value: {np.round(below.max(),3)} "
              f"(relative gap {gap:.1%}) -> gap {'found' if has_gap else 'not found'}")
    else:
        has_gap = False                            # nothing outside the window, no gap defined
        gap = np.nan
        print("No values outside the runaway window, gap size set to NaN.")

    prob = close.mean()
    return (prob, gap) if return_gap else prob
            
# ------------------------------------------------------------------
# MLE estimation & GR law
# ------------------------------------------------------------------
def MLE_b_value(magnitudes, cutoffMagnitude):
    MLE_b = np.log10(np.exp(1))/(np.mean(magnitudes) - cutoffMagnitude)
    MLE_b = MLE_b * (len(magnitudes) - 1)/len(magnitudes) # small sample bias correction
    MLE_a = np.log10(len(magnitudes)) + MLE_b * cutoffMagnitude # get a from MLE estimation for b
    
    return MLE_b, MLE_a

def MLE_powerlaw_exponent(data, cutoffPoint): # From Clauset, Shahizi and Newman (2009)
    cutoffData = data[data > cutoffPoint]
    
    n = len(cutoffData)
    sumTerm = np.sum(np.log(cutoffData/(cutoffPoint - 1/2)))
    alpha = 1 + n / sumTerm
    
    alpha_error = (alpha - 1) / np.sqrt(n)
    
    return alpha, alpha_error

def GRlaw(mags, cutoffMag, b):
    return 10**(-b*(mags-cutoffMag))

# ------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------
def plotHistogram(data, bins, figLabels = None, figTitle = None, figSubTitle = None, textCount = None):
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(1,1,1)
    ax.hist(data, bins)
    #ax.set_yscale("log", base=10)
    
    fig.supxlabel(figLabels[0])
    fig.supylabel(figLabels[1])
    
    if figTitle != 'None':
        fig.suptitle(figTitle, fontsize=16)
        plt.title(figSubTitle, fontsize = 12)
        
    #Add text
    if textCount != 'None':
        ax.text(
            0.5, 0.8,
            f"Nr. of interrupted large events: {textCount}",
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top')
        
    plt.grid()
    
# make a cumulative loglog plot of some data
def plotHistoCum(data, bins, figLabels = None, figTitle = None, figSubTitle = None, bVal = None):
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(1,1,1)
    for iData, currentData in enumerate(data): 
        ax.ecdf(currentData, complementary=True, linewidth = 3)
    ax.set_yscale("log", base=10)
    
    fig.supxlabel(figLabels[0])
    fig.supylabel(figLabels[1])
    
    if figTitle != 'None':
        fig.suptitle(figTitle, fontsize=16)
        plt.title(figSubTitle, fontsize = 12)
        
    #Add text
    #if bVal != 'None':
    #    ax.text(
    #        0.5, 0.8,
    #        f"MLE b-value estimation: {bVal}",
    #        transform=ax.transAxes,
    #        fontsize=8,
    #        verticalalignment='top')
        
    #ax.legend(["Predicted", "Simulated"])
    plt.grid()
    return ax
    
# ------------------------------------------------------------------
# Code starts here
# ------------------------------------------------------------------

# dataPath1 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/1006/100asps_Rinis/infofiles/asperityMagnitudes/"
# dataPath2 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/1006/100asps_Rinis/datafiles/simulationMagnitudes/"

#dataPath1 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/1006/200asps_4point9MPa_Rinis/infofiles/asperityMagnitudes/"
#dataPath2 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/1006/200asps_4point9MPa_Rinis/datafiles/simulationMagnitudes/"

# dataPath1 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/3_bins_5_15_25_35_50/infofiles/above10percent/"
# dataPath2 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/3_bins_5_15_25_35_50/datafiles/above10percent/"

# dataPath1 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/0923_stresstest/datafiles/data_4/4/"
# dataPath2 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/0923_stresstest/datafiles/data_3/3/"

targetfile = 30
dataPath1 = f"/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/phase_transition_1016/infofiles/{targetfile}/"
dataPath2 = f"/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/phase_transition_1016/datafiles/data_{targetfile}/{targetfile}/"

# dataPath1 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/IA2005/Plotting_with_Python/fieldAnalysis/eventSizes/"
# dataPath2 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/IA2005/Plotting_with_Python/fieldAnalysis/eventMoments/"

#dataPath_allInfo = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/IA2005/Plotting_with_Python/savedDcMaps/1005/orderedData/"
#dataPath_allInfo = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/1006/200asps_4point9MPa_Rinis/infofiles/"
dataPath_allInfo = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/1006/100asps_Rinis/infofiles/"
m_c = 1.0928
correctionTerm = 0.12
#correctionTerm = 0
#asperityCount = 5
mapCount = 10

# %% Load and treat predicted & simulated magnitudes
dataList1 = getMagnitudes(dataPath1, mapCount)
dataList2 = getMagnitudes(dataPath2, mapCount)
    
dataArray1 = np.concatenate(dataList1)
dataArray2 = np.concatenate(dataList2)

dataArray2 = dataArray2 - correctionTerm

## hotfix for messed up predicted magnitudes - event magnitudes relation ##
dataArray1 = dataArray1[-len(dataArray2):]

order = np.argsort(dataArray1)
dataArray1_ordered = dataArray1[order]
dataArray2_ordered = dataArray2[order]

mean1 = np.mean(dataArray1)
mean2 = np.mean(dataArray2)

diffs = np.abs(dataArray1 - dataArray2)

P_runaway, gapSize = runaway_probability(dataArray2, return_gap = True)
print(f"Runaway probability: {np.round(P_runaway,3)}")

# %% Load full area covers and b-values from file
allAreaRatios = np.loadtxt(dataPath_allInfo + 'allAreaRatios.txt', delimiter=" ")
allbValues = np.loadtxt(dataPath_allInfo + 'allbValues.txt', delimiter=" ")
allAsperityCounts = np.loadtxt(dataPath_allInfo + 'allAsperityCounts.txt', delimiter=" ")
            
                       
MLE_estimate_1 = MLE_b_value(dataArray1, m_c)
b1 = MLE_estimate_1[0]
a1 = MLE_estimate_1[1]

MLE_estimate_2 = MLE_b_value(dataArray2, m_c)
b2 = MLE_estimate_2[0]
a2 = MLE_estimate_2[1]

#powerlaw_estimate_1 = MLE_powerlaw_exponent(dataArray1, 10)

histo_mags1 = np.linspace(np.min(dataArray1), np.max(dataArray1), 100)
histo_mags2 = np.linspace(np.min(dataArray2), np.max(dataArray2), 100)
#histoBins = np.linspace(1, 3.5, 10) # for magnitudes
histoBins = np.linspace(1, 2, 10) #
histoLabels = ["eventsize [el. count]", "N"]
histoTitle = 'Eventsize Histogram'
infoString = f"S=0.02, cover 81%, P_r = {np.round(P_runaway,3)}"

# plotHistogram(dataArray1, histoBins, figLabels = histoLabels, figTitle = histoTitle)

ax0 = plotHistoCum([dataArray1, dataArray2], histoBins, figLabels = histoLabels, figTitle = 'Magnitude exceedance curve \n' + infoString)
ax0.plot(histo_mags1, GRlaw(histo_mags1, m_c, b1))
ax0.plot(histo_mags2, GRlaw(histo_mags2, m_c, b2))
ax0.legend([f"Predicted, b={np.round(b1, 3)}", f"Simulated, b={np.round(b2, 3)}"], prop={'size': 15})
#ax0.set_xscale('log')
plt.show()

# fig4, ax4 = plt.subplots()
# fig4.set_size_inches(10,6)
# ax4.scatter(np.arange(int(asperityCount*mapCount)), np.sort(dataArray1, axis = None))
# ax4.scatter(np.arange(int(asperityCount*mapCount)), np.sort(dataArray2, axis = None))
# ax4.legend(["Predicted", "Simulated"])
# ax4.set_title("Predicted vs. simulated magnitudes")
# ax4.set_xlabel("index")
# ax4.set_ylabel("magnitude")

fig6, ax6 = plt.subplots()
fig6.set_size_inches(10,6)
ax6.scatter(dataArray1_ordered, dataArray1_ordered)
ax6.scatter(dataArray1_ordered, dataArray2_ordered)
ax6.legend(["Predicted", "Simulated"], prop={'size': 15})
ax6.set_title("Magnitude comparison: " + infoString)
ax6.set_xlabel("predicted magnitude")
ax6.set_ylabel("magnitude")
plt.grid()

# %%
# fig5, ax5 = plt.subplots()
# fig5.set_size_inches(10,6)
# sc = ax5.scatter(allAreaRatios, allAsperityCounts, s = 110, c = allbValues, cmap = 'plasma')
# # ax5.axvline(0.05)
# # ax5.axvline(0.1)
# # ax5.axvline(0.15)
# ax5.set_ylim([90, 110])
# cbar = plt.colorbar(sc)
# cbar.set_label('b-value')
# ax5.set_title("100 maps: analysis")
# ax5.set_xlabel("area cover ratio")
# ax5.set_ylabel("asperity count")