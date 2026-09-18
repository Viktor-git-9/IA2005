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
            if file.endswith('.txt'):
                dataList.append(np.loadtxt(dataPath + file, delimiter=" "))
                
    return dataList
            
# ------------------------------------------------------------------
# MLE estimation & GR law35712339244214043571233924421404
# ------------------------------------------------------------------
def MLE_b_value(magnitudes, cutoffMagnitude):
    MLE_b = np.log10(np.exp(1))/(np.mean(magnitudes) - cutoffMagnitude)
    MLE_b = MLE_b * (len(magnitudes) - 1)/len(magnitudes) # small sample bias correction
    MLE_a = np.log10(len(magnitudes)) + MLE_b * cutoffMagnitude # get a from MLE estimation for b
    
    return MLE_b, MLE_a

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

#dataPath1 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/9Maps_5/asperityMagnitudes/"
#dataPath2 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/9Maps_5/simulationMagnitudes/"

dataPath1 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/3_bins_5_15_25_35_50/infofiles/above10percent/"
dataPath2 = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/3_bins_5_15_25_35_50/datafiles/above10percent/"
dataPath_allInfo = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/3_bins_5_15_25_35_50/infofiles/"
m_c = 1.0928
correctionTerm = 0.12
#asperityCount = 5
mapCount = 10

# %% Load and treat predicted & simulated magnitudes
dataList1 = getMagnitudes(dataPath1, mapCount)
dataList2 = getMagnitudes(dataPath2, mapCount)
    
dataArray1 = np.concatenate(dataList1)
dataArray2 = np.concatenate(dataList2)

dataArray2 = dataArray2 - correctionTerm

mean1 = np.mean(dataArray1)
mean2 = np.mean(dataArray2)

diffs = np.abs(dataArray1 - dataArray2)

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

GRlaw = lambda m: 10**(-b1*(m-m_c)) # Gutenberg-Richter law

histo_mags = np.linspace(np.min(dataArray1), np.max(dataArray2), 100)
histoBins = np.linspace(1, 3.5, 10)
histoLabels = ["magnitude", "N"]
histoTitle = 'Magnitude Histogram'
plotHistogram(dataArray1, histoBins, figLabels = histoLabels, figTitle = histoTitle)

ax0 = plotHistoCum([dataArray1, dataArray2], histoBins, figLabels = histoLabels, figTitle = 'Magnitude exceedance curve: sparse maps')
#ax0.plot(histo_mags, GRlaw(histo_mags))
ax0.legend([f"Predicted, b={np.round(b1, 3)}", f"Simulated, b={np.round(b2, 3)}"])
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
ax6.scatter(np.sort(dataArray1, axis = None), np.sort(dataArray1, axis = None))
ax6.scatter(np.sort(dataArray1, axis = None), np.sort(dataArray2, axis = None))
ax6.legend(["Predicted", "Simulated"])
ax6.set_title("Predicted vs. simulated magnitudes")
ax6.set_xlabel("predicted magnitude")
ax6.set_ylabel("magnitude")
plt.grid()

# %%
fig5, ax5 = plt.subplots()
fig5.set_size_inches(10,6)
sc = ax5.scatter(allAreaRatios, allAsperityCounts, s = 110, c = allbValues, cmap = 'plasma')
ax5.axvline(0.05)
ax5.axvline(0.1)
ax5.axvline(0.15)
ax5.set_ylim([0, 50])
cbar = plt.colorbar(sc)
cbar.set_label('b-value')
ax5.set_xlabel("area cover ratio")
ax5.set_ylabel("asperity count")