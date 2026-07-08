#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:32:52 2026

@author: viktor
"""
# %% Import what's needed

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
import pickle
from plotting.plotContour import plotContours # import plotting functions

# %% Define functions to build asperity map and make a cumulative plot

# create an asperity on a grid with specified location, radius, and value
def createAsperity(DcMap, asperityDc, asperityRadius, asperityLocation):  
    for i in range(int(asperityLocation[0]-asperityRadius)-1, int(asperityLocation[0]+asperityRadius)+1, 1):
        for j in range(int(asperityLocation[1]-asperityRadius)-1, int(asperityLocation[1]+asperityRadius)+1, 1):
            rad = np.sqrt((i-asperityLocation[0])**2 + (j-asperityLocation[1])**2)
            if rad <= asperityRadius:
                if DcMap[i,j] > asperityDc:
                    DcMap[i,j] = asperityDc
                    
    return DcMap
    

# create an asperity map by placing many asperities on a grid. No cuts at the boundaries allowed
def DcMapGenerator(ixmax, DcVals, radii):    
    DcMap = np.zeros((ixmax, ixmax)) + 1000
    X0    = []
    Y0    = []
    radii = np.sort(radii)[::-1]
    DcVals = np.sort(DcVals)[::-1]
    
    for irad, rad in enumerate(radii):
        asperityDc = DcVals[irad]
        check = 0
        skipCount = 0
        while check == 0:
            proposedLocation = rand.uniform(0, ixmax, 2)
            if proposedLocation[0] - rad < 1 or proposedLocation[0] + rad > ixmax or proposedLocation[1] - rad < 1 or proposedLocation[1] + rad > ixmax:
                skipCount = skipCount + 1
            else:
                location = proposedLocation
                DcMap = createAsperity(DcMap, asperityDc, rad, location)
                check = 1
                X0.append(proposedLocation[0])
                Y0.append(proposedLocation[1])
            
    X0 = np.array(X0)
    Y0 = np.array(Y0)
    return DcMap, X0, Y0
        

# make a cumulative loglog plot of some data
def plotHistoCum(data, bins, figLabels = None, figTitle = None, figSubTitle = None, textCount = None):
    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(1,1,1)
    ax.ecdf(data, complementary=True)
    ax.set_yscale("log", base=10)
    
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
    return ax

# find total area covered by asperities
def findBreakableArea(data, backgroundValue):
    totalArea      = data.size
    backgroundArea = np.sum(np.count_nonzero(data == backgroundValue))
    breakableArea  = totalArea - backgroundArea
    areaRatio = breakableArea/totalArea
    
    return breakableArea, areaRatio

# %% Setup constants and functions

asperityCount = 100 # how many asperities are to be placed
b = 1.5 # b-value of underlying magnitude distribution. To be recovered by the samples!
r_c = 20 # minimums asperity radius in meter. Divide by 4 to get element count
delta_sigma = 3*1e6 # stress drop, for conversion from radius to magnitude and back
doSave = True

moment2mag  = lambda M: (np.log10(M) - 9.1) / 1.5 # https://gfzpublic.gfz.de/pubman/faces/ViewItemOverviewPage.jsp?itemId=item_4015
rad2moment  = lambda r: 16/7 * delta_sigma * r**3 # Madariaga&Ruiz (2016)

m_c = moment2mag(rad2moment(r_c)) # minimum (completeness) magnitude
a = np.log10(asperityCount)+b*m_c # a-value of prescribed magnitude distribution

GRlaw = lambda m: 10**(-b*(m-m_c)) # Gutenberg-Richter law
pdf = lambda r: 2*b/r_c * (r/r_c)**(-2*b-1) # pdf of GR law
ccdf = lambda r: (r/r_c)**(-2*b) # complementary cdf of GR law
cdf_inverse = lambda r: r_c * (1-r)**(-1/(2*b)) # inverse cdf of GR law. Use this for sampling


r = np.linspace(0.1, 1000, 10000) # some radii for distribution plotting
magnitudes = moment2mag(rad2moment(r)) # magnitudes from these radii

#pdf_vals = pdf(r)
#cdf_vals = 1 - ccdf(r)

# %% Draw samples, make maps

mapCount = 1 # how many samples / maps
areaRatios = [] # prepare some lists
bVals = []
aVals = []

# loop over all samples
for iMap in range(0, mapCount):
    
    u = np.random.rand(int(asperityCount))      # inverse transform sampling w/ cdf^-1
    r_samples = cdf_inverse(u)
    
    radii = r_samples # get rid of asperities that will be too large for the map. potential cause for bias
    correctedAsperities = len(radii[radii>400])
    if correctedAsperities != 0:
        print(f"Map {iMap}: Correcting {correctedAsperities} large asperities with r>100el." )
        radii[radii>400] = 400 # artificially get rid of too large asperities. Fix this by introducing a rmax?

    artificial_moments = rad2moment(radii) # Get moments and magnitudes from radius sample
    artificial_magnitudes = moment2mag(artificial_moments)
    artificial_magnitudes_sorted = np.sort(artificial_magnitudes)[::-1] # sort magnitudes, for plotting
    logOfCounts = np.log10(np.arange(1, len(artificial_magnitudes_sorted)+1))
    #artificial_b, artificial_a = np.polyfit(artificial_magnitudes_sorted, logOfCounts, 1)
    
    MLE_artificial_b = np.log10(np.exp(1))/(np.mean(artificial_magnitudes) - m_c) # do MLE estimation on magnitude sample to recover b
    MLE_artificial_b = MLE_artificial_b * (asperityCount - 1)/asperityCount
    MLE_artificial_a = np.log10(asperityCount) + MLE_artificial_b * m_c # get a from MLE estimation for b
    bVals.append(MLE_artificial_b)
    aVals.append(MLE_artificial_a)
        
    ixmax = 256
    DcVals = 0.05*radii/4 # calculate values for asperities, linear dependance
    
    DcMap, X0, Y0 = DcMapGenerator(ixmax, DcVals, radii/4) # build asperity map
    if doSave is True:
        file1 = open("hetero.bin", mode="wb")
        np.save(file1, DcMap)
        
        file2 = open("X0.bin", mode="wb")
        np.save(file2, X0)
        
        file3 = open("Y0.bin", mode="wb")
        np.save(file3, Y0)

    breakableArea, areaRatio = findBreakableArea(np.array(DcMap), 1000) # how much area do the asperities cover?
    areaRatios.append(areaRatio)

artificialbMean = np.mean(bVals) # calculate means of a, b, coverage over all samples
artificialbVar = np.var(bVals)
artificialaMean = np.mean(aVals)
artificialaVar = np.var(aVals)
coverageMean = np.mean(areaRatios)
coverageVar  = np.var(areaRatios)

print(f"Coverage ratio average: {np.round(coverageMean, 3)}")
print(f"Coverage ratio variance: {np.round(coverageVar, 3)}")

# %% Plots

histoBins = np.linspace(1, 3.5, 10)
histo_mags = np.linspace(np.min(artificial_magnitudes), np.max(artificial_magnitudes), 100)
histoLabels = ["magnitude", "N"]
histoTitle = 'Magnitude Histogram'

#fig1, ax1 = plt.subplots()
#ax1.plot(r, pdf_vals)
#ax1.set_title("PDF")

#fig2, ax2 = plt.subplots()
#ax2.plot(r, cdf_vals)
#ax2.set_title("CDF")

fig3, ax3 = plt.subplots()
ax3.scatter(np.arange(int(asperityCount)), radii/4)
ax3.set_title("Radii after correction")
ax3.set_xlabel("index")
ax3.set_ylabel("radii [el.]")

fig4, ax4 = plt.subplots()
ax4.scatter(np.arange(int(asperityCount)), artificial_magnitudes)
ax4.set_title("Artificial Magnitudes")
ax4.set_xlabel("index")
ax4.set_ylabel("magnitude")

fig5, ax5 = plt.subplots()
ax5.scatter(artificial_magnitudes_sorted, logOfCounts)
ax5.plot(magnitudes, a-magnitudes*b)
ax5.plot(magnitudes, MLE_artificial_a - magnitudes * MLE_artificial_b)
plt.xlim([0,4])
plt.ylim([-1, 3])
plt.legend(["Sampled magnitudes", "Prescribed distribution", "Fit from sample"])
ax5.set_title("Fit test for last generated Dc map")
ax5.set_xlabel("magnitude")
ax5.set_ylabel("log(count)")

fig6, ax6 = plt.subplots()
ax6.plot(magnitudes, a-magnitudes*b)
ax6.plot(magnitudes, artificialaMean - magnitudes*artificialbMean)
plt.xlim([0,4])
plt.ylim([-1, 3])
plt.legend(["Prescribed distribution", "Averaged fit from all samples"])
ax6.set_title("Fit test for all averaged samples")
ax6.set_xlabel("magnitude")
ax6.set_ylabel("log(count)")

ax0 = plotHistoCum(artificial_magnitudes, histoBins, figLabels = histoLabels, figTitle = 'Test: Magnitude exceedance curve \n from radii')
ax0.plot(histo_mags, GRlaw(histo_mags))
plt.show()

fig10 = plotContours([DcMap], clims=[0.25, 5.5], cbarLabels=["Dc [mm]"], titles=[f"Asperity count = {len(radii)}, Area ratio = {round(areaRatio,3)}"], globalTitle="Dc Map")
