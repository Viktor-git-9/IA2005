#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:32:52 2026

@author: viktor
"""
# %%

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.sampling as sampling
import scipy.integrate as integrate
import numpy.random as rand
from plotting.plotContour import plotContours # import plotting functions

# %%
def createAsperity(DcMap, asperityDc, asperityRadius, asperityLocation):  
    for i in range(int(asperityLocation[0]-asperityRadius)-1, int(asperityLocation[0]+asperityRadius)+1, 1):
        for j in range(int(asperityLocation[1]-asperityRadius)-1, int(asperityLocation[1]+asperityRadius)+1, 1):
            rad = np.sqrt((i-asperityLocation[0])**2 + (j-asperityLocation[1])**2)
            if rad <= asperityRadius:
                if DcMap[i,j] > asperityDc:
                    DcMap[i,j] = asperityDc
                    
    return DcMap
    

def DcMapGenerator(ixmax, DcVals, radii):    
    DcMap = np.zeros((ixmax, ixmax)) + 1000
    
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
                
    return DcMap
        

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

def findBreakableArea(data, backgroundValue):
    totalArea      = data.size
    backgroundArea = np.sum(np.count_nonzero(data == backgroundValue))
    breakableArea  = totalArea - backgroundArea
    areaRatio = breakableArea/totalArea
    
    return breakableArea, areaRatio

# %%
a = 1.5
b = 1
delta_sigma = 3*1e6
K = 10**(a + 2*b/3*(9.1 - np.log10(16/7*delta_sigma)))
#K = 3000
rmin = K**(1/(2*b))
#rmin = 4

GRlaw = lambda M: 10**a * 10**(-b*M)
pdf = lambda r: 2*b*K * r**(-2*b-1)
ccdf = lambda r: K * r**(-2*b)
cdf_inverse = lambda r: (K / (1-r) )**(1/(2*b))
moment2mag  = lambda m: (np.log10(m) - 9.1) / 1.5

r = np.linspace(rmin, 1000, 10000)

pdf_vals = pdf(r)
cdf_vals = 1 - ccdf(r)

# %%

mapCount = 10
areaRatios = []
bVals = []

for iMap in range(0, mapCount):
    
    u = np.random.rand(int(10**a))      # Uniform(0,1)
    r_samples = cdf_inverse(u)

    artificial_moments = 16/7 * delta_sigma * r_samples**3
    artificial_magnitudes = moment2mag(artificial_moments)
    artificial_magnitudes_sorted = np.sort(artificial_magnitudes)[::-1]
    logOfCounts = np.log10(np.arange(1, len(artificial_magnitudes_sorted)+1))
    artificial_b, artificial_lna = np.polyfit(artificial_magnitudes_sorted, logOfCounts, 1)
    bVals.append(artificial_b)
    
    radii = r_samples/4
    radii[radii>100] = 100 # artificially get rid of too large asperities. Gotta fix this by introducing a rmax!
    ixmax = 256
    DcVals = 0.1*radii
    
    DcMap = DcMapGenerator(ixmax, DcVals, radii)
    breakableArea, areaRatio = findBreakableArea(np.array(DcMap), 1000)
    
    areaRatios.append(areaRatio)

artificialbMean = np.mean(bVals)
artificialbVar = np.var(bVals)
coverageMean = np.mean(areaRatios)
coverageVar  = np.var(areaRatios)

print(f"Coverage ratio average: {np.round(coverageMean, 3)}")
print(f"Coverage ratio variance: {np.round(coverageVar, 3)}")

# %%
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

#fig3, ax3 = plt.subplots()
#ax3.scatter(np.arange(int(10**a)), r_samples/4)
#ax3.set_title("Radii")

fig4, ax4 = plt.subplots()
ax4.scatter(np.arange(int(10**a)), artificial_magnitudes)
ax4.set_title("Magnitudes")

fig5, ax5 = plt.subplots()
ax5.scatter(artificial_magnitudes_sorted, logOfCounts)
ax5.set_title("Scatter plot of magnitude exceedance")

ax0 = plotHistoCum(artificial_magnitudes, histoBins, figLabels = histoLabels, figTitle = 'Test: Magnitude exceedance curve \n from radii')
ax0.plot(histo_mags, GRlaw(histo_mags))
plt.show()

fig10 = plotContours([DcMap], clims=[0.25, 10], cbarLabels=["Dc [mm]"], titles=[f"Asperity count = {len(radii)}, Area ratio = {round(areaRatio,3)}"], globalTitle="Dc Map")
