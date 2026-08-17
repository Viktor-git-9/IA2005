#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 14:20:43 2026

@author: viktor
"""

import numpy as np
import matplotlib.pyplot as plt
from plotting.plotContour import plotContours # import plotting functions

def loadBinArray(filename, shape, dtype=np.float64):
    """
    Load a Fortran-written binary stream file into a numpy array.

    Arguments
    ----------
    filename : str
        Path to the .bin file.
    timeIndex: int
        Stage-specific time step at which data is to be loaded
    shape : tuple of ints
        Desired shape of the array
    dtype : numpy dtype, optional
        Data type. Default is float64 (Fortran double precision).

    Returns
    -------
    numpy array with given shape, Fortran memory order.
    """
    arr = np.fromfile(filename, dtype=dtype) \
                .reshape(shape, order='F')
    arr = np.array(arr)
    return arr

def loadTXTarray(filename):
    data = np.loadtxt(filename)

    # sort by x first, then y
    data_sorted = data[np.lexsort((data[:,1], data[:,0]))]
    
    x = np.unique(data_sorted[:,0])
    y = np.unique(data_sorted[:,1])
    
    nx = len(x)
    ny = len(y)
    
    V = data_sorted[:,2].reshape(nx, ny)
    return V

def findBreakableArea(data, backgroundValue):
    totalArea      = data.size
    backgroundArea = np.sum(np.count_nonzero(data == backgroundValue))
    breakableArea  = totalArea - backgroundArea
    areaRatio = breakableArea/totalArea
    
    return breakableArea, areaRatio

#datapath = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/plot_heterogeneity/data/small_maps/smoother/"
datapath = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/single_asperities/r_148/hetero/"
#filename2 = "full_hetero.bin"
filename1 = "hetero.bin"
filename3 = "X0.bin"
filename4 = "Y0.bin"
#shape1 = [4096, 4096]
shape1 = [256, 256]
shape2 = [4096, 4096]
shape3 = 1

showStop = 255
#scatterInds = [20, 55, 58, 97, 80, 23, 78, 115, 52]
#scatterInds = [45, 43, 29, 17, 38, 5, 7, 24, 11]
#scatterInds = np.array([18, 27, 8, 23, 52, 6, 13, 44, 53])-1 # 4_4
#scatterInds = np.array([62, 13, 69, 15, 61, 23, 2, 33, 95])-1 # 6_4
#scatterInds = np.array([14, 18, 16, 2, 6, 27, 5, 9, 23])-1 # 2_4
scatterInds = np.array(range(0, 100, 1))

#hetero_renorm_off = [loadBinArray(datapath + filename1, shape1)[96:160,96:160]]
hetero_renorm_off = [loadBinArray(datapath + filename1, shape1)[0:showStop, 0:showStop]]
#full_hetero  = [loadBinArray(datapath + filename2, shape2)[0:showStop, 0:showStop]]
x0  = [loadBinArray(datapath + filename3, shape3)]
#x0[0] = x0[0] - 1920 # shift to center of full DcMap
y0  = [loadBinArray(datapath + filename4, shape3)]
#y0[0] = y0[0] - 1920

breakableArea, areaRatio = findBreakableArea(np.array(hetero_renorm_off), 1000)

#hetero_renorm_off[0] = np.transpose(hetero_renorm_off[0])
fig1 = plotContours(hetero_renorm_off, clims=[0.25, 3], cbarLabels=["Dc [mm]"], titles=["n_hypo = 100, Species = 4"], globalTitle="Dc Map")
ax1 = fig1.axes[0]
for i in scatterInds:
    #print(i)
    ax1.scatter(x0[0][i], y0[0][i], color='red', s=5, zorder=100)
    ax1.text(x0[0][i]*(1+0.01), y0[0][i]*(1+0.01), str(i+1), color='red')
plt.show()

# fig2 = plotContours(full_hetero, clims=[0,4], cbarLabels=["Dc [mm]"], titles=["Dc heterogeneity"], globalTitle="Renormalization ON")
# ax2 = fig2.axes[0]
# #ax2.scatter(shape2[0]/2, shape2[1]/2, color='red', s=20, zorder=100)
# plt.show()

print(np.min(x0))
print(np.max(x0))

