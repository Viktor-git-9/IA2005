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

#datapath = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/Kram/plot_heterogeneity/data/fullsize_with_different_densities/"
#datapath = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/Kram/plot_heterogeneity/data/4_4_DcMaps_from_center/"
#datapath = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/Kram/plot_heterogeneity/data/256/cut_from_center/"
datapath = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/plot_heterogeneity/data/4_7_256/heteros/"
filename1 = "heterogeneity0_4.bin"
#filename1 = "hetero_4_4.bin"
# filename2 = "x0_4_4.bin"
# filename3 = "y0_4_4.bin"
shape1 = [256, 256]
#shape1 = [64, 64]
shape2 = 64

showStop = 255

hetero1 = [loadBinArray(datapath + filename1, shape1)[0:showStop, 0:showStop]]
# x0  = [loadBinArray(datapath + filename2, shape2)]
# x0[0] = x0[0] - 1920 # shift to center of full DcMap
# y0  = [loadBinArray(datapath + filename3, shape2)]
# y0[0] = y0[0] - 1920

fig1 = plotContours(hetero1, clims=[0,4], cbarLabels=["Dc [mm]"], titles=["Density = 4, Scales = 4, n_hypo = 64"], globalTitle="Dc Map")
#plt.scatter(x0[0], y0[0], color='red', s=5, zorder=100)
ax1 = fig1.axes[0]
plt.show()

