#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 10:37:48 2026

@author: viktor
"""

import numpy as np
from pathlib import Path
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

def findBreakableArea(data, backgroundValue):
    totalArea      = data.size
    backgroundArea = np.sum(np.count_nonzero(data == backgroundValue))
    breakableArea  = totalArea - backgroundArea
    areaRatio = breakableArea/totalArea
    
    return breakableArea, areaRatio



directory = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/plot_heterogeneity/data/small_maps/smoother/varied_A")
shape1 = [256, 256]

breakableArea = []
areaRatio     = []
fileNames     = []

for file in sorted(directory.iterdir()):
    print(file)
    heteroMap = [loadBinArray(file, shape1)]
    currentBreakableArea, currentAreaRatio = findBreakableArea(np.array(heteroMap), 1000)
    
    fileNames.append(file)
    breakableArea.append(currentBreakableArea)
    areaRatio.append(currentAreaRatio)
    
#breakableArea.sort()
#areaRatio.sort()

#xData = np.array([3, 5, 7, 9, 11])
xData = np.array([1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0])
fig, ax1 = plt.subplots()

ax1.scatter(xData, areaRatio)
plt.grid()
plt.title("D = 2, N_0 = 64, r_0 = 5.65 m")
ax1.set_xlabel("A")
ax1.set_ylabel("R")


    