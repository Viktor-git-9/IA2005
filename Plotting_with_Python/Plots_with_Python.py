# -*- coding: utf-8 -*-
"""
Spyder-Editor

Dies ist eine temporäre Skriptdatei.
"""
import numpy as np
import matplotlib.pyplot as plt

nx = 64
ny = 64
nt = 501
timeIndex = 70

ruptureTimes = np.fromfile("ruptureTimes3.bin", dtype=np.float64) \
      .reshape((nx, ny, nt), order='F')
print(np.min(ruptureTimes[:, :, timeIndex]))

rtSlice = ruptureTimes[:, :, timeIndex]

fig=plt.figure()
plt.contourf(ruptureTimes[:, :, timeIndex])
plt.colorbar()
