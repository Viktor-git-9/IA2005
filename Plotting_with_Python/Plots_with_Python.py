# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

nx = 64
ny = 64
nt = 501
timeIndex = 70

ruptureTimes = np.fromfile("ruptureTimes1.bin", dtype=np.float64) \
      .reshape((nx, ny, nt), order='F')
slipHistories = np.fromfile("slipHistories1.bin", dtype=np.float64) \
      .reshape((nx, ny, nt), order='F')

rtSlice = ruptureTimes[:, :, timeIndex]

fig=plt.figure(1)
plt.contourf(ruptureTimes[:, :, timeIndex])
plt.title("Rupture time")
plt.colorbar()

fig=plt.figure(2)
plt.contourf(ruptureTimes[:, :, timeIndex])
plt.title("Slip")
plt.colorbar()

