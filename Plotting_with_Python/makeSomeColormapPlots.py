#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 11:23:39 2026

@author: viktor
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

def plotMap(x, y, Z, axesLabels=None, colorbarLabel=None, title=None,
            showPoints=True, resolution=300, method="cubic"):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    Z = np.asarray(Z, dtype=float)      # shape (len(y), len(x))

    # Interpolator on the original, possibly uneven grid.
    # Z is indexed [y, x], so the axes are passed as (y, x).
    interp = RegularGridInterpolator((y, x), Z, method=method)

    # Fine, evenly spaced grid covering the same range
    xf = np.linspace(x[0], x[-1], resolution)
    yf = np.linspace(y[0], y[-1], resolution)
    XF, YF = np.meshgrid(xf, yf)
    Zf = interp(np.stack([YF.ravel(), XF.ravel()], axis=-1)).reshape(XF.shape)

    fig, ax = plt.subplots()
    im = ax.imshow(
        Zf,
        cmap="plasma",
        interpolation="nearest",    # already smooth, no need for further interpolation
        origin="lower",
        extent=[x[0], x[-1], y[0], y[-1]],
        aspect="auto",
        vmin=Z.min(), vmax=Z.max(), # keep the color range tied to the real data
    )

    if showPoints:
        X, Y = np.meshgrid(x, y)
        ax.scatter(X.ravel(), Y.ravel(), s=20, facecolor="white",
                   edgecolor="white", linewidth=1, zorder=3)

    if axesLabels is not None:
        ax.set_xlabel(axesLabels[0])
        ax.set_ylabel(axesLabels[1])
    ax.set_title(title)
    fig.colorbar(im, label=colorbarLabel)
    plt.show()

# Manually defined square array (3x3 here, but any N x N works)
Z1 = np.array([
    [0.0, 0.0, 0.09, 0.474, 0.836, 0.941],
    [0.0, 0.0, 0.47, 0.753, 0.894, 0.98],
    [0.0, 0.327, 0.75, 0.897, 0.979, 0.99],
    [0.323, 0.865, 0.95, 0.959, 1.0, 1.0],
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
])


NaN = 0
Z2 = np.array([
    [NaN, NaN, NaN, 0.396, 0.537, 0.696],
    [NaN, NaN, 0.30, 0.456, 0.406, 0.687],
    [NaN, 0.165, 0.41, 0.533, 0.589, 0.677],
    [0.123, 0.466, 0.54, 0.546, NaN, NaN],
    [NaN, NaN, NaN, NaN, NaN, NaN]
])

# Manually defined axes (one value per column / row)
x = np.array([0.15, 0.26, 0.40, 0.49, 0.69, 0.81]) # cover ratio
y = np.array([0.67, 0.43, 0.25, 0.11, 0.02]) # S-value

plotMap(x, y, Z1, ["cover ratio", "S-value"], "P(runaway)", "Colormap of runaway probability")
plotMap(x, y, Z2, ["cover ratio", "S-value"], "Relative gap", "Colormap of relative gap size")

