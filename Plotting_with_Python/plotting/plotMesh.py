#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 11:57:33 2026

@author: viktor
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

def plotMeshSlices(
    data_list,
    slice_dim="x",        # "x" or "y"
    slice_index=0,        # index of the fixed spatial dimension
    spatial_coords=None, # x or y coordinates
    t_coords=None,       # time coordinates
    titles=None,
    xlabel="Time",
    ylabel=None,
    suptitle=None,
    elev=35,
    azim=60
):
    """
    Plot stacked 3D mesh plots (Matlab-like) from 3D arrays.

    Each array must have shape (nx, ny, nt).
    """

    nplots = len(data_list)
    fig = plt.figure(figsize=(20, 4 * nplots))

    for i, data in enumerate(data_list):

        ax = fig.add_subplot(nplots, 1, i + 1, projection="3d")

        if slice_dim == "x":
            slice_data = data[:, slice_index, :]   # (nx, nt)
            spatial = spatial_coords if spatial_coords is not None else np.arange(data.shape[0])
            spatial_label = "x"

        elif slice_dim == "y":
            slice_data = data[slice_index, :, :]   # (ny, nt)
            spatial = spatial_coords if spatial_coords is not None else np.arange(data.shape[1])
            spatial_label = "y"

        else:
            raise ValueError("slice_dim must be 'x' or 'y'")

        time = t_coords if t_coords is not None else np.arange(slice_data.shape[1])

        T, S = np.meshgrid(time, spatial)

        ax.plot_wireframe(
            T, S, slice_data,
            rstride=1,
            cstride=1,
            linewidth=0.2,
            color="0.3"
        )

        ax.view_init(elev=elev, azim=azim)

        # --- Individual axis labels ---
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel if ylabel is not None else spatial_label)
        ax.set_zlabel("")

        if titles is not None:
            ax.set_title(titles[i])

    if suptitle is not None:
        fig.suptitle(suptitle, y=0.98)

    plt.tight_layout()
    plt.show()