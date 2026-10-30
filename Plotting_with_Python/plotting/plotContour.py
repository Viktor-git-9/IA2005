#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 11:08:19 2025

@author: viktor
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Annulus
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plotContours(data, titles=None, cbarLabels=None, clims = None, cmap='inferno', globalTitle='none'):
    """
    Plot a list of 2D arrays in a 2 x (N/2) grid.
    If N is odd, add one extra plot to the right.

    Parameters
    ----------
    data : list of 2D numpy arrays
    titles : list of str, optional
    cmap : str
    globalTitle : str
    """

    N = len(data)
    if N == 0:
        raise ValueError("Data list must not be empty.")

    if titles is None:
        titles = [f"Array {i+1}" for i in range(N)]

    if len(titles) != N:
        raise ValueError("Length of titles must match length of data.")
        
    if cbarLabels is None:
        cbarLabels = [None] * N
    if len(cbarLabels) != N:
        raise ValueError("Length of cbar_labels must match length of data.")

    ncols = N // 2          # number of columns in main grid
    odd = (N % 2 == 1)
    
    total_cols = ncols + (1 if odd else 0)
    
    fig = plt.figure(figsize=(5*total_cols, 8))
    gs = GridSpec(2, total_cols, figure=fig)
    
    axes = []
    
    # Create main grid
    for i in range(2):
        for j in range(ncols):
            ax = fig.add_subplot(gs[i, j])
            axes.append(ax)
            
    # extra subplot on the right if necessary
    if odd:
        extra_ax = fig.add_subplot(gs[:, -1])
        axes.append(extra_ax)
        
    # plot plot plot 
    for ax, arr, title, cbarLabel in zip(axes, data, titles, cbarLabels):
        if clims is not None:
            im = ax.imshow(arr.T, cmap=cmap, vmin=clims[0], vmax=clims[1], origin='lower')
        else:
            im = ax.imshow(arr.T, cmap=cmap, origin='lower')
        ax.set_aspect('equal')      # <- enforces square axes
        ax.set_title(title)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax, format='%.0e')
        
        annulus = Annulus((170.79, 156.34), r=48, width=2, color='c', linestyle = '-.')
        #ax.add_patch(annulus)
        
        if cbarLabels is not None:
            cbar.set_label(cbarLabel)

    # Global labels and title
    if globalTitle != 'none':
        fig.suptitle(globalTitle, fontsize=16)

    fig.supxlabel("x [element counts]")
    fig.supylabel("y [element counts]")

    plt.tight_layout()
    #plt.show()
    
    return fig


def plotSingleContour(arr, title=None, cbarLabel=None, clims=None, cmap='inferno', globalTitle='none'):
    """
    Plot a single 2D array as an image/contour plot with a colorbar.

    Parameters
    ----------
    arr : 2D numpy array
    title : str, optional
    cbarLabel : str, optional
    clims : tuple (vmin, vmax), optional
    cmap : str
    globalTitle : str
    """

    if arr.ndim != 2:
        raise ValueError("Input array must be 2D.")

    if title is None:
        title = "Array"

    fig, ax = plt.subplots(figsize=(6, 5))

    # plot
    if clims is not None:
        im = ax.imshow(arr.T, cmap=cmap, vmin=clims[0], vmax=clims[1], origin='lower')
    else:
        im = ax.imshow(arr.T, cmap=cmap, origin='lower')

    ax.set_aspect('equal')
    ax.set_title(title)
    ax.set_xlabel("x [element counts]")
    ax.set_ylabel("y [element counts]")

    # colorbar matched to the axes height
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, format='%.0e')
    if cbarLabel is not None:
        cbar.set_label(cbarLabel)

    # global title
    if globalTitle != 'none':
        fig.suptitle(globalTitle, fontsize=16)

    fig.tight_layout()

    return fig
    
