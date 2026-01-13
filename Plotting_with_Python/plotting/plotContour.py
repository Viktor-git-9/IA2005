#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 11:08:19 2025

@author: viktor
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def plotContours(data, titles=None, cbarLabels=None, cmap='inferno', globalTitle='none'):
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
        im = ax.imshow(arr, cmap=cmap, origin='lower')
        ax.set_aspect('equal')      # <- enforces square axes
        ax.set_title(title)
        cbar = plt.colorbar(im, ax=ax)
        
        if cbarLabels is not None:
            cbar.set_label(cbarLabel)

    # Global labels and title
    if globalTitle != 'none':
        fig.suptitle(globalTitle, fontsize=16)

    fig.supxlabel("x [element counts]")
    fig.supylabel("y [element counts]")

    plt.tight_layout()
    plt.show()

def plotContour2x2(data1, data2, data3, data4, titles=None, cmap='inferno', globalTitle='none'):
    """
    Plot four 2D arrays in a 2x2 grid with individual colorbars.

    Parameters
    ----------
    A, B, C, D : 2D numpy arrays
        Arrays of identical shape.
    titles : list of str, optional
        A list of four titles for the subplots.
    cmap : str
        Matplotlib colormap name (default 'viridis').
    """
    
    if titles is None:
        titles = ["A", "B", "C", "D"]
        
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    arrays = [data1, data2, data3, data4]
    
    for ax, arr, title in zip(axes.ravel(), arrays, titles):
        im = ax.imshow(arr, cmap=cmap, origin='lower', aspect='auto')
        
        ax.set_title(title)
        plt.colorbar(im, ax=ax)
        
    if globalTitle is not None:
        fig.suptitle(globalTitle, fontsize=16)
        
    fig.supxlabel("x [element counts]")
    fig.supylabel("y [element counts]")
    plt.tight_layout()
    plt.show()
    
def plotContour2x2plus1(data1, data2, data3, data4, data5=None, titles=None, cmap='inferno', globalTitle='none'):
    """
    Plot four 2D arrays in a 2x2 grid with individual colorbars.

    Parameters
    ----------
    A, B, C, D : 2D numpy arrays
        Arrays of identical shape.
    titles : list of str, optional
        A list of four titles for the subplots.
    cmap : str
        Matplotlib colormap name (default 'viridis').
    """
    
    if titles is None:
        titles = ["A", "B", "C", "D", "E"]
        
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    arrays = [data1, data2, data3, data4]
    
    for ax, arr, title in zip(axes.ravel(), arrays, titles):
        im = ax.imshow(arr, cmap=cmap, origin='lower', aspect='auto')
        
        ax.set_title(title)
        plt.colorbar(im, ax=ax)
    
    if data5 is not None:
        extraax = fig.add_axes([1.05, 0.26, 0.4, 0.4])
        extraim = extraax.imshow(data5, cmap=cmap, origin='lower', aspect='auto')
        extraax.set_title(titles[-1])
        plt.colorbar(extraim, ax=extraax)
        
    if globalTitle is not None:
        fig.suptitle(globalTitle, fontsize=16)
        
    fig.supxlabel("x [element counts]")
    fig.supylabel("y [element counts]")
    plt.tight_layout()
    plt.show() 