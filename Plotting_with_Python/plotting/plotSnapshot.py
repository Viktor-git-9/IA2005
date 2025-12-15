#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 11:08:19 2025

@author: viktor
"""
import numpy as np
import matplotlib.pyplot as plt

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
    

    
def plotProfiles2x2(data1, data2, data3, data4, plotDim, plotInd, titles=None, globalTitle='none'):
    """
    Description goes here

    Parameters
    ----------
    data1 : TYPE
        DESCRIPTION.
    data2 : TYPE
        DESCRIPTION.
    data3 : TYPE
        DESCRIPTION.
    data4 : TYPE
        DESCRIPTION.
    titles : TYPE, optional
        DESCRIPTION. The default is None.
    globalTitle : TYPE, optional
        DESCRIPTION. The default is 'none'.

    Returns
    -------
    None.

    """
    if titles is None:
        titles = ["A", "B", "C", "D"]
      
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    arrays = [data1, data2, data3, data4]
    
    for ax, arr, title in zip(axes.ravel(), arrays, titles):
        match plotDim:
            case "x":
                profile = arr[:, plotInd]
            case "y":
                profile = arr[plotInd,:]
            case _:
                raise ValueError("Unknown plotting direction.")
            
        im = ax.plot(profile)
        
        ax.set_title(title)
        
    if globalTitle is not None:
        fig.suptitle(globalTitle, fontsize=16)
        
    fig.supxlabel("x [element counts]")
    fig.supylabel("y [element counts]")
    plt.tight_layout()
    plt.show()
    
def plotProfiles(data, labels, title, plotDim, plotInd):
    print("help")
    #something goes here
    
    
    