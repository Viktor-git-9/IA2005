#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 15:11:49 2025

@author: viktor
"""
import numpy as np
import matplotlib.pyplot as plt

def getProfiles(data, profileDim, profileInd):
    match profileDim:
        case "x": # In the Fortran data, the first index corresponds to x...
            profile = data[:, profileInd]
        case "y": # ... and the second to y.
            profile = data[profileInd, :]
        case _:
            raise ValueError("Unknown plotting direction.")
    return profile

def plotProfiles2x2(datalist, titles=None, globalTitle='none'):
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
    
    for ax, profile, title in zip(axes.ravel(), datalist, titles):
            
        im = ax.plot(profile)        
        ax.set_title(title)
        
    if globalTitle is not None:
        fig.suptitle(globalTitle, fontsize=16)
        
    fig.supxlabel("x [element counts]")
    fig.supylabel("y [element counts]")
    plt.tight_layout()
    plt.show()
    
def plotProfiles2x2plus1(datalist, titles=None, globalTitle='none'):
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
    
    for ax, profile, title in zip(axes.ravel(), datalist, titles):     
        im = ax.plot(profile)
        ax.set_title(title)
        
    extraax = fig.add_axes([1.05, 0.26, 0.4, 0.4])
    extraim = extraax.plot(datalist[-1])
    extraax.set_title(titles[-1])
        
    if globalTitle is not None:
        fig.suptitle(globalTitle, fontsize=16)
        
    fig.supxlabel("x [element counts]")
    fig.supylabel("y [element counts]")
    plt.tight_layout()
    plt.show()
    
def plotProfiles(data, axesLabels, lineLabels, globalTitle='none', stretchFactor = 'none'):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    
    #xvals:
    if stretchFactor is not None:
        N = len(data[0])
        xvals = np.linspace(1, stretchFactor * N, N)
    else:
        xvals = np.linspace(1, N, N)
    
    for profile, linelabel in zip(data, lineLabels):              
        ax.plot(xvals, profile, label=linelabel, linewidth=3)
    plt.grid()
    plt.legend(fontsize="large")
    plt.show
    
    ax.set_xlabel(axesLabels[0])
    ax.set_ylabel(axesLabels[1])
       
    if globalTitle is not None:
        fig.suptitle(globalTitle, fontsize=16)