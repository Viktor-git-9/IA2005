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
    
def plotProfiles(data, axesLabels, lineLabels, globalTitle=None, stretchFactor=None):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    
    #xvals:
    N = len(data[0])
    if stretchFactor is not None:
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

def plot2Profiles(arrays, xlabel_bottom="Index (Array 1)", xlabel_top="Index (Array 2)",
                    ylabel="Value", title=None, labels=None, lineInd=None):
    """
    Plot two 1D numpy arrays with different lengths using two x-axes.

    Parameters
    ----------
    arrays : list of numpy.ndarray
        List containing exactly two 1D arrays.
    xlabel_bottom : str
        Label for the bottom x-axis (array 1).
    xlabel_top : str
        Label for the top x-axis (array 2).
    ylabel : str
        Label for the y-axis.
    title : str
        Plot title.
    labels : list of str
        Labels for the two lines (for legend).
    """

    if len(arrays) != 2:
        raise ValueError("arrays must contain exactly two numpy arrays")

    arr1, arr2 = arrays

    if labels is None:
        labels = ["Array 1", "Array 2"]

    x1 = np.arange(len(arr1))
    x2 = np.arange(len(arr2))

    fig, ax1 = plt.subplots()

    # Bottom axis (array 1)
    line1, = ax1.plot(x1, arr1, color="blue", label=labels[0])
    ax1.set_xlabel(xlabel_bottom)
    ax1.set_ylabel(ylabel)

    # Top axis (array 2)
    ax2 = ax1.twiny()
    line2, = ax2.plot(x2, arr2, color="red", label=labels[1])
    ax2.set_xlabel(xlabel_top)
    
    if lineInd:
        plt.axvline(lineInd)

    if title:
        ax1.set_title(title)

    # Combined legend
    lines = [line1, line2]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper left")

    plt.tight_layout()
    plt.show()