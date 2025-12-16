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
        case "x":
            profile = data[:, profileInd]
        case "y":
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
    
def plotProfiles(data, labels, plotDim, plotInd, globalTitle='none',):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    for arr, linelabel in zip(data, labels):
        match plotDim:
            case "x":
                profile = arr[:,plotInd]
            case "y":
                profile = arr[plotInd, :]
            case _:
                raise ValueError("Unknown plotting direction.")
                
        ax.plot(profile, label=linelabel, linewidth=3)
    plt.grid()
    plt.legend()
    plt.show
        
    if globalTitle is not None:
        fig.suptitle(globalTitle, fontsize=16)