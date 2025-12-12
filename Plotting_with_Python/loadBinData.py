#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 13:07:23 2025

@author: viktor
"""

import numpy as np
import matplotlib.pyplot as plt

def loadBinArray(filename, timeIndex, shape, dtype=np.float64):
    """
    Load a Fortran-written binary stream file into a numpy array.

    Arguments
    ----------
    filename : str
        Path to the .bin file.
    timeIndex: int
        Stage-specific time step at which data is to be loaded
    shape : tuple of ints
        Shape of the array as (nx, ny) or (nx, ny, nz).
    dtype : numpy dtype, optional
        Data type. Default is float64 (Fortran double precision).

    Returns
    -------
    numpy array with given shape, Fortran memory order.
    """
    if len(shape) == 2:
        arr = np.fromfile(filename, dtype=dtype) \
                .reshape(shape, order='F')
    elif len(shape) == 3:
        arr = np.fromfile(filename, dtype=dtype) \
            .reshape(shape, order='F')
        arr = arr[:,:,timeIndex]
    return arr
    
    
    
def loadSnapshot(scaleIndex, timeIndex, baseNames, shapes):
    """
    Load the four arrays corresponding to one simulation time index.

    Arguments
    ----------
    time_index : int
        Snapshot number (0, 1, 2, ...).
    base_names : list of str
        List of 4 base filenames, e.g.
        ["firstFilename", "secondFilename", "thirdFilename", "fourthFilename"]
    shape : list of tuples of ints
        List elements are arrays shape (nx, ny) or (nx, ny, nz).

    Returns
    -------
    A tuple (A, B, C, D) of numpy arrays.
    """
    arrays = []
    for base, shape in zip(baseNames, shapes):
        fname = f"{base}{scaleIndex}.bin"
        arr = loadBinArray(fname, timeIndex, shape)
        arrays.append(arr)

    return tuple(arrays)

def loadParameters(filename):
    """
    Load parameters (nx, ny, nt, offset) from a .txt file where the values are
    stored in a column (one value per line).

    Parameters
    ----------
    filename : str
        Path to the .txt file.

    Returns
    -------
    tuple
        A tuple (nx, ny, nt, offset) containing the parameter values.
    """
    with open(filename, 'r') as f:
        params = [float(line.strip()) for line in f if line.strip()]
        
        if len(params) != 4:
            raise ValueError(f"Expected 4 parameters in file, but found {len(params)}.")
        return tuple(params)
