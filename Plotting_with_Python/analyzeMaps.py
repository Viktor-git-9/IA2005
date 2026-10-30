#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 12:43:07 2026
Script to read Rinis from infofiles and reconstruct the asperity sizes.

@author: viktor
"""

import numpy as np
import matplotlib.pyplot as plt
from plotting.plotContour import plotContours # import plotting functions
from pathlib import Path

def loadBinArray(filename, shape, dtype=np.float64):
    """
    Load a Fortran-written binary stream file into a numpy array.

    Arguments
    ----------
    filename : str
        Path to the .bin file.
    timeIndex: int
        Stage-specific time step at which data is to be loaded
    shape : tuple of ints
        Desired shape of the array
    dtype : numpy dtype, optional
        Data type. Default is float64 (Fortran double precision).

    Returns
    -------
    numpy array with given shape, Fortran memory order.
    """
    arr = np.fromfile(filename, dtype=dtype) \
                .reshape(shape, order='F')
    arr = np.array(arr)
    return arr

root = Path("/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/smooth_distributions/1006/100asps_3MPa/infofiles/")

shape = 100 # asperity count of maps
correctionFactor = 3/2 # Get from Rinis to asperity radii. Pay mind to this!
asperitySizes = {}

for subfolder in sorted(root.iterdir()):
    if not subfolder.is_dir():
        continue

    bin_path = subfolder / "Rinis.bin"
    if not bin_path.exists():
        print(f"Skipping {subfolder.name}: no Rinis.bin found")
        continue

    asperitySizes[subfolder.name] = loadBinArray(str(bin_path), shape) * 3/2
    
#%% Plotting

names = list(asperitySizes.keys())
vectors = [np.asarray(asperitySizes[n]).ravel() for n in names]

# Shared bin edges spanning all data
all_values = np.concatenate(vectors)
n_bins = 12
edges = np.linspace(all_values.min(), all_values.max(), n_bins + 1)
centers = 0.5 * (edges[:-1] + edges[1:])

# Histogram each vector with the same edges
counts = np.vstack([np.histogram(v, bins=edges)[0] for v in vectors])  # (n_realizations, n_bins)

# One distinct color per vector (tab10 has exactly 10 colors)
colors = plt.cm.tab10(np.arange(len(vectors)))

fig, ax = plt.subplots(figsize=(9, 5))

# Every realization, faint
for v, name, color in zip(vectors, names, colors):
    counts, _ = np.histogram(v, bins=edges)
    ax.plot(centers, counts, color=color, lw=1.8, marker="o", ms=4, label=name)

# Median and spread across realizations
# median = np.median(counts, axis=0)
# lo, hi = np.percentile(counts, [5, 95], axis=0)
# ax.fill_between(centers, lo, hi, color="steelblue", alpha=0.25, label="5-95th percentile")
# ax.plot(centers, median, color="navy", lw=2, label="Median")

ax.set_title("Asperity distributions: 10 Maps with 10 asperities")
ax.set_xlabel("Asperity Radius [els.]")
ax.set_ylabel("Count")
ax.legend(title="Subfolder", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.grid()
plt.tight_layout()
plt.show()