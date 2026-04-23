"""
simulation_data.py
------------------
Data structure for managing Fortran-generated binary arrays across many loop
iterations.  Arrays are lazy-loaded (read on first access, then freed) unless
the user explicitly pins certain indices into a persistent in-memory cache.

Expected file naming convention on disk:
    moment_1.bin, momentRate_1.bin, magnitude_1.bin
    moment_2.bin, momentRate_2.bin, magnitude_2.bin
    ...
    moment_N.bin, momentRate_N.bin, magnitude_N.bin
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator
from plotting.plotProfiles import plotProfiles

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# One step (all three arrays for a single loop index)
# ---------------------------------------------------------------------------

@dataclass
class SimStep:
    """Holds the three arrays for one loop iteration."""
    index: int
    moment: np.ndarray
    momentRate: np.ndarray
    magnitude: np.ndarray

    # ------------------------------------------------------------------
    # Convenience analysis helpers on a single step
    # ------------------------------------------------------------------
    def max_values(self) -> dict[str, float]:
        return {"moment": float(self.moment.max()),
                "momentRate": float(self.momentRate.max()),
                "magnitude": float(self.magnitude.max())}

    def min_values(self) -> dict[str, float]:
        return {"moment": float(self.moment.min()),
                "momentRate": float(self.momentRate.min()),
                "magnitude": float(self.magnitude.min())}

    def summary(self) -> str:
        lines = [f"Step {self.index}:"]
        for name in ("moment", "momentRate", "magnitude"):
            arr = getattr(self, name)
            lines.append(
                f"  {name}: shape={arr.shape}  min={arr.min():.4g}"
                f"  max={arr.max():.4g}  mean={arr.mean():.4g}"
            )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main data-management class
# ---------------------------------------------------------------------------

class SimulationData:
    """
    Manages lazy-loading of binary arrays produced by a Fortran code.

    Parameters
    ----------
    data_dir : str | Path
        Directory that contains all the .bin files.
    n_steps : int
        Total number of loop iterations (files are indexed 1 … n_steps).
    dtype : np.dtype
        NumPy dtype matching the Fortran REAL kind (default: float64).
    shape : tuple[int, ...] | None
        Array shape expected for each variable.  Pass None to keep a 1-D
        flat array (you can always reshape later).
    var_names : tuple[str, str, str]
        The three variable prefixes used in file names (default: moment/B/C).
    """

    def __init__(
        self,
        data_dir: str | Path,
        n_steps: int,
        dtype: np.dtype = np.float64,
        shape: tuple[int, ...] | None = None,
        var_names: tuple[str, str, str] = ("moment", "momentRate", "magnitude"),
    ):
        self.data_dir = Path(data_dir)
        self.n_steps = n_steps
        self.dtype = dtype
        self.shape = shape
        self.var_names = var_names

        # Indices that are kept permanently in memory.
        self._pinned: dict[int, SimStep] = {}

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _bin_path(self, var: str, index: int) -> Path:
        return self.data_dir / f"{var}_{index}.bin"

    def _load_array(self, var: str, index: int) -> np.ndarray:
        path = self._bin_path(var, index)
        arr = np.fromfile(path, dtype=self.dtype)
        if self.shape is not None:
            arr = arr.reshape(self.shape)
        return arr

    def _load_step(self, index: int) -> SimStep:
        arrays = {v: self._load_array(v, index) for v in self.var_names}
        return SimStep(index=index, **arrays)

    # ------------------------------------------------------------------
    # Pinning API
    # ------------------------------------------------------------------

    def pin(self, indices: list[int]) -> None:
        """
        Load and permanently cache the steps at the given indices.
        Useful for steps you want to inspect / plot repeatedly without
        re-reading from disk.
        """
        for i in indices:
            if i not in self._pinned:
                self._pinned[i] = self._load_step(i)
        print(f"Pinned indices: {sorted(self._pinned)}")

    def unpin(self, indices: list[int] | None = None) -> None:
        """
        Release pinned steps from memory.
        Pass None (default) to unpin everything.
        """
        if indices is None:
            self._pinned.clear()
            print("All indices unpinned.")
        else:
            for i in indices:
                self._pinned.pop(i, None)
            print(f"Unpinned {indices}. Remaining: {sorted(self._pinned)}")

    def pinned_indices(self) -> list[int]:
        return sorted(self._pinned)

    # ------------------------------------------------------------------
    # Access API
    # ------------------------------------------------------------------

    def __getitem__(self, index: int) -> SimStep:
        """
        Return the SimStep for `index`.
        Uses the pinned cache if available; otherwise loads from disk
        and immediately discards after return (not cached).
        """
        if index in self._pinned:
            return self._pinned[index]
        return self._load_step(index)

    def iter_steps(
        self,
        start: int = 1,
        stop: int | None = None,
        step: int = 1,
    ) -> Iterator[SimStep]:
        """
        Iterate over steps [start, stop) with a given stride.
        Pinned steps are returned from cache; others are loaded and freed.

        Parameters
        ----------
        start : int  First index (default 1).
        stop  : int  Exclusive upper bound (default n_steps + 1).
        step  : int  Stride (default 1).
        """
        if stop is None:
            stop = self.n_steps + 1
        for i in range(start, stop, step):
            yield self[i]

    # ------------------------------------------------------------------
    # Global analysis  (stream through disk, never hold all data at once)
    # ------------------------------------------------------------------

    def global_max(self, var: str) -> float:
        """Maximum value of `var` across ALL steps."""
        return max(
            float(self._load_array(var, i).max())
            for i in range(1, self.n_steps + 1)
        )

    def global_min(self, var: str) -> float:
        """Minimum value of `var` across ALL steps."""
        return min(
            float(self._load_array(var, i).min())
            for i in range(1, self.n_steps + 1)
        )

    def stats_over_steps(
        self,
        var: str,
        start: int = 1,
        stop: int | None = None,
    ) -> dict[str, np.ndarray]:
        """
        Compute per-step max, min, and mean of `var` for steps [start, stop).
        Returns a dict with keys 'indices', 'max', 'min', 'mean'.
        """
        if stop is None:
            stop = self.n_steps + 1
        indices, maxs, mins, means = [], [], [], []
        for i in range(start, stop):
            arr = self._load_array(var, i)
            indices.append(i)
            maxs.append(arr.max())
            mins.append(arr.min())
            means.append(arr.mean())
        return {
            "indices": np.array(indices),
            "max":     np.array(maxs),
            "min":     np.array(mins),
            "mean":    np.array(means),
        }

    # ------------------------------------------------------------------
    # Dunder niceties
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return self.n_steps

    def __repr__(self) -> str:
        return (
            f"SimulationData(n_steps={self.n_steps}, "
            f"dtype={self.dtype}, shape={self.shape}, "
            f"pinned={self.pinned_indices()})"
        )
    
# ------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------
    
def plotHistogram(data, bins, figLabels = None, figTitle = None, figSubTitle = None, textCount = None):
    fig = plt.figure(figsize=(5, 8))
    ax = fig.add_subplot(1,1,1)
    ax.hist(data, bins)
    
    fig.supxlabel(figLabels[0])
    fig.supylabel(figLabels[1])
    
    if figTitle != 'None':
        fig.suptitle(figTitle, fontsize=16)
        plt.title(figSubTitle, fontsize = 12)
        
    #Add text
    if textCount != 'None':
        ax.text(
            0.5, 0.8,
            f"Nr. of interrupted large events: {textCount}",
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top')
        
    plt.grid()
    



# ===========================================================================
# Example usage (run this file directly to see it in action with dummy data)
# ===========================================================================

def _make_dummy_files(data_dir: Path, n_steps: int, shape: tuple) -> None:
    """Create small random .bin files for demonstration purposes."""
    data_dir.mkdir(exist_ok=True)
    rng = np.random.default_rng(0)
    for i in range(1, n_steps + 1):
        for var in ("moment", "momentRate", "magnitude"):
            arr = rng.random(shape).astype(np.float64)
            arr.tofile(data_dir / f"{var}_{i}.bin")


if __name__ == "__main__":
    #import tempfile
    data_path = "/home/viktor/Dokumente/Doktor/ENS_BRGM/Code/data/asperity_statistics/6_4_1000_rn"
    

    N_STEPS = 116
    SHAPE   = (501)   # adjust to match your actual array dimensions


    data_dir = Path(data_path)

    #print("=== Creating dummy .bin files ===")
    #_make_dummy_files(data_dir, N_STEPS, SHAPE)

    # -----------------------------------------------------------
    # 1. Initialise
    # -----------------------------------------------------------
    sim = SimulationData(data_dir, n_steps=N_STEPS, shape=SHAPE)
    print(sim)

    # -----------------------------------------------------------
    # 2. Pin a few steps you want to keep in memory
    # -----------------------------------------------------------
    sim.pin([116])

    # -----------------------------------------------------------
    # 3. Quick access to a pinned step (no disk read)
    # -----------------------------------------------------------
    step116 = sim[116]
    print("\n--- Single step access ---")
    print(step116.summary())
    print("Max values:", step116.max_values())

    
    # -----------------------------------------------------------
    # 4. My own stuff
    # -----------------------------------------------------------
    eventMagnitudes = []
    unfinishedCount = 0
    for step in sim.iter_steps():
        step_max = step.magnitude.max()
        eventMagnitudes.append(step_max)
        
        if step.magnitude[-1] != 0:
            unfinishedCount = unfinishedCount + 1
        
    histoBins = np.linspace(0, 3.5, 10)
    histoLabels = ["magnitude", "N"]
    histoTitle = 'Magnitude Distribution'
    plotHistogram(eventMagnitudes, histoBins, figLabels = histoLabels, figTitle = histoTitle, figSubTitle = "Timesteps: " + str(SHAPE), textCount = unfinishedCount)
    
    plotProfiles([eventMagnitudes], ["ihypo", "Magnitudes"], "mag")
    
    
    
# -------------------------------------------------------------
# Some cool things
# -------------------------------------------------------------

# # -----------------------------------------------------------
# # 1. Iterate over ALL steps, computing a running maximum.
# #    Only pinned steps stay in memory; the rest are freed
# #    automatically after each loop body.
# # -----------------------------------------------------------
# print("\n--- Streaming loop over all steps ---")
# running_max = -np.inf
# for step in sim.iter_steps():
#     step_max = step.moment.max()
#     if step_max > running_max:
#         running_max = step_max
# print(f"Global max of moment (streaming): {running_max:.6f}")

# # -----------------------------------------------------------
# # 2. Convenience global methods (no data stays in RAM)
# # -----------------------------------------------------------
# print("\n--- Global statistics ---")
# print(f"global_max('moment') = {sim.global_max('moment'):.6f}")
# print(f"global_min('momentRate') = {sim.global_min('momentRate'):.6f}")

# # -----------------------------------------------------------
# # 3. Per-step statistics for a range (good for plotting)
# # -----------------------------------------------------------
# print("\n--- Per-step stats for magnitude, steps 1-10 ---")
# stats = sim.stats_over_steps("magnitude", start=1, stop=11)
# for i, mx, mn, mu in zip(
#     stats["indices"], stats["max"], stats["min"], stats["mean"]
# ):
#     print(f"  step {i:3d}:  max={mx:.4f}  min={mn:.4f}  mean={mu:.4f}")

# # -----------------------------------------------------------
# # 4. Access a pinned step for in-depth work (e.g. plotting)
# # -----------------------------------------------------------
# print("\n--- Accessing pinned step 20 for in-depth analysis ---")
# s20 = sim[20]
# # Example: histogram of moment values
# counts, edges = np.histogram(s20.moment, bins=8)
# print("moment histogram (8 bins):")
# for c, lo, hi in zip(counts, edges[:-1], edges[1:]):
#     bar = "#" * (c // 2)
#     print(f"  [{lo:.3f}, {hi:.3f})  {bar} ({c})")

# # -----------------------------------------------------------
# # 5. Unpin to free memory when no longer needed
# # -----------------------------------------------------------
# sim.unpin([1, 5])
# print(f"\nAfter unpinning 1 and 5: {sim}")

        
    
