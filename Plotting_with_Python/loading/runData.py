"""
Class to generalize the loading of model runs.
Contains run name, run data directory, and the resolution in time and space used 
for a specific run.


Created on Mon Dec 15 10:55:11 2025

@author: viktor
"""
import numpy as np
from pathlib import Path
from .loaders import loadBinArray

class runData:
    def __init__(self, name, runDir, variables, nx, ny, nt, offset, runTime=None, kernelTime=None):
        """
        Parameters
        ----------
        name : str
            Identifier for the run (e.g. 'run_01')
        runDir : Path or str
            Directory where this run's data is stored
        shape : tuple
            Array shape (nx, ny) or (nx, ny, nz)
        """
        self.name = name
        self.runDir = Path(runDir)
        self.variables = variables
        self.sizes = {
            "nx": nx,
            "ny": ny,
            "nt": nt}
        self.offset = offset
        if runTime is not None:
            self.runTime = runTime
        if kernelTime is not None:
            self.kernelTime = kernelTime
        
    def filePath(self, var, scale):
        return self.runDir / f"{var}{scale}.bin"
    
    def shapeFromSpec(self, spec):
        return tuple(self.sizes[d] for d in spec.dims)
        
    def load(self, var, scale, timeInd = None, dtype=np.float64):
        """
        Load a variable for a given scale and time if applicable

        Parameters
        ----------
        var : str
            Variable name
        scale : int
            Scale index
        time_index : int, optional
            If provided and variable is time-dependent, returns a 2D slice
        dtype : np.dtype
            Data type for numpy

        Returns
        -------
        np.ndarray
            2D slice if time_index is given, else full array
        """
        spec = self.variables[var]
        path = self.filePath(var, scale)
        shape = self.shapeFromSpec(spec)
        

        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
            
        arr = loadBinArray(path, shape, dtype)

        # for 3D array: return 2D slice if requested and variable has time dimension
        if len(shape) == 3: 
            if spec.timeDependent and timeInd is not None:
                return arr[:, :, timeInd]
            elif spec.timeDependent and timeInd is None:
                return arr[:, :, :]
    
        return arr