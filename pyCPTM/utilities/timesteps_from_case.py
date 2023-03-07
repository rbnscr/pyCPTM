import numpy as np
from pathlib import Path
import os

def timesteps_from_case(pathToCase=[], numeric=True):
    pathToCase = Path(pathToCase)
    if numeric:
        timesteps = np.array([])
    elif not numeric:
        timesteps = np.array([], dtype=np.dtype("U"))

    for subdata_dir in pathToCase.iterdir():
        subdata_dir = str(subdata_dir)
        subdata_dir = subdata_dir.split(os.path.sep)[-1]
        try:
            float(subdata_dir)
            if numeric:
                timestep = float(subdata_dir)
            elif not numeric:
                timestep = subdata_dir
            timesteps = np.append(timesteps, timestep)
        except ValueError:
            # print(f"{subdata_dir} is not a timestep.")
            pass
    return timesteps
