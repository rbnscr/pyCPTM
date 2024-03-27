from pathlib import Path
import numpy as np

from pyCPTM.utilities import find_matches
from pyCPTM.utilities import timesteps_from_case_as_string


def load_phi(
    pathToCase=[],
    chosenPhase=[],
    chosenTimesteps=[],
    maxFaceNumber=2.5 * 10**6,
    precision="single",
):
    # chosenTimesteps as string

    if not isinstance(maxFaceNumber, int):
        maxFaceNumber = int(maxFaceNumber)

    if isinstance(pathToCase,str):
        # String to constant/polyMesh folder
        data_dir = Path(pathToCase)
    else:
        # Must be a Path instance directing to the case folder
        data_dir = pathToCase

    # Gather timesteps as array from a case directory as array of string
    caseTimestepsString = timesteps_from_case_as_string(data_dir)

    if chosenTimesteps:
        usedTimesteps = find_matches(chosenTimesteps, caseTimestepsString)
    elif not chosenTimesteps:
        usedTimesteps = caseTimestepsString
    # Remove timestep 0
    if any(usedTimesteps == "0"):
        usedTimesteps = np.delete(usedTimesteps, usedTimesteps == "0")

    if chosenPhase:
        phases = [f"phi.{phase}" for phase in chosenPhase]
    elif not chosenPhase:
        phases = ["phi"]

    # Only extract data from within brackets:
    # bracketDecider == 0 outside of brackets
    # bracketDecider == 1 inside of brackets
    bracketDecider = 0

    # init phi array
    phi = np.zeros(
        shape=(maxFaceNumber, len(phases), len(usedTimesteps)),
        dtype=np.dtype(precision),
    )
    timestepPhase = []  # for logging

    # Determine phi* pathes/files to look at
    for timestep in usedTimesteps:
        for phase in phases:
            phiPath = data_dir / f"{timestep}" / f"{phase}"
            timestepPhase.append(f"{timestep}_{phase}")
            i = 0
            if phiPath.exists() and timestep != "0":
                with open(phiPath) as f:
                    for line in f:
                        if line == "(\n" or line == ")\n":
                            bracketDecider += 1 * (line == "(\n") - 1 * (line == ")\n")
                        elif bracketDecider == 1:
                            line = np.array([line], dtype=np.dtype(precision))
                            phi[
                                i,
                                phases.index(phase),
                                np.where(usedTimesteps == timestep)[0][0],
                            ] = line
                            i += 1
    try:
        phi = phi[:i, :, :]
    except UnboundLocalError:
        print(
            "Something went wrong. Make sure, that there are timesteps in the case folder."
        )
        return
    return phi, timestepPhase
