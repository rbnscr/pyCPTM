# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 11:10:43 2022

@author: schro22
"""
import numpy as np
from pyCPTM.utilities import timesteps_from_case_as_string
from pyCPTM.io.write_steady_scalar_to_case import _write_header, _write_foot
import os

def write_transient_scalar_to_case(pathToCase=[], scalarName=[], scalar=np.array([])):
    # Number of entries, that go into the file
    numberOfEntries = len(scalar[:,])
    # Get timesteps from case
    timestepsFromCase = timesteps_from_case_as_string(pathToCase)
    # Check, if there are as many timesteps in the case as arrays to write to file
    if numberOfEntries != len(timestepsFromCase):
        return print("Number of timesteps does not match size of array.")
    # Create and write to new text file
    timestep_iter = 0
    for timestep in timestepsFromCase:
        with open(pathToCase + f"{os.path.sep}{timestep}{os.path.sep}{scalarName}", "w") as f:
            f.write(_write_header(scalarName, timestep, len(scalar[timestep_iter, :])))
            for row in scalar[timestep_iter, :]:
                # for i in range(0, len(scalar[timestep_iter, :])):
                f.write(f"{row}\n")
                # f.write(f"{scalar[timestep_iter, i]}\n")
            f.write(_write_foot())
            f.close()
        print(
            f"Timestep {timestep_iter+1} / {numberOfEntries} done."
        )  # Might be deleted at a later stage, or toggled with option.
        timestep_iter += 1


# test_array = np.array(
#     [
#         [1, 3, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#         [1, 2, 3, 4],
#     ]
# )

# write_transient_scalar_to_case(
#     pathToCase="M:/Schroeder/Code/pyCPTM/data/passiveScalarPitzDaily",
#     scalarName="Test",
#     scalar=test_array,
# )

# write_transient_scalar_to_case(
#     pathToCase=absolutePath, scalarName="cpt_tracer", scalar=cx
# )
