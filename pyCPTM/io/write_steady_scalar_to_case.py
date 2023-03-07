# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 11:10:43 2022

@author: schro22
"""
import numpy as np
from pyCPTM.utilities import timesteps_from_case_as_string
import os

def write_steady_scalar_to_case(pathToCase=[], scalarName=[], scalar=np.array([])):
    # Number of entries, that go into the file
    numberOfEntries = len(scalar)
    # Get timesteps from case
    timestepsFromCase = timesteps_from_case_as_string(pathToCase)
    # Create and write to new text file
    for timestep in timestepsFromCase:
        with open(pathToCase + f"{os.path.sep}{timestep}{os.path.sep}{scalarName}", "w") as f:
            f.write(_write_header(scalarName, timestep, numberOfEntries))
            for row in scalar:
                f.write(f"{row}\n")
            f.write(_write_foot())
            f.close()


def _write_header(objectName, timestep, numberOfEntries):
    header_string = (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "  =========                 |\n"
        "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n"
        "   \\\\    /   O peration     | Website:  https://openfoam.org \n"
        "    \\\\  /    A nd           | Version:  dev\n"
        "     \\\\/     M anipulation  |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n"
        "{\n"
        "    format      ascii;\n"
        "    class       volScalarField::Internal;\n"
        f'    location    "{timestep}";\n'
        f"    object      {objectName};\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
        "\n"
        "dimensions      [0 0 0 0 0 0 0];\n"
        "\n"
        "value   nonuniform List<scalar>\n"
        f"{numberOfEntries}\n"
        "(\n"
    )

    return header_string


def _write_foot():
    foot_string = (
        ")\n"
        ";\n"
        "\n"
        "\n"
        "// ************************************************************************* //\n"
    )

    return foot_string


# write_scalar_to_case(
#     pathToCase="M:/Schroeder/Code/pyCPTM/data/stirredTankDynamicMesh",
#     scalarName="test",
#     scalar=np.array([1, 2, 3, 4, 5]),
# )
