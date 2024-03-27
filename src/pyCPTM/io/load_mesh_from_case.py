from pathlib import Path
import re
import numpy as np

def load_mesh_from_case(pathToCase=[], maxFaceNumber=10**6):
    if isinstance(pathToCase,str):
        # String to constant/polyMesh folder
        data_dir = Path(pathToCase + r"\constant\polyMesh")
    else:
        # has to be path to polymesh then
        data_dir = pathToCase

    startFacePattern = re.compile(r"\bstartFace\b")
    curlyBracketPattern = re.compile(r"\{")

    # Only extract data from within brackets:
    # bracketDecider == 0 outside of brackets
    # bracketDecider == 1 inside of brackets
    bracketDecider = 0

    # owner = np.array([],dtype=np.dtype('uintc'))
    # neighbour = np.array([],dtype=np.dtype('uintc'))

    cellArray = np.array(np.arange(maxFaceNumber), dtype=np.dtype("uintc"))

    lines = []
    boundaryNames = []
    startFace = []

    for file in data_dir.iterdir():
        if file.name == "boundary":
            with open(file) as f:
                for line in f:
                    if line == "(\n" or line == ")\n":
                        bracketDecider += 1 * (line == "(\n") - 1 * (line == ")\n")
                    lines.append(line)
                    if re.search(curlyBracketPattern, line) and bracketDecider == 1:
                        boundaryNames.append(re.split("\n| ", lines[-2])[-2])
                    if re.search(startFacePattern, line):
                        startFace.append(int(re.split(";| ", line)[-2]))
                        # print(f'Boundary:{boundaryNames[-1]}, startFace: {startFace[-1]}')

    # Iterator i
    i = 0

    # Export Owner and Neighbour
    for file in data_dir.iterdir():
        if file.name == "owner" or file.name == "neighbour":
            lineCounter = 0
            with open(file) as f:
                for line in f:
                    # if lineCounter % 100 == 0:
                    #     print(lineCounter)
                    lineCounter += 1
                    if line == "(\n" or line == ")\n":
                        bracketDecider += 1 * (line == "(\n") - 1 * (line == ")\n")
                    elif bracketDecider == 1:
                        cellArray[i] = np.uintc(line)
                        # if file.name == "owner":
                        # owner = np.append(owner, np.uintc(line))
                        # if file.name == "neighbour":
                        # neighbour = np.append(neighbour, np.uintc(line))
                        i += 1
            if file.name == "owner":
                cellArray = cellArray[:i]
                i = 0
                owner = cellArray
                cellArray = np.array(np.arange(maxFaceNumber), dtype=np.dtype("uintc"))
            elif file.name == "neighbour":
                cellArray = cellArray[:i]
                i = 0
                neighbour = cellArray
                cellArray = np.array(np.arange(maxFaceNumber), dtype=np.dtype("uintc"))

    return owner, neighbour, (boundaryNames, startFace)
