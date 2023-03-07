#%% Load modules

from pyCPTM.cptm import CPTMSimple

import pandas as pd

import numpy as np

import pytest

testPath = "test/testCase"

owner = pd.read_csv(
    testPath + "/testMesh.csv", usecols=["owner"], sep=";"
).to_numpy()

neighbour = pd.read_csv(
    testPath + "/testMesh.csv", usecols=["neighbour"], sep=";"
).to_numpy()

internalMesh = np.array([owner[: len(neighbour)], neighbour]).T[0, :, :]

phi = pd.read_csv(
    testPath + "/testMesh.csv", usecols=["phi"], sep=";", decimal=","
).to_numpy()[:, 0]

coordinates = pd.read_csv(
    testPath + "/testData.csv", usecols=["cx", "cy", "cz"], sep=";"
).to_numpy()

coordinates = coordinates[:16, :]

cell_volume = pd.read_csv(
    testPath + "/testData.csv", usecols=["V"], sep=";"
).to_numpy()

cell_volume = cell_volume[:16, :]

cellIndex = pd.read_csv(
    testPath + "/testData.csv", usecols=["idx"], sep=";"
).to_numpy()

cellIndex = cellIndex[:16, :]

cellData = pd.read_csv(
    testPath + "/testData.csv", usecols=["U_x", "U_y", "U_z"], sep=";", decimal=","
)



cm = CPTMSimple(internalMesh, phi, cell_volume, coordinates, cellData=cellData, cptIdx=cellIndex[:, 0])
cm.projection()
testIdx = [0,3,1,1,3,3,1,1,2,5,6,6,7,7,4,6]
# testIdx might change depending on renumbering after projection or not


checkerboardIdx = np.array([0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0])
checkerboardModel = CPTMSimple(internalMesh, phi, cell_volume, coordinates, cellData=cellData, cptIdx=checkerboardIdx)
checkerboardModel.projection()

checkerboardTestIdx = [0,8,1,9,10,2,11,3,4,12,5,13,14,6,15,7]



def test_answer():
    assert cm.nCpt == 8
    assert cm.nCell == 16
    assert (cm.cptIdx == testIdx).all()
    
    assert checkerboardModel.nCpt == 16
    assert (checkerboardModel.cptIdx == checkerboardTestIdx).all()