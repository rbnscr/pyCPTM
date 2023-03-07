# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 11:28:28 2022

@author: schro22
"""
# https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/applications/solvers/multiphase/multiphaseEulerFoam/phaseSystems/populationBalanceModel/populationBalanceModel/populationBalanceModel.C

import numpy as np
from .breakupModel import BreakupModel
from .daughterSizeDistribution import DaughterSizeDistribution

class PopulationBalanceModel:

    def __init__(self, sizeGroups, phaseSystemDict):
        self.sizeGroups = sizeGroups # diameter [m]

    def birth_by_coalescence(self,i,j):
        pass

    def death_by_coalescence(self,i,j):
        pass

    def birth_by_breakup(self,i,j):
        pass

    def death_by_breakup(self,i,breakupDict):
        pass

    def calc_delta(self):
        pass

    def eta(self,i,v):
        x0 = self.sizeGroups[0]
        xi = self.sizeGroups[i]
        xm = self.sizeGroups[-1]

        lowerBoundary = x0
        upperBoundary = xm

        if i != 0:
            lowerBoundary = self.sizeGroups[i-1]
        if i != len(self.sizeGroups)-1:
            upperBoundary = self.sizeGroups[i+1]
        if np.logical_or(np.logical_and(i == 0, v < x0), np.logical_and(i == len(self.sizeGroups)-1, v > xm)):
            return v/xi
        elif np.logical_or(v < lowerBoundary, v > upperBoundary):
            return 0
        elif v == xi:
            return 1
        elif v > xi:
            return (upperBoundary-v)/(upperBoundary-xi)
        else:
            return (v-lowerBoundary)/(xi-lowerBoundary)

    def solve(self):
        pass

