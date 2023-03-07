# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 12:11:40 2022

@author: schro22
"""
import numpy as np
from scipy.special import erfc
import math

# openfoam implementation https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/applications/solvers/multiphase/multiphaseEulerFoam/phaseSystems/populationBalanceModel/breakupModels/Laakkonen/Laakkonen.C
# Laakkonen et al. (2007) https://doi.org/10.1016/j.ces.2006.10.006

class BreakupModel:
    breakupModels = {
        "laakkonen": "laakkonenBreakupRate",
        "test": "testBreakupRate",
        "keyToCall": "methodToCall",
    }

    def __init__(self, v, modelDict):
        self.v = v
        self.modelDict = modelDict
        modelName = modelDict["model"]
        self.call_model(modelName)

    def call_model(self, arg):
        return getattr(self, self.breakupModels[arg])(self.v, self.modelDict)

    def laakkonenBreakupRate(self, v, modelDict):
        # Read keys from modelDict
        sigma = modelDict["sigma"]
        epsilonCont = modelDict["epsilonCont"]
        muCont = modelDict["muCont"]
        rhoCont = modelDict["rhoCont"]
        rhoDisp = modelDict["rhoDisp"]

        # Set defaults
        c1default = 2.25
        c2default = 0.04
        c3default = 0.01
        c1 = modelDict.get("c1")
        c2 = modelDict.get("c2")
        c3 = modelDict.get("c3")
        if c1 == None:
            c1 = c1default
        if c2 == None:
            c2 = c2default
        if c3 == None:
            c3 = c3default

        # volume based

        # breakupRate g(v) | breakage frequency [1/s]
        # \sigma      |  Surface tension [N/m]
        # v_i         |  Volume of mother bubble i [m3]
        # \epsilon_c  |  Turbulent dissipation rate of continuous phase [m^2/s^3]
        # \mu_c       |  Molecular dynamic viscosity of liquid phase [Pa s]
        # \rho_c      |  Density of continuous phase [kg/m^3]
        # \rho_d      |  Density of disperse phase [kg/m^3]

        # Property     | Description             | Required      | Default value
        # C1           | coefficient C1          | no            | 2.25
        # C2           | coefficient C2          | no            | 0.04
        # C3           | coefficient C3          | no            | 0.01
        # daughterSizeDistributionModel | inh. from breakupModel | inherited

        breakupRate = (
            c1
            * np.cbrt(epsilonCont)
            * erfc(
                np.sqrt(
                    (c2 * sigma)
                    / (rhoCont * epsilonCont ** (2 / 3) * (6 * v / math.pi) ** (5 / 9))
                    + (c3 * muCont)
                    / (
                        np.sqrt(rhoCont * rhoDisp)
                        * np.cbrt(epsilonCont)
                        * (6 * v / math.pi) ** (4 / 9)
                    )
                )
            )
        )
        self.breakupRate = breakupRate

    def testBreakupRate(self, v, modelDict):
        self.breakupRate = 1


# dictToClass = {"model":"laakkonen",
#                 "sigma":0.01,
#                 "epsilonCont":0.01,
#                 "muCont":0.01,
#                 "rhoCont":0.01,
#                 "rhoDisp":0.01
#                 }
#
# print(BreakupModel(0.01, dictToClass).breakupRate)


# startSize = 0.0005
# endSize = 0.002
# steps = 20

# v = np.arange(startSize, endSize, (endSize - startSize) / steps)

# sigma = 0.072
# epsilonC = 10
# muC = 1e-3
# rhoC = 1000
# rhoD = 1.225

# test = laakkonenBreakupRate(v, sigma, epsilonC, muC, rhoC, rhoD,c1 = 1)
# print(test)
