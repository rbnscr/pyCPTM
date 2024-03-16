# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:23:32 2022

@author: schro22
"""

import math
import numpy as np

class CoalescenceModel:
    coalescenceModels = {
        "luo": "luoCoalescenceModel",
        "keyToCall": "methodToCall",
    }
    
    def __init__(self, vi, vj, modelDict):
        self.vi = vi
        self.vj = vj
        self.modelDict = modelDict
        modelName = modelDict["model"]
        self.call_model(modelName)
    
    def call_model(self, arg):
        return getattr(self, self.coalescenceModels[arg])(self.vi, self.vj, self.modelDict)

    def luoCoalescenceModel(self,
        vi, vj, modelDict
    ):
        # Read keys from modelDict
        sigma = modelDict["sigma"]
        epsilonCont = modelDict["epsilonCont"]
        rhoCont = modelDict["rhoCont"]
        rhoDisp = modelDict["rhoDisp"]
        Cvm = modelDict["Cvm"]
        
        # Set defaults
        C1default = 1
        betadefault = 2
        C1 = modelDict.get("C1")
        beta = modelDict.get("beta")
        if C1 == None:
            C1 = C1default
        if beta == None:
            beta = betadefault
        # https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/applications/solvers/multiphase/multiphaseEulerFoam/phaseSystems/populationBalanceModel/coalescenceModels/Luo/Luo.H
        
        # d_i         |  Diameter of bubble i [m]
        # d_j         |  Diameter of bubble j [m]
        # u_{ij}      |  Mean approach velocity [m/s]
        # \xi_{ij}    |  Bubble size ratio [-]
        # \rho_d      |  Density of dispersed phase [kg/m^3]
        # \rho_c      |  Density of continuous phase [kg/m^3]
        # \sigma      |  Surface tension [N/m]
        # C_{vm}      |  Virtual mass coefficient [-]
        # C_1         |  Coefficient [-]
        # \beta       |  Coefficient [-]
        # \epsilon_c  |  Continuous phase turbulent dissipation rate [m^2/s^3]
    
        # beta         | Coefficient beta        | no          | 2.0
        # C1           | Coefficient C1          | no          | 1.0
    
        di = volumeToDiameter(vi)
        dj = volumeToDiameter(vj)
        xi = di / dj
        uij = meanApproachVelocity(di, epsilonCont, xi, beta)
    
        # SPECIFIC coalesence Rate [m³/s] ? Oder ist es [1/(m³*s)]?
    
        coalescenceRate = (
            math.pi
            / 4
            * (di + dj) ** 2
            * uij
            * np.exp(
                -C1
                * np.sqrt(0.75 * (1 + xi ** 2) * (1 + xi ** 3))
                / (np.sqrt(rhoDisp / rhoCont + Cvm) * (1 + xi) ** 3)
                * np.sqrt(rhoCont * di * uij ** 2 / sigma)
            )
        )
    
        self.coalescenceRate = coalescenceRate
    

def volumeToDiameter(sphereVolume):
    return np.cbrt(6 * sphereVolume / math.pi)


def meanApproachVelocity(d, epsilonCont, xi, beta):
    return np.sqrt(beta) * np.cbrt(epsilonCont * d) * np.sqrt(1 + xi ** (-2 / 3))


# def luoCoalescenceModelMatrix(
#     bubbleSizes, sigma, epsilonCont, rhoCont, rhoDisp, Cvm, C1=1, beta=2
# ):
#     outputMatrix = np.zeros((len(bubbleSizes), len(bubbleSizes)))
#     for i in range(len(bubbleSizes)):
#         for j in range(len(bubbleSizes)):
#             outputMatrix[i, j] = luoCoalescenceModel(
#                 bubbleSizes[i],
#                 bubbleSizes[j],
#                 sigma,
#                 epsilonCont,
#                 rhoCont,
#                 rhoDisp,
#                 Cvm,
#                 C1=1,
#                 beta=2,
#             )
#     return outputMatrix


# Tests
# size = volume
# startSize = 0.005
# endSize = 0.02
# steps = 20

# sizes = np.arange(startSize, endSize, (endSize - startSize) / steps)



# dictToClass = {"model":"luo",
#                 "sigma":0.072,
#                 "epsilonCont":0.7,
#                 "rhoCont":1000,
#                 "rhoDisp":1.225,
#                 "Cvm":1
#                 }

# print(CoalescenceModel(0.01, 0.02, dictToClass).coalescenceRate)
