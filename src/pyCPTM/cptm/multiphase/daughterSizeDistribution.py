# -*- coding: utf-8 -*-

import numpy as np

# openfoam implementation https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/applications/solvers/multiphase/multiphaseEulerFoam/phaseSystems/populationBalanceModel/daughterSizeDistributionModels/LaakkonenDaughterSizeDistribution/LaakkonenDaughterSizeDistribution.C


class DaughterSizeDistribution:
    dsdModels = {
        "laakkonen": "laakkonenDaughterSizeDistribution",
        "keyToCall": "methodToCall",
    }

    def __init__(self, i, k, modelDict):
        self.i = i
        self.k = k
        self.modelDict = modelDict
        modelName = modelDict["model"]
        self.call_model(modelName)

    def call_model(self, arg):
        return getattr(self, self.dsdModels[arg])(self.i, self.k, self.modelDict)

    def antiderivativeLaakkonen(self, xk, v, bndr, rang, c):
        nBubbles = 4 / 3 + c / 3
        first_term = (xk ** (-c - 3)) * ((xk - v) ** c) * (v - xk)
        second_term = (
            (c + 1) * (c + 2) * (c + 3) * v ** 3
            - (c + 1) * (c + 2) * (bndr * (c + 4) - 3 * xk) * v ** 2
            - 2 * v * xk * (c + 1) * (bndr * (c + 4) - 3 * xk)
            - 2 * bndr * c * xk ** 2
            + 6 * xk ** 3
            - 8 * bndr * xk ** 2
        )
        third_term = 2 * rang * (c + 4)
        return nBubbles * first_term * second_term / third_term

    def laakkonenDaughterSizeDistribution(self, i, k, modelDict):

        sizeClasses = modelDict["sizeClasses"]
        cdefault = 18.25
        c = modelDict.get("c")
        if c == None:
            c = cdefault
        # Calculate contribution to sizeGroup i due to breakup in sizeGroup k
        # volume based

        x0 = sizeClasses[0]
        xi = sizeClasses[i] - x0
        xk = sizeClasses[k] - x0

        if i == 0:
            xii = sizeClasses[i + 1] - x0

            if k == 0:
                self.dsd = 1
                return
            self.dsd = self.antiderivativeLaakkonen(
                xk, xi, xii, xii - xi, c
            ) - self.antiderivativeLaakkonen(xk, xii, xii, xii - xi, c)
            return
        elif i == k:
            x = sizeClasses[i - 1] - x0
            self.dsd = self.antiderivativeLaakkonen(
                xk, xi, x, xi - x, c
            ) - self.antiderivativeLaakkonen(xk, x, x, xi - x, c)
            return
        else:
            x = sizeClasses[i - 1] - x0
            xii = sizeClasses[i + 1] - x0
            self.dsd = (
                self.antiderivativeLaakkonen(xk, xi, xii, xii - xi, c)
                - self.antiderivativeLaakkonen(xk, xii, xii, xii - xi, c)
                + self.antiderivativeLaakkonen(xk, xi, x, xi - x, c)
                - self.antiderivativeLaakkonen(xk, x, x, xi - x, c)
            )
            return


# def laakkonenDSDMatrix(sizeClasses, c):
#     outputMatrix = np.zeros((len(sizeClasses), len(sizeClasses)))
#     for i in range(len(sizeClasses) - 1):
#         for j in range(len(sizeClasses)):
#             if i <= j:
#                 outputMatrix[i, j] = laakkonenDaughterSizeDistribution(sizeClasses, i, j, c)
#     return outputMatrix


# Tests
# size = volume
# startSize = 0.0005
# endSize = 0.002
# steps = 20

# sizes = np.arange(startSize, endSize, (endSize - startSize) / steps)

# dsdDict = {"model": "laakkonen", "sizeClasses": sizes}

# print(DaughterSizeDistribution(1, 10, dsdDict).dsd)

# test = laakkonenDSDMatrix(sizes,2)

# print(np.sum(test, axis=0))
