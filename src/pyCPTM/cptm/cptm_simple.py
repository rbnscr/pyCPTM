import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import bsr_matrix
from scipy.sparse import csc_matrix
from scipy.sparse.csgraph import connected_components

from scipy.integrate import solve_ivp
from pyCPTM.cptm.basic_models import tracer_model_sparse
import pandas as pd

from pyCPTM.io import load_mesh_from_case
from pyCPTM.io import load_phi
from pyCPTM.io import write_steady_scalar_to_case

from pyCPTM.utilities import sort_unsort

import time

class CPTMSimple:
    """Baseline compartment modeling class with all the needed meshing and compartmentalization methods
    """
    def __init__(
        self,
        internalMesh,
        internalPhi,
        cellVolume: np.ndarray,
        cellCoordinates: np.ndarray,
        cellData=[],
        cptIdx=[],
    ):
        self.internalMesh = np.copy(internalMesh)
        self.internalPhi = np.copy(internalPhi[: len(internalMesh)])
        self.cellVolume = cellVolume
        self.cellCoordinates = cellCoordinates
        self.nCell = len(cellVolume)

        assert isinstance(cellData, pd.DataFrame), "cellData must be a pandas.DataFrame"
        self.cellData = cellData
        self.__cptIdx = np.copy(cptIdx)
        if len(cptIdx):
            self.calc_n_cpt()


        if not hasattr(self, "cptIdxCache"):
            self.init_idxCache()

    @property
    def cptIdx(self):
        return self.__cptIdx

    @cptIdx.setter
    def cptIdx(self, cptIdx):
        if len(cptIdx) == self.nCell:
            self.__cptIdx = cptIdx

    def calc_cpt_volume(self, phaseAlpha: np.ndarray = []):
        """Calculates the volume of all compartments.

        Args:
            phaseAlpha (np.ndarray, optional): Volume-portion of a given phase per cell. Only useful with multiphase approaches, which are not yet implemented. Defaults to [], which sets phaseAlpha to 1 for each cell.
        """
        if self.cptIdx.size == 0:
            # Returns if there are no compartments defined yet.
            print("Compartment indices have to be defined.")
        else:
            self.calc_n_cpt()
            if phaseAlpha == []:
                phaseAlpha = np.ones(self.nCell)
            cptVolume: np.ndarray = np.zeros(self.nCpt) # Initialize volume array

            # Calculation of the compartment volumes
            for i in range(self.nCpt):
                cptVolume[i] = self.cellVolume[self.cptIdx == i].T.dot(phaseAlpha[self.cptIdx == i])
            self.cptVolume = cptVolume
            print(f"Volume(s) of {self.nCpt} compartment(s) calculated.")

    def calc_n_cpt(self):
        """Calculate the number of compartments.
        """
        self.nCpt = len(np.unique(self.cptIdx))
        print(f"Number of compartments: {self.nCpt}")
        # If the compartment indices don't line up with the compartment number, the compartment are renumbered.
        if max(self.cptIdx + 1) != self.nCpt:
            self.renumber()

    def renumber(self):
        """Renumber the compartment indices. Followed by a recalculation of the number of compartments.
        """
        newIndex = 0
        for i in np.unique(self.cptIdx):
            self.cptIdx[self.cptIdx == i] = newIndex
            newIndex += 1
        # self.nCpt = len(np.unique(self.cptIdx))
        self.calc_n_cpt()
        print(f"{self.nCpt} compartments renumbered.")

    def projection(self):
        """Spatially projects the compartments. Cells, that belong to the same compartment, but are not continuously connected, are assigned new compartment indices. This method does not work with dynamic mesh approach.
        """
        print(f"Number of compartments before projection: {self.nCpt}")

        # New compartment indice array is initialized. All unassigned cells have the value -1. The indexing starts at 0.
        idxCpt = np.zeros(len(self.cellCoordinates), dtype=np.int32) 
        idxCpt = idxCpt - 1 
        newIndex = 0
        # check = np.zeros(len(self.internalMesh), dtype=np.int32)

        # Looping over all present compartment indices. 
        # The projected compartments must have at least the same amount of compartments as before the projection.
        # All unprojected compartments are called _cluster_ in the following loop.
        for cluster in range(self.nCpt):
            internalMesh = np.copy(self.internalMesh)
            cellsInCluster = np.where(self.cptIdx == cluster)

            # All Faces, that are within the cell cluster are labeled as TRUE
            facesInCluster = np.logical_and(
                np.isin(internalMesh[:, 0], cellsInCluster),
                np.isin(internalMesh[:, 1], cellsInCluster),
            )

            # "New" Mesh for the creation of a sparse matrix
            internalMesh = internalMesh[facesInCluster]

            # Create graph with only ones as weights
            graph = bsr_matrix(
                (np.ones(len(internalMesh)), (internalMesh[:, 1], internalMesh[:, 0])),
                shape=(self.nCell, self.nCell),
            )

            # Analyze the connected components of the sparse graph _graph_
            _nComponents, labels = connected_components(
                csgraph=graph, directed=False, return_labels=True
            )

            # Trim and renumber labels
            labels = labels[cellsInCluster]
            labels_ = np.copy(labels)
            for i in np.unique(labels):
                labels[labels_ == i] = newIndex
                newIndex += 1

            # Integrate into idxCpt
            idxCpt[cellsInCluster] = labels
        self.cptIdx = idxCpt
        self.calc_n_cpt()
        self.calc_cpt_volume()

    def planar_slices(self, nSlices, sliceStart, sliceEnd, dim):

        # ISSUE: would be nicer, if all slice-methods would be cached and could then be consolidated at the end.
        # Initialize new compartment indices
        idxCpt = np.zeros(len(self.cellCoordinates), dtype=np.int32)
        idxCpt = idxCpt - 1
        cellCoordinates = np.copy(self.cellCoordinates)

        # Exclude outliers
        detectedOutliers = np.logical_or(
            cellCoordinates[:, dim] < sliceStart, cellCoordinates[:, dim] >= sliceEnd
        )

        cellCoordinates[detectedOutliers, :] = np.nan

        # Thickness of each slice
        sliceThickness = abs(sliceStart - sliceEnd) / nSlices

        # Slicing
        for i in range(nSlices):
            cellsInSlice = np.logical_and(
                cellCoordinates[:, dim] >= sliceThickness * i + sliceStart,
                cellCoordinates[:, dim] < sliceThickness * (i + 1) + sliceStart,
            )
            idxCpt[cellsInSlice] = i

            idxCpt[idxCpt == -1] = max(idxCpt)
            idxCpt[detectedOutliers] = -1
            self.cptIdx = idxCpt

    def plot_idx(self, compartment=[]):
        """Plotting indices of a given compartment(s).

        Args:
            compartment (list, optional): Compartment indice to be plotted. Defaults to [].
        """
        if self.cptIdx.size == 0:
            print("Compartment indices have to be defined.")
        elif compartment == []:
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")
            im = ax.scatter(
                self.cellCoordinates[:, 0],
                self.cellCoordinates[:, 1],
                self.cellCoordinates[:, 2],
                c=self.cptIdx,
            )
            fig.colorbar(im, ax=ax)
        else:
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")
            im = ax.scatter(
                self.cellCoordinates[self.cptIdx == compartment, 0],
                self.cellCoordinates[self.cptIdx == compartment, 1],
                self.cellCoordinates[self.cptIdx == compartment, 2],
            )

    def conv_exchange_flow(self):
        """Calculated the convective transitioning matrix based on Delafosse (2014) (http://dx.doi.org/10.1016/j.ces.2013.11.033)
        """

        if not hasattr(self, "nCpt"):
            self.calc_n_cpt()

        # Flip internalMesh, so that all flows are positive
        _internalMesh = np.copy(self.internalMesh)
        _internalMesh[:, [0, 1]] = _internalMesh[:, [1, 0]]

        internalMeshFlipped = np.copy(self.internalMesh)
        internalPhiFlipped = np.copy(self.internalPhi)

        swapOwnerNeighbour = self.internalPhi < 0
        internalMeshFlipped[swapOwnerNeighbour, :] = _internalMesh[
            swapOwnerNeighbour, :
        ]
        internalPhiFlipped = abs(self.internalPhi)

        # Sort owner and neighbour
        ownerSorted, ownerUnsortIndex = sort_unsort(internalMeshFlipped[:, 0])
        neighbourSorted, neighbourUnsortIndex = sort_unsort(internalMeshFlipped[:, 1])

        # Setup for the _first_ loop
        firstEntryOwner = 0
        firstEntryNeighbour = 0
        # return ownerSorted

        # Adding 25 is kind of a hack to (drastiacally) speed up the algorithm. It is assumed that no cell has more than 25 faces. 
        # Calculating the max number of faces per cell however is also slow.
        for i in range(self.nCell + 1):
            try:
                lastEntryOwner = (
                    np.argwhere(
                        ownerSorted[
                            firstEntryOwner : min(
                                len(ownerSorted), firstEntryOwner + 25
                            )
                        ]
                        == i
                    )[-1][0]
                    + firstEntryOwner
                    + 1  # Add 1 for indexing
                )
                ownerSorted[firstEntryOwner:lastEntryOwner] = self.cptIdx[i]
                firstEntryOwner = lastEntryOwner
            except IndexError:
                lastEntryOwner = firstEntryOwner
            try:
                lastEntryNeighbour = (
                    np.argwhere(
                        neighbourSorted[
                            firstEntryNeighbour : min(
                                len(neighbourSorted), firstEntryNeighbour + 25
                            )
                        ]
                        == i
                    )[-1][0]
                    + firstEntryNeighbour
                    + 1  # Add 1 for indexing
                )
                neighbourSorted[firstEntryNeighbour:lastEntryNeighbour] = self.cptIdx[i]
                firstEntryNeighbour = lastEntryNeighbour
            except IndexError:
                lastEntryNeighbour = firstEntryNeighbour
        # Unsort indices
        ownerUnsorted = ownerSorted[ownerUnsortIndex]
        neighbourUnsorted = neighbourSorted[neighbourUnsortIndex]

        compartmentMesh = np.array([ownerUnsorted, neighbourUnsorted]).T

        # Ignore diagonal
        # bool, where compartments differ
        nonDiagonalEntries = compartmentMesh[:, 0] != compartmentMesh[:, 1]
        # Compartment Mesh and Phi, where compartments differ
        nonDiagonalCompartmentMesh = compartmentMesh[nonDiagonalEntries, :]
        nonDiagonalInternalPhi = internalPhiFlipped[nonDiagonalEntries]

        # Unique compartment pairs
        uniqueCompartmentPairs = np.unique(compartmentMesh, axis=0)
        # nonDiagonalUniqueEntries = (
        #     uniqueCompartmentPairs[:, 0] != uniqueCompartmentPairs[:, 1]
        # )
        # nonDiagonalUniqueCompartmentPairs = uniqueCompartmentPairs[nonDiagonalUniqueEntries]
        nonDiagonalUniqueCompartmentPairs = uniqueCompartmentPairs

        # initialize new phi for 'uniqueCompartmentPairs'
        nonDiagonalCompartmentPhi = np.zeros(len(nonDiagonalUniqueCompartmentPairs))

        # Sum up all phi for each pair
        for i in range(len(nonDiagonalUniqueCompartmentPairs)):
            entriesForCompartmentPair = np.logical_and(
                nonDiagonalCompartmentMesh[:, 0]
                == nonDiagonalUniqueCompartmentPairs[i, 0],
                nonDiagonalCompartmentMesh[:, 1]
                == nonDiagonalUniqueCompartmentPairs[i, 1],
            )

            nonDiagonalCompartmentPhi[i] = sum(
                nonDiagonalInternalPhi[entriesForCompartmentPair]
            )
        # Create sparse matrix
        exchangeFlowGraph: csc_matrix = csc_matrix(
            (
                nonDiagonalCompartmentPhi,
                (
                    nonDiagonalUniqueCompartmentPairs[:, 1],
                    nonDiagonalUniqueCompartmentPairs[:, 0],
                ),
            ),
            shape=(self.nCpt, self.nCpt),
        )

        # Diagonal with negative values
        graphDiag = -exchangeFlowGraph.sum(1)[:]
        graphDiag = np.squeeze(np.asarray(graphDiag))
        exchangeFlowGraph.setdiag(graphDiag)

        self.exchangeFlowGraph = exchangeFlowGraph
        # return exchangeFlowGraph

    def init_tracer(self, compartments: list, tracerAmount: float):
        """Initializes conservative tracer in compartments

        Args:
            compartments (list): List of compartments that the tracer concentration shall be initialized
            tracerAmount (float): Tracer concentration that is initialized
        """
        initialConcentration = np.zeros(self.nCpt)
        initialConcentration[compartments] = tracerAmount
        self.initialTracerConcentration = np.asarray(initialConcentration)

    def init_timespan(self, startTime: float, endTime: float, timeStep: float):
        """Initialize the timespan of the virtual tracer simulation. Wrapper for numpy.linspace.

        Args:
            startTime (float): Start time
            endTime (float): endTime
            timeStep (float): timeStep
        """
        self.t = np.linspace(startTime, endTime, timeStep)

    def solve_tracer(self):
        """Starts and solves the virtual tracer simulation.
        """
        if not (
            hasattr(self, "cptVolume")
            and hasattr(self, "exchangeFlowGraph")
            and hasattr(self, "initialTracerConcentration")
            and hasattr(self, "t")
        ):
            print(f"One or more attribute(s) is missing: \ncptVolume {hasattr(self, "cptVolume")}, \nexchangeFlowGraph {hasattr(self, "exchangeFlowGraph")}, \ninitialTracerConcentration {hasattr(self, "initialTracerConcentration")}, \nt {hasattr(self, "t")}\nExiting program.")
            sys.exit(0)
        else:

            exchangeFlowGraph = self.exchangeFlowGraph.multiply(1 / self.cptVolume)

            self.sol = solve_ivp(
                tracer_model_sparse,
                t_span=(0, max(self.t)),
                y0=self.initialTracerConcentration,
                args=([exchangeFlowGraph]),
                t_eval=self.t,
                method="Radau",
            )
            self.normalize_tracer()

    def normalize_tracer(self):
        # Needs some work
        self.normalizedTracer = (
            self.sol.y.T
            # * self.initialTracerConcentration.dot(self.cptVolume)
            / (self.cptVolume)
            * sum(self.cptVolume)
        )

    def plot_tracer(self, ylim):
        # Passt noch nicht ganz. Müsste nach Masse normalisiert werden
        plt.plot(self.t, self.normalizedTracer)
        plt.ylim(ylim)

    def cell_by_cell(self, varNames, numberSplits, autoConsolidate = True):
        entries = len(varNames)
        idxCpt = np.zeros(self.nCell)
        idxCptCache = np.zeros((self.nCell, entries))
        idxCptCache = idxCptCache - 1
        for name in varNames:
            entryIndex = varNames.index(name)
            entrySplits = numberSplits[entryIndex]

            varToSplit = np.asarray(self.cellData.apply(name))
            minVarToSplit = varToSplit.min()
            maxVarToSplit = varToSplit.max()

            margin = (maxVarToSplit - minVarToSplit) / entrySplits

            for i in range(entrySplits):
                idxCptCache[varToSplit >= minVarToSplit + margin * i, entryIndex] += 1
        uniquePairs = np.unique(idxCptCache, axis=0)

        for pair in range(len(uniquePairs)):
            pairMatchingIndex = np.where(
                (idxCptCache == uniquePairs[pair, :]).all(axis=1)
            )
            idxCpt[pairMatchingIndex] = pair

        print(f"{len(np.unique(idxCpt))} sections found by cell-by-cell algorithm.")

        if autoConsolidate == True:
            self.init_idxCache(1)
            self.set_idxCache(np.copy(idxCpt))
            self.consolidate_idx()
        elif autoConsolidate == False:
            self.set_idxCache(np.copy(idxCpt))

        # self.cptIdx = idxCpt
        # self.renumber()
        # self.projection()

    def calc_coef_var(self):
        # Calculates the volume averaged coefficient of variation
        # Haringa (2018)
        volumeSum = sum(self.cptVolume)
        meanTracer = self.normalizedTracer.dot(self.cptVolume) / volumeSum
        diffTracerMean = self.normalizedTracer.T - meanTracer
        # return diffTracerMean, meanTracer
        divTracer = ((diffTracerMean / meanTracer) ** 2).T
        # return divTracer, self.cptVolume
        numerator = divTracer.dot(self.cptVolume)
        # return numerator, volumeSum
        weightedCov = np.sqrt(numerator / volumeSum)
        self.tracerCov = weightedCov
        if any(np.gradient(self.tracerCov)>0):
            print("Numerical error might have caused \"de-mixing\".")



        idxKukukova = np.argwhere(np.diff(np.sign(self.tracerCov - 0.05))).flatten()
        idxHartmann = np.argwhere(np.diff(np.sign(self.tracerCov - 0.0283))).flatten()
        try:
            print(f"Mixing time [s]:\n {self.t[idxKukukova][0]}\n {self.t[idxHartmann][0]}")
        except IndexError:
            print("No full mixing occured.")

    def plot_coef_var(self):
        plt.plot(self.t, self.tracerCov)
        plt.ylabel("Coefficient of Mixing [-]")
        plt.xlabel("Time [s]")
        plt.yscale("log")
        plt.ylim(bottom = 10e-4)
        plt.grid(which="minor", linestyle="--", linewidth=0.3)
        plt.grid()
        # Intersection points
        plt.plot(
            [min(self.t), max(self.t)],
            [0.05, 0.05],
            color="k",
            linestyle="--",
            linewidth=1,
        )
        plt.plot(
            [min(self.t), max(self.t)],
            [0.0283, 0.0283],
            color="k",
            linestyle="--",
            linewidth=1,
        )

        plt.text(0.1, 0.06, "95% Mixing (Kukuková et al., 2008)", fontsize="small")
        plt.text(0.1, 0.02, "95% Mixing (Hartmann et al., 2006)", fontsize="small")

        idxKukukova = np.argwhere(np.diff(np.sign(self.tracerCov - 0.05))).flatten()
        plt.plot(
            self.t[idxKukukova], self.tracerCov[idxKukukova], "ko", fillstyle="none"
        )

        idxHartmann = np.argwhere(np.diff(np.sign(self.tracerCov - 0.0283))).flatten()
        plt.plot(
            self.t[idxHartmann], self.tracerCov[idxHartmann], "ko", fillstyle="none"
        )

    # @classmethod
    # def instance_from_case(cls, pathToCase: str, dataFileName: str, timeStep):
    #     owner, neighbour, boundary = load_mesh_from_case(
    #         pathToCase=absolutePath, maxFaceNumber=5 * 10 ** 6
    #     )
    #     internalMesh = np.array([owner[: len(neighbour)], neighbour]).T
    #     phi, b = load_phi(pathToCase=absolutePath, chosenTimesteps=['11'], maxFaceNumber=10 ** 7)

    #     CPTMSimple()

    def get_cpt_at_coordinate(self, coordinates):
        if self.cptIdx.size != 0:
            return int(self.cptIdx[
                np.argmin(np.sum((self.cellCoordinates - coordinates) ** 2, axis=1))
            ])
        else:
            print("Compartment indices not defined.")

    def merge_cpt(self, mergeThreshold):
        pass

    def hist_vol(self, nBins=100):
        plt.hist(self.cptVolume, bins=nBins)
        plt.ylabel("Frequency [-]")

    def write_idx_to_case(self, pathToCase):
        write_steady_scalar_to_case(
            pathToCase=pathToCase, scalarName="cptIdx", scalar=self.cptIdx
        )

    def init_idxCache(self, nCache=10):
        cptIdxCache = np.zeros((self.nCell,nCache),dtype=np.int64)
        self.cptIdxCache = np.copy(cptIdxCache)
        self.cptIdxCacheCounter = 0
        # print("Initialization successful.")
        # print(self.cptIdxCache)
        # return self.cptIdxCache

    def consolidate_idx(self, applyProjection = True):


        #############################
        # THIS WHOLE FUNCTION IS SLOW
        #############################


        # This is quite similar to renumber(self)
        print("Start consolidation.")
        print(f"{len(np.unique(self.cptIdxCache, axis=0))} unique combinations found.")
        newIndex = 0
        idxCpt = np.zeros(self.nCell)

        if self.cptIdxCache.shape[1] != 1:

            # start = time.time()

            idxCpt = np.zeros(self.nCell)

            for i in np.unique(self.cptIdxCache, axis=0):
                for j in range(self.nCell):
                    if all(self.cptIdxCache[j,:] == i):
                        idxCpt[j] = newIndex
                newIndex += 1
                print(f"Cluster {i} done.")

            # for j in range(self.nCell):
            #     for i in np.unique(self.cptIdxCache, axis=0):
            #         if all(self.cptIdxCache[j,:] == i):
            #             idxCpt[j] = newIndex
            #             break
            #         newIndex += 1
            #     print(f"Cell {j} done.")
            #     newIndex = 0

            # end = time.time()
            # print(end - start)

        else:
            idxCpt = np.copy(self.cptIdxCache[:,0])
        self.cptIdx = idxCpt
        self.calc_n_cpt()
        # self.renumber() # not necessarily needed
        print("Finished consolidation.")
        if applyProjection == True:
            self.projection()
        # else:
        #     self.calc_cpt_volume()
        self.init_idxCache()



    def set_idxCache(self, cptIdxToSet, setAll = False):
        if len(cptIdxToSet) != self.nCell:
            print(f"cptIdxToSet requires length of {self.nCell}.")
        else:
            if setAll == False:
                # Default option to add layers to the compartment index cache
                self.cptIdxCache[:,self.cptIdxCacheCounter] = np.squeeze(cptIdxToSet)
                self.cptIdxCacheCounter += 1
            elif setAll == True:
                # This is for clearing (clear_idx())
                newIndex = np.max(self.cptIdxCache)
                self.cptIdxCache[cptIdxToSet == 1,:] = newIndex + 1

    def clear_idx(self, region):
        # Clears a region from any compartment indices / Creates a new single compartment at the location
        self.set_idxCache(region, setAll = True)