import numpy as np
from mpi4py import MPI
import PySolverInterface as precice

class MBDynAdapter:
    def __init__(self, mbdHelper, configFileName = "preCICE.xml"):
        self.mbd = mbdHelper
        self.interface = precice.PySolverInterface("Structure_Solver", 0, 1) # proc no, nprocs
        self.interface.configure(configFileName)
        self.dim = self.interface.getDimensions()
        nodes = self.mbd.getNodes()

        if self.dim == 2:
            self.nodesid = np.where(nodes[:,2] < 1e-6)
            nodes = nodes[self.nodesid]

        self.nnodes = len(nodes)

        nmeshID = self.interface.getMeshID("Structure_Nodes")
        self.nodeVertexIDs = self.nnodes*[0.0]
        self.interface.setMeshVertices(nmeshID, self.nnodes, np.ravel(nodes[:,:self.dim]), self.nodeVertexIDs)
        self.displacements = np.array(self.dim*self.nnodes*[0.0])

        ccs = self.mbd.getCellCenters()
        self.ncells = len(ccs)
        cmeshID = self.interface.getMeshID("Structure_CellCenters")
        self.cellVertexIDs = self.ncells*[0.0]
        self.interface.setMeshVertices(cmeshID, self.ncells, np.ravel(ccs[:,:self.dim]), self.cellVertexIDs)
        self.stresses = np.array(self.dim*self.ncells*[0.0])

        self.displacementsID = self.interface.getDataID("Displacements", nmeshID)
        self.stressesID = self.interface.getDataID("Stresses", cmeshID)

        self.dt = self.interface.initialize()
        self.mbd.controlDict["timeStep"] = self.dt
        self.mbd.initializeMBDyn()


        if (self.interface.isActionRequired(precice.PyActionWriteInitialData())):
            self.interface.writeBlockVectorData(self.displacementsID, self.nnodes, self.nodeVertexIDs, self.displacements)
            self.interface.fulfilledAction(precice.PyActionWriteInitialData())

        self.interface.initializeData();

        if (self.interface.isReadDataAvailable()):
            self.interface.readBlockVectorData(self.stressesID, self.ncells, self.cellVertexIDs, self.stresses)

    def runPreCICE(self):
        iteration = 0
        previousDisplacements = self.mbd.getDisplacements()
        while (self.interface.isCouplingOngoing()):
            if (self.interface.isActionRequired(precice.PyActionWriteIterationCheckpoint())):
                self.interface.fulfilledAction(precice.PyActionWriteIterationCheckpoint())

            if self.dim == 2:
                s = np.zeros((self.ncells,3))
                s[:,:self.dim] = np.reshape(self.stresses,(-1,self.dim))
            else:
                s = np.reshape(self.stresses,(-1,3))
            self.mbd.setLoads(0,s)

            if self.mbd.solve(False):
                break
            displacements = self.mbd.getDisplacements()
            relDisplacements = displacements - previousDisplacements
            if self.dim == 2:
                relDisplacements = relDisplacements[self.nodesid][:,:self.dim]

            self.interface.writeBlockVectorData(self.displacementsID, self.nnodes, self.nodeVertexIDs, np.ravel(relDisplacements))
            self.interface.advance(self.dt)
            self.interface.readBlockVectorData(self.stressesID, self.ncells, self.cellVertexIDs, self.stresses)

            if (self.interface.isActionRequired(precice.PyActionReadIterationCheckpoint())): # i.e. not yet converged
                self.interface.fulfilledAction(precice.PyActionReadIterationCheckpoint())
            else:
                previousDisplacements = displacements.copy()
                iteration += 1

                if self.mbd.solve(True):
                    break
                if iteration % self.mbd.controlDict['output frequency']  == 0:
                    self.mbd.writeVTK(iteration)
        mbd.finalize()
