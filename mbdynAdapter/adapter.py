import numpy as np
from mpi4py import MPI
import precice

class MBDynAdapter:
    def __init__(self, mbdHelper, configFileName = "../precice-config.xml"):
        self.mbd = mbdHelper
        self.interface = precice.Interface("Structure_Solver", 0, 1) # proc no, nprocs
        self.interface.configure(configFileName)
        self.dim = self.interface.get_dimensions()
        nodes = self.mbd.getNodes()

        # if self.dim == 2:
            # self.nodesid = np.where(nodes[:,2] < 1e-6)
            # nodes = nodes[self.nodesid]

        self.nnodes = len(nodes)

        nmeshID = self.interface.get_mesh_id("Structure_Nodes")
        self.nodeVertexIDs = self.nnodes*[0.0]
        self.interface.set_mesh_vertices(nmeshID, self.nnodes, np.ravel(nodes[:,:self.dim]), self.nodeVertexIDs)
        self.displacements = np.array(self.dim*self.nnodes*[0.0])

        # ccs = self.mbd.getCellCenters()
        # self.ncells = len(ccs)
        # cmeshID = self.interface.get_mesh_id("Structure_CellCenters")
        # self.cellVertexIDs = self.ncells*[0.0]
        # self.interface.set_mesh_vertices(cmeshID, self.ncells, np.ravel(ccs[:,:self.dim]), self.cellVertexIDs)
        self.force = np.array(self.dim*self.nnodes*[0.0])

        self.displacementsID = self.interface.get_data_id("DisplacementDelta", nmeshID)
        self.forceID = self.interface.get_data_id("Force", nmeshID)

        self.dt = self.interface.initialize()
        self.mbd.controlDict["timeStep"] = self.dt
        self.mbd.initializeMBDyn()


        if (self.interface.is_action_required(precice.action_write_initial_data())):
            self.interface.write_block_vector_data(self.displacementsID, self.nnodes, self.nodeVertexIDs, self.displacements)
            self.interface.fulfilled_action(precice.action_write_initial_data())

        self.interface.initialize_data()

        if (self.interface.is_read_data_available()):
            self.interface.read_block_vector_data(self.forceID, self.nnodes, self.nodeVertexIDs, self.force)

    def runPreCICE(self):
        iteration = 0
        previousDisplacements = self.mbd.getDisplacements()
        while (self.interface.is_coupling_ongoing()):
            if (self.interface.is_action_required(precice.action_write_iteration_checkpoint())):
                self.interface.fulfilled_action(precice.action_write_iteration_checkpoint())

            if self.dim == 2:
                f = np.zeros((self.nnodes,3))
                f[:,:self.dim] = np.reshape(self.force,(-1,self.dim))
            else:
                f = np.reshape(self.force,(-1,3))
            self.mbd.setForces(f)

            if self.mbd.solve(False):
                break
            displacements = self.mbd.getDisplacements()
            relDisplacements = displacements - previousDisplacements

            self.interface.write_block_vector_data(self.displacementsID, self.nnodes, self.nodeVertexIDs, np.ravel(relDisplacements))
            self.interface.advance(self.dt)
            self.interface.read_block_vector_data(self.forceID, self.nnodes, self.nodeVertexIDs, self.force)

            if (self.interface.is_action_required(precice.action_read_iteration_checkpoint())): # i.e. not yet converged
                self.interface.fulfilled_action(precice.action_read_iteration_checkpoint())
            else:
                previousDisplacements = displacements.copy()
                iteration += 1

                if self.mbd.solve(True):
                    break
                if iteration % self.mbd.controlDict['output frequency']  == 0:
                    self.mbd.writeVTK(iteration)
        self.mbd.finalize()
