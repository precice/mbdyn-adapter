import itertools
import numpy as np
import os
from subprocess import Popen, PIPE
from mbc_py_interface import mbcNodal

class Mesh:
    pass

class MBDynHelper:
    def __init__(self):
        self.initialized = False
        self.loadsChanged = False

    def initializeMBDyn(self):
        self.writeMBDynScript()
        self.membraneForces = []
        self.nodalForces = []
        path = '{}.sock'.format(self.caseName)
        host = ''
        port = 0
        timeout = -1 # forever
        verbose = 0
        data_and_next = 1
        refnode = 0 # no reference node
        labels = 0
        rot = 0 # orientation vector
        accels = 0
        self.structurelogfile = open("log.mbdyn",'w')
        self.process = Popen(['mbdyn','-f', self.mbdscript], stdout=self.structurelogfile, stderr=self.structurelogfile)
        self.process.stdin = ''
        # self.process = Popen(['mbdyn','-f', self.mbdscript, '-o', self.caseName])
        self.nodal = mbcNodal(path, host, port, timeout, verbose, data_and_next, refnode, len(self.mesh.nodes), labels, rot, accels);
        self.nodal.negotiate()
        self.nodal.recv()
        self.initialized = True
        print('MBDyn initialized\ncontrolDict {} \nmaterialDict {}'.format(self.controlDict, self.materialDict))

    def readMsh(self, fileName):
        self.caseName = os.path.splitext(os.path.basename(fileName))[0]
        names = []
        nodes = None
        membranes = []
        edges = []
        with open(fileName) as fin:
            readNames = False
            readNodes = False
            readElements = False
            nlines = 0
            for line in fin:
                if '$PhysicalNames' in line:
                    readNames = True
                elif '$Nodes' in line:
                    readNodes = True
                elif '$Elements' in line:
                    readElements = True
                elif readNames:
                    nlines = int(line)
                    lines_gen = itertools.islice(fin,0,nlines)
                    for l in lines_gen:
                        cells = l.split()
                        names.append((int(cells[0]), int(cells[1]), cells[2].replace('"','')))
                    readNames = False
                elif readNodes:
                    nlines = int(line)
                    nodes = np.genfromtxt(itertools.islice(fin,0,nlines),dtype=float)
                    readNodes = False
                elif readElements:
                    nlines = int(line)
                    lines_gen = itertools.islice(fin,0,nlines)
                    for e in lines_gen:
                        cells = map(int,e.split())
                        if cells[1] == 1:
                            edges.append(cells)
                        elif cells[1] == 3:
                            membranes.append(cells)
                        else:
                            print('Skipped element, type {}'.format(cells[1]))
                    edges = np.array(edges)
                    membranes = np.array(membranes)

                    self.mesh = Mesh()
                    self.mesh.nodes = np.array(nodes)[:,1:]
                    self.mesh.edges = edges[:,-2:] - 1
                    self.mesh.edgeNames = edges[:,3]
                    self.mesh.membranes = membranes[:,-4:] - 1
                    self.mesh.membraneNames = membranes[:,3]
                    self.mesh.names = names
                    readElements = False
        self.indexing = 10**(int(np.log10(len(self.mesh.nodes))+1))
        self.nnodes = len(self.mesh.nodes)
        print('File: {}'.format(fileName))
        print('Named regions')
        print('{0:<15} {1:<10} {2:<10}'.format('Name','Beams','Membranes'))
        for n in self.mesh.names:
            nedges = len(self.mesh.edges[np.where(self.mesh.edgeNames == n[1])])
            nmembranes = len(self.mesh.membranes[np.where(self.mesh.membraneNames == n[1])])
            print('{:<15} {:<10} {:<10}'.format(n[2],nedges,nmembranes))
        print('{:-^36}'.format(''))
        print('{:<15} {:<10} {:<10}'.format('Sum',len(self.mesh.edges),len(self.mesh.membranes)))
        self.materialDict = {'C':0, 'E':0, 'nu':0, 't':0, 'rho:':0}
        self.initialized = False
        self.preStress = []

    def getAreas(self, projection = False):
        if self.initialized:
            n = np.reshape(self.nodal.n_x, (-1,3))
        else:
            n = self.mesh.nodes
        tris = np.append(self.mesh.membranes[:,:3], \
                         self.mesh.membranes[:,[0,2,3]],axis = 0)
        proj = np.cross((n[tris[:,1]] - n[tris[:,0]]),\
                        (n[tris[:,2]] - n[tris[:,0]]))
        areas = 0.5 * np.linalg.norm(proj,axis = 1)
        if projection:
            return areas, proj
        else:
            return areas

    def getCellCenters(self):
        if self.initialized:
            n = np.reshape(self.nodal.n_x,(-1,3))
        else:
            n = self.mesh.nodes
        tris = np.append(self.mesh.membranes[:,:3], self.mesh.membranes[:,[0,2,3]],axis = 0)
        cellCenters = np.sum(n[tris], axis=1) / 3.0
        return cellCenters

    def getNodes(self):
        return self.mesh.nodes

    def getDisplacements(self):
        return np.reshape(self.nodal.n_x, (-1,3)) - self.mesh.nodes

    def setForces(self, force):
        self.force = force
        self.nodal.n_f[:] = np.ravel(force)

    def setDistributedLoads(self, pressure = 0, stresses = 0):
        self.pressure = pressure
        self.stresses = stresses
        self.loadsChanged = True

    def calcLoads(self):
        pressure = self.pressure
        stresses = self.stresses

        forces = np.zeros(self.mesh.nodes.shape)
        areas, proj = self.getAreas(True)

        triforces = proj * pressure / 2.0 + stresses * areas[:,None]
        memforces = triforces[:len(self.mesh.membranes)] + triforces[len(self.mesh.membranes):]

        # Preferrably no looping
        for i,m in enumerate(self.mesh.membranes):
            forces[m] += memforces[i] / 4.0

        self.membraneForces = memforces
        self.nodalForces = forces

        # Must be the same pointer
        self.nodal.n_f[:] = np.ravel(forces)
        return forces

    def solve(self, converged = True):
        # self.calcLoads()
        stop = self.nodal.send(converged)
        if stop:
            self.writeVTK('final')
            return True
        stop = self.nodal.recv()
        if stop:
            self.writeVTK('final')
            return True
        return False

    def solveStatic(self, tolerance = 1e-6, maxIterations = 10000, write=True):
        previousX = 0
        for i in range(maxIterations):
            self.calcLoads()
            if self.solve(True):
                return
            currentX = self.getDisplacements()
            twoNormX = np.linalg.norm(currentX - previousX)
            previousX = currentX
            print('Finished iteration: {}/{}, displacement two-norm diff: {}/{}'.format(i,maxIterations, twoNormX, tolerance))
            if twoNormX < tolerance:
                print('Converged in {}/{} iterations'.format(i, maxIterations))
                if write:
                    self.writeVTK('final')
                return
        print('No convergence in {} iterations'.format(maxIterations))
        return True

    def finalize(self):
        self.nodal.destroy()

    def writeVTK(self, extension=False):
        nnodes = len(self.mesh.nodes)
        nmembranes = len(self.mesh.membranes)
        nodesstr = '\n'.join(['{} {} {}'.format(*n) for n in self.mesh.nodes])
        membranesstr = '\n'.join(['4 {} {} {} {}'.format(*m) for m in self.mesh.membranes])
        typestr = '\n'.join(np.array(9*np.ones(len(self.mesh.membranes),dtype=int),dtype=str))
        displacements = self.getDisplacements()
        displacementsstr = '\n'.join(['{} {} {}'.format(*d) for d in (displacements)])
        stressesstr = ''
        if os.path.isfile('{}.pla'.format(self.caseName)):
            memStresses = []
            with open('{}.pla'.format(self.caseName)) as fin:
                lines = fin.readlines()
                for l in lines[:-len(self.mesh.membranes) - 1:-1]:
                    memStresses.append(l.split()[1:])
                memStresses = np.reshape(np.array(memStresses[::-1],dtype=float),(-1,4,3))
                stresses = np.zeros(self.mesh.nodes.shape)
                for i,m in enumerate(self.mesh.membranes):
                    stresses[m] += memStresses[i]
            stressesstr = '\n'.join(['{} {} {}'.format(*s) for s in (stresses)])
        if len(self.nodalForces):
            nodalForcesstr = '\n'.join(['{} {} {}'.format(*f) for f in (self.nodalForces)])
            forcesstr = '\n'.join(['{} {} {}'.format(*f) for f in (self.membraneForces)])
        if extension:
            fileName = '{}_{:0>5}.vtk'.format(self.caseName,extension)
        else:
            fileName = '{}.vtk'.format(self.caseName)
        print("Writing output: {}".format(fileName))
        with open(fileName,'w') as fout:
            fout.write('''# vtk DataFile Version 2.0
{}
ASCII
DATASET UNSTRUCTURED_GRID
POINTS {nnodes} float
{nodesstr}
CELLS {nmembranes} {}
{membranesstr}
CELL_TYPES {nmembranes}
{typestr}

POINT_DATA {nnodes}
VECTORS displacements float
{displacementsstr}
'''.format(self.caseName, 5*nmembranes, **locals()))
            if stressesstr:
                fout.write('''
VECTORS stresses float
{stressesstr}
'''.format(**locals()))
            if len(self.nodalForces):
                fout.write('''
VECTORS forces float
{nodalForcesstr}

CELL_DATA {nmembranes}
VECTORS forces float
{forcesstr}
'''.format(**locals()))


    def writeMBDynScript(self, fileName = ''):
        if not fileName:
            self.mbdscript = self.caseName + '.mbd'
        else:
            self.mbdscript = filename
        nnodes = self.nnodes
        nmembranes = len(self.mesh.membranes)
        njoints = nnodes
        materialDict = {'C':0, 'E':0, 'nu':0, 't':0, 'rho':0}
        materialDict.update(self.materialDict)
        if 'rho' in self.materialDict and self.materialDict['rho']:
            nmasses = len(self.mesh.nodes)
            mass = self.MBMass()
        else:
            nmasses = 0
            mass = '' 
        if 'C' in self.materialDict and self.materialDict['C']:
            damping = self.MBDamping()
            nnodes += 1
            njoints += (nnodes + 1)
        else:
            damping = ''
        
        with open(self.mbdscript, 'w') as fout:
            fout.write('''
begin: data;
    problem: initial value;
end: data;

begin: initial value;
    initial time: {initialTime};
    final time: {finalTime};
    time step: {timeStep};
    method: ms, 0.6;
    tolerance: 1e-6;
    max iterations: 300;
    linear solver: umfpack;
    output: iterations; 
    threads: disable;
    modify residual test;
end: initial value;

begin: control data;
    structural nodes: {nnodes};
    beams: 0;
    plates: {nmembranes};
    joints: {njoints};
    rigid bodies: {nmasses};
    forces: 1;
    default output: none, plates;
    output frequency: {output frequency};
end: control data;
'''.format(nnodes = nnodes, njoints = njoints, nmembranes = nmembranes, nmasses=nmasses, **self.controlDict))

            fout.write('''
set: real C   = {C};
set: real Et  = {E};
set: real nut = {nu};
set: real t   = {t};
set: real rho = {rho};
'''.format(**materialDict))

            fout.write('''
begin: nodes;
{}
end: nodes;

begin: elements;
{}
{}
{}
{}
{}
end: elements;
'''.format(self.MBNodes(), self.MBJoints(), damping, mass, self.MBMembranes(), self.MBForces()))

    def setPreStress(self, preStress):
        self.preStress = preStress
        


    def MBNodes(self):
        nodesstr = "structural: {}, dynamic, {}, {}, {}, eye, null, null;"
        string = '\n'.join([nodesstr.format(i,n[0],n[1],n[2]) for i, n in enumerate(self.mesh.nodes)])
        if 'C' in self.materialDict and self.materialDict['C']:
            string += '\nstructural: {}, static, 0.0, 0.0, -1.0, eye, null, null, output, no;'.format(len(self.mesh.nodes))

        return string

    def MBDamping(self):
        index = 2 * self.indexing 
        clampnode = len(self.mesh.nodes)
        string = 'joint: {}, clamp, {}, node, node;\n'.format(index, clampnode)
        index += 1
        # string += '\n'.join(['joint: {}, deformable displacement joint, {}, null, {}, null, linear viscoelastic isotropic, K, C;'.format(index+i, clampnode, i) for i in range(len(self.nodes))])
        string += '\n'.join(['joint: {}, deformable displacement joint, {}, null, {}, null, linear viscous, C;'.format(index+i, clampnode, i) for i in range(len(self.mesh.nodes))])
        string += '\n\n'
        return string

    def MBJoints(self):
        index = 3 * self.indexing
        jointsstr = "joint: {}, total pin joint, {}, position, null, position orientation, eye, rotation orientation, eye,\nposition, {}, {}, {}, position orientation, eye, rotation orientation, eye,\nposition constraint, {}, {}, {}, null, orientation constraint, active, active, active, null;\n"

        mbjoints = []

        dofs = np.zeros(self.mesh.nodes.shape, dtype=bool)

        for name in self.mesh.names:
            if 'fix' in name[2]:
                iedges = self.mesh.edges[np.where(self.mesh.edgeNames == name[1])]
                imembranes = self.mesh.membranes[np.where(self.mesh.membraneNames == name[1])]
                idofs = np.unique(np.append(np.ravel(iedges),np.ravel(imembranes)))
                if 'All' in name[2]:
                    dofs[idofs] = True
                if 'X' in name[2]:
                    dofs[idofs, 0] = True
                if 'Y' in name[2]:
                    dofs[idofs, 1] = True
                if 'Z' in name[2]:
                    dofs[idofs, 2] = True
                if 'XY' in name[2]:
                    dofs[idofs, 0] = True
                    dofs[idofs, 1] = True
                if 'XZ' in name[2]:
                    dofs[idofs, 0] = True
                    dofs[idofs, 2] = True
                if 'YZ' in name[2]:
                    dofs[idofs, 1] = True
                    dofs[idofs, 2] = True
        for i,d in enumerate(dofs):
            dof = ['inactive','inactive','inactive']
            n = self.mesh.nodes[i]
            if d[0]:
                dof[0] = 'active'
            if d[1]:
                dof[1] = 'active'
            if d[2]:
                dof[2] = 'active'
            mbjoints.append(jointsstr.format(index + i, i, n[0],n[1],n[2], *dof))
        string = '\n'.join(mbjoints)
        return string

    def MBMass(self):
        index  = 6 * self.indexing
        massstr = "body: {}, {}, {}*rho*t, null, diag, 0.0, 0.0, 0.0;"

        masses = np.zeros(len(self.mesh.nodes))

        areas = self.getAreas()

        memforces = areas[:len(self.mesh.membranes)] + areas[len(self.mesh.membranes):]

        # Preferrably no looping
        for i,m in enumerate(self.mesh.membranes):
            masses[m] += areas[i] / 4.0

        string = '\n'.join([massstr.format(index + i, i, m) for i, m in enumerate(masses)])
        return string

        
    def MBMembranes(self):
        index = 4 * self.indexing
        if len(self.preStress):
            membranesstr = "membrane4eas: {}, {}, {}, {}, {}, isotropic, E, Et, nu, nut, thickness, t, prestress, {}, {}, {};"
            string = '\n'.join([membranesstr.format(index + i, m[0], m[1], m[2], m[3], *ps) for i, (m, ps) in enumerate(zip(self.mesh.membranes, self.preStress))])
        else:
            membranesstr = "membrane4eas: {}, {}, {}, {}, {}, isotropic, E, Et, nu, nut, thickness, t;"
            string = '\n'.join([membranesstr.format(index + i, m[0], m[1], m[2], m[3]) for i, m in enumerate(self.mesh.membranes)])
        return string 

    def MBForces(self):
        index = 5 * self.indexing
        nodesstr = ','.join(map(str,range(len(self.mesh.nodes))))
        string = '''
force: {}, external structural, socket, create, yes, path, "{}.sock", coupling, tight, orientation, none, accelerations, no, {}, {};'''.format(index, self.caseName, len(self.mesh.nodes), nodesstr)
        return string

if __name__ == "__main__":
    pass
