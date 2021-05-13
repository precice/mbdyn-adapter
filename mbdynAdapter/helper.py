#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from subprocess import Popen
import os
from mbc_py_interface import mbcNodal
from bs4 import BeautifulSoup
import precice
import numpy as np
from scipy.spatial.transform import Rotation as R
import logging

# create logger
module_logger = logging.getLogger('adapter.helper')

class MBDynHelper:
    def __init__(self, mesh):
        self.initialized = False
        self.process = None
        self.nodal = None
        self.log_file = None
        self.log_file_path = '../mbdyn.log'
        self.input_file_name = 'shell.mbd'
        self.mesh = mesh
        self.load_changed = False
        self.pressure = 0
        self.stresses = 0
        self.node_forces = 0
        #self.cell_forces = None
        self._debug_samples = [0]

    def _equidistant_samples(self, num_samples=6):
        self._debug_samples.append(0)
        num_nodes = self.mesh.number_of_nodes()
        interval = num_nodes/(num_samples-1)
        for i in range(1, num_samples):
            self._debug_samples.append(int(interval*i)-1)


    def initialize(self, case='shell'):


        self.input_file_name = case + '.mbd'
        self.log_file = open(self.log_file_path, 'w')
        self.process = Popen(['mbdyn', '-f', self.input_file_name],
                             stdout=self.log_file,
                             stderr=self.log_file)
        self.process.stdin = ''

        path = '{name}.sock'.format(name=case)
        module_logger.debug('socket path: %s' % path)
        host = ''
        port = 0
        timeout = -1
        verbose = 1
        data_and_next = 1
        refnode = 0
        nodes = self.mesh.number_of_nodes()
        labels = 0 # 16
        rot = 0 # for rotvec 256, rotmat 512, euler 1024; see mbc.h enum MBCType
        accels = 0
        self.nodal = mbcNodal(path, host, port, timeout, verbose,
                              data_and_next, refnode, nodes, labels,
                              rot, accels)
        self.nodal.negotiate()
        print(self.nodal.recv())
        self.initialized = True

    def finalize(self):
        try:
            self.nodal.destroy()
            self.log_file.close()
        except AttributeError:
            print('Warning: Could not close log file or destroy mbc.')

    #TODO: Testing, bug: forces are saved as stresses
    def write_output_vtk(self, extension=False):
        num_nodes = self.mesh.number_of_nodes()
        num_shells = self.mesh.number_of_shells()

        nodes_str = '\n'.join(
            ['{} {} {}'.format(*n) for n in self.get_nodes()])
        shells_str = '\n'.join(
            ['4 {} {} {} {}'.format(*m) for m in self.mesh.shells])
        type_str = '\n'.join(
            np.array(9*np.ones(len(self.mesh.shells), dtype=int), dtype=str))
        displacements_str = '\n'.join(
            ['{} {} {}'.format(*d) for d in self.get_absolute_displacement()])
        stresses_str = ''

        if os.path.isfile('{}.pla'.format(
                os.path.splitext(self.mesh.name)[0])):
            shell_stresses = []
            with open('{}.pla'.format(
                    os.path.splitext(self.mesh.name)[0])) as fin:
                lines = fin.readlines()
                for line in lines[:-len(self.mesh.shells) - 1:-1]:
                    shell_stresses.append(line.split()[1:])
                shell_stresses = np.reshape(
                    np.array(shell_stresses[::-1], dtype=float), (-1, 4, 3))
                stresses = np.zeros(self.mesh.nodes.shape)
                for i, nodes in enumerate(self.mesh.shells):
                    stresses[nodes] += shell_stresses[i]
            stresses_str = '\n'.join(
                ['{} {} {}'.format(*s) for s in stresses])

        if isinstance(self.node_forces, np.ndarray):
            node_forces_str = '\n'.join(
                ['{} {} {}'.format(*f) for f in self.node_forces])
            # cell_forces_str = '\n'.join(
            #    ['{} {} {}'.format(*f) for f in self.cell_forces])

        if extension:
            if 'init' in str(extension):
                file_name = '{}_{}.vtk'.format(
                    os.path.splitext(self.mesh.name)[0], extension)
            else:
                file_name = '{}_{:0>5}.vtk'.format(
                    os.path.splitext(self.mesh.name)[0], extension)
        else:
            file_name = '{}.vtk'.format(os.path.splitext(self.mesh.name)[0])

        print("Writing output: {}".format(file_name))

        def vtk_header(title, num_points, points,
                       num_cells, cells, cell_types):
            header = '# vtk DataFile Version 2.0\n{name}\nASCII\n'.format(
                name=title)

            geometry = 'DATASET UNSTRUCTURED_GRID\n'
            key_points = 'POINTS {npts} float\n{pts}\n'.format(
                npts=num_points, pts=points)
            key_cells = 'CELLS {ncll} {size}\n{cll}\n'.format(
                ncll=num_cells, size=5*num_cells, cll=cells)
            key_cell_types = 'CELL_TYPES {ncll}\n{ctype}\n'.format(
                ncll=num_cells, ctype=cell_types)

            geometry += key_points + key_cells + key_cell_types

            return header + geometry

        def vtk_vector(data_name, data_str, data_type='float'):
            out_str = '''VECTORS {dname} {dtype}\n{dstr}'''
            return out_str.format(
                dname=data_name, dtype=data_type, dstr=data_str)

        with open(file_name, 'w') as output_file:
            output_file.write(
                vtk_header(self.mesh.name, num_nodes, nodes_str,
                           num_shells, shells_str, type_str))
            output_file.write('POINT_DATA {}\n'.format(num_nodes))
            output_file.write(vtk_vector('displacements',
                                         displacements_str))

            if stresses_str:
                output_file.write(vtk_vector('stresses', stresses_str))
            if isinstance(self.node_forces, np.ndarray):
                output_file.write(vtk_vector('forces', node_forces_str))
                #output_file.write('CELL_DATA {}\n'.format(num_shells))
                #output_file.write(vtk_vector('forces', cell_forces_str))

    def get_absolute_displacement(self):
        return self.get_nodes() - self.mesh.nodes

    def get_nodes(self):
        if self.initialized:
            return np.reshape(self.nodal.n_x, (-1, 3))
        else:
            return self.mesh.nodes

    def get_forces(self):
        return np.reshape(self.nodal.n_f, (-1, 3))

    def set_forces(self, forces):
        self.node_forces = forces
        self.nodal.n_f[:] = np.ravel(forces)

    #TODO: create option for stresses
    def set_pressure(self, pressure):
        self.pressure = float(pressure)
        self.load_changed = True

    #TODO: something about the greyed out part breaks things
    def calc_pressure_forces(self, forces=0, relaxation=1, limiting=10):

        module_logger.debug('rotvec from mbdyn: \n %s' % self.nodal.n_theta)

        shell_normals = self.mesh.calc_shell_normals(n_x=self.get_nodes(), invert=1)
        node_normals_weighted = np.zeros((self.mesh.number_of_nodes(), 3))
        for i, nodes in enumerate(self.mesh.shells):
            node_normals_weighted[nodes] += 0.25 * shell_normals[i]

        pressure_forces = node_normals_weighted * self.pressure

        if not isinstance(limiting, type(None)):
            max_value_pressure = np.max(np.linalg.norm(pressure_forces, axis=1))
            if max_value_pressure > limiting:
                pressure_forces = self.node_forces
            if not isinstance(self.node_forces, (int, float)):
                max_value_fluid = np.max(np.linalg.norm(forces, axis=1))
                if max_value_fluid > limiting:
                    forces = forces / max_value_fluid * 0.3

        if relaxation != 1:
            new_forces = self.node_forces + \
                (pressure_forces + forces - self.node_forces) * relaxation
        else:
            new_forces = pressure_forces + forces

        # self.node_forces = np.multiply(fix_force, self.pressure)
        forces_norm = np.linalg.norm(new_forces, axis=1)
        module_logger.debug(
            'min, max, sum forces after pressure applied:\n{}, {}, {}'.format(
                np.min(forces_norm), np.max(forces_norm),
                np.sum(forces_norm)))
        module_logger.debug(
            'forces after pressure applied sample:\n{}'.format(
                new_forces[self._debug_samples,:]))

        self.set_forces(new_forces)


    def solve(self, converged=False):
        if self.nodal.send(converged):
            self.write_output_vtk('final')
            module_logger.debug('on send')
            return True
        if self.nodal.recv():
            module_logger.debug('on recv')
            self.write_output_vtk('final')
            return True
        return False

    #TODO
    def solve_static(self, tolerance=1e-6, max_iterations=10000,
                     write=True):
        previous_position = 0
        for i in range(max_iterations):
            self.calc_pressure_forces(relaxation=0.3, limiting=20000)
            if self.solve(True):
                return True
            current_position = self.get_absolute_displacement()
            two_norm_diff = np.linalg.norm(
                current_position - previous_position)
            previous_position = current_position
            module_logger.debug('Finished iteration: {}/{}, displacement two-norm diff: {}/{}'.format(
                i, max_iterations, two_norm_diff, tolerance))
            if write:
                    self.write_output_vtk(str(i))
            if two_norm_diff < tolerance and i > 500:
                print('Converged in {}/{} iterations'.format(
                    i, max_iterations))
                if write:
                    self.write_output_vtk('static-final')
                return True
        print('No convergence in {} iterations'.format(max_iterations))
        return False

    def solve_initial(self, tolerance=5e-6, max_iterations=10000,
                     write=True):

        previous_position = 0

        # calculate static force magnitude on each node
        self.calc_pressure_forces()
        node_forces_mag = np.linalg.norm(self.node_forces, axis=1)[:, np.newaxis]

        # calculate node normals
        def node_normals(xyz=self.get_nodes()):
            shell_normals = self.mesh.calc_shell_normals(n_x=xyz, normalize=True)
            normals = np.zeros((self.mesh.number_of_nodes(), 3))
            for i, nodes in enumerate(self.mesh.shells):
                normals[nodes] += shell_normals[i]
            return normalize_vectors(normals)

        # calculate static force in new direction
        def new_force():
            return node_normals() * node_forces_mag

        for i in range(max_iterations):
            if self.solve(True):
                return True

            current_position = self.get_absolute_displacement()

            two_norm_diff = np.linalg.norm(
                current_position - previous_position)

            previous_position = current_position

            module_logger.debug('Finished iteration: {}/{}, displacement two-norm diff: {}/{}'.format(
                i, max_iterations, two_norm_diff, tolerance))

            update = new_force()
            update *= ((i+1)/200) if i < 200 else 1
            self.node_forces = update
            self.set_forces(update)

            forces_norm = np.linalg.norm(update, axis=1)
            module_logger.debug(
                'min, max, sum forces after pressure applied:\n{}, {}, {}'.format(
                    np.min(forces_norm), np.max(forces_norm),
                    np.sum(forces_norm)))
            module_logger.debug(
                'forces after pressure applied sample:\n{}'.format(
                    update[self._debug_samples,:]))

            if write and i % 5 == 0:
                    self.write_output_vtk('init_{:0>5}'.format(i))

            if two_norm_diff < tolerance and i > 500:
                module_logger.debug('Converged in {}/{} iterations'.format(
                    i, max_iterations))
                if write:
                    self.write_output_vtk('init_final')
                return True

        module_logger.debug('No convergence in {} iterations'.format(max_iterations))

        return False

    #TODO
    def get_cell_areas(self):
        if self.initialized:
            areas = self.mesh.calc_shell_areas(
                vertices=np.reshape(self.nodal.n_x, (-1, 3)))
        else:
            areas = self.mesh.calc_shell_areas()
        return areas

    #TODO
    def get_cell_centers(self):
        if self.initialized:
            cell_centers = self.mesh.calc_shell_centers(
                n_x=np.reshape(self.nodal.n_x, (-1, 3)))
        else:
            cell_centers = self.mesh.calc_shell_centers()
        return cell_centers

    # TODO
    def get_node_normals(self):
        if self.initialized:
            cell_normals = self.mesh.calc_shell_normals(
                n_x=self.get_nodes(),normalize=True)
        else:
            cell_normals = self.mesh.calc_shell_normals(normalize=True)

        normals = np.zeros((self.mesh.number_of_nodes(),3))
        for i, nodes in enumerate(self.mesh.shells):
            normals[nodes] += cell_normals[i]
        normals = normalize_vectors(normals)

        return normals


class PreciceHelper:
    def __init__(self, path):
        self.interface = None
        self.config_path = path
        self.dimensions = 0
        self.num_vertices = 0
        self.vertex_ids = 0
        self.quad_ids = 0
        self.displacement_id = 0
        self.displacement = None
        self.force_id = 0
        self.force = None
        self.time_step = 0

    def setup_interface(self, solver_name='Structure_Solver'):
        print(solver_name, self.config_path)
        self.interface = precice.Interface(
            solver_name, str(self.config_path), 0, 1)

    def configure_interface(self, nodes, grid_name='Structure_Nodes',
                            quads=None):
        self.num_vertices = len(nodes)
        self.dimensions = self.interface.get_dimensions()

        mesh_id = self.interface.get_mesh_id(grid_name)
        vertices = nodes

        self.displacement = np.zeros((self.num_vertices, self.dimensions))
        self.force = np.zeros((self.num_vertices, self.dimensions))

        self.vertex_ids = self.interface.set_mesh_vertices(
            mesh_id, vertices)

        module_logger.debug('precice vertex ids:\n %s' % str(self.vertex_ids))

        if not isinstance(quads, type(None)):
            for ids in quads:
                self.quad_ids = self.interface.set_mesh_quad_with_edges(
                    mesh_id, ids[0], ids[1], ids[2], ids[3])

        self.displacement_id = self.interface.get_data_id(
            'DisplacementDelta', mesh_id)
        self.force_id = self.interface.get_data_id(
            'Force', mesh_id)

        self.time_step = self.interface.initialize()

        if self.interface.is_read_data_available():
            self.interface.read_block_vector_data(self.force_id,
                                                  self.vertex_ids)

    def initialize_data(self):
        if self.interface.is_action_required(
                precice.action_write_initial_data()):
            self.interface.write_block_vector_data(self.displacement_id,
                                                   self.vertex_ids,
                                                   self.displacement)
            self.interface.mark_action_fulfilled(
                precice.action_write_initial_data())

        self.interface.initialize_data()

    def get_participant_name_from_xml(self):
        with open(self.config_path, 'r', encoding='utf8') as file:
            content = file.read()
            soup = BeautifulSoup(content, 'xml')
            participants = soup.find_all("participant")
            for solver in participants:
                if 'Force' in solver.find('read-data').attrs['name']:
                    name = solver.attrs['name']
        return name

    def get_mesh_name_from_xml(self):
        with open(self.config_path, 'r', encoding='utf8') as file:
            content = file.read()
            soup = BeautifulSoup(content, 'xml')
            participants = soup.find_all("participant")
            for solver in participants:
                if 'Force' in solver.find('read-data').attrs['name']:
                    mesh_names = solver.find_all('use-mesh')
                    for mesh in mesh_names:
                        if 'provide' in mesh.attrs:
                            name = mesh.attrs['name']
        return name

    def advance_time(self):
        print("MBDyn Adapter: Advancing in time")
        self.time_step = self.interface.advance(self.time_step)

    def read_data(self):
        if self.interface.is_read_data_available():
            self.force = self.interface.read_block_vector_data(
                self.force_id, self.vertex_ids)

    def write_data(self, write_data):
        if self.interface.is_write_data_required(self.time_step):
            self.interface.write_block_vector_data(
                self.displacement_id, self.vertex_ids, write_data)


class Mesh:
    def __init__(self):
        self.name = ''
        self.nodes = np.array(None)
        self.node_constraints = np.array(None)
        self.node_orientations = np.array(None)
        self.edges = np.array(None)
        self.edge_names = np.array(None)
        self.shells = np.array(None)
        self.shell_names = np.array(None)


    def constraints_from_edge_names(self):
        self.node_constraints = np.full((len(self.nodes), 6), False,
                                        dtype='?')
        for idx in range(len(self.edges)):
            cur_constraints = np.full(6, False, dtype='?')
            cur_name = self.edge_names[idx].casefold()
            if 'fix' in cur_name:
                cur_name = cur_name.replace('fix', '')
                if cur_name == 'all':
                    cur_constraints[:] = True
                else:
                    if 'x' in cur_name:
                        cur_constraints[0] = True
                    if 'y' in cur_name:
                        cur_constraints[1] = True
                    if 'z' in cur_name:
                        cur_constraints[2] = True
                    if 'a' in cur_name:
                        cur_constraints[3] = True
                    if 'b' in cur_name:
                        cur_constraints[4] = True
                    if 'c' in cur_name:
                        cur_constraints[5] = True
            for node in self.edges[idx]:
                self.node_constraints[node, :] += cur_constraints

    def number_of_nodes(self):
        return len(self.nodes)

    def number_of_shells(self):
        return len(self.shells)

    def calc_shell_normals(self, n_x=None, normalize=False, invert=1):
        if isinstance(n_x, type(None)):
            n_x = self.nodes

        triangles = np.append(self.shells[:, :3],
                              self.shells[:, [0, 2, 3]], axis=0)
        triangle_proj = 0.5 * np.cross(
            (n_x[triangles[:, 1]] - n_x[triangles[:, 0]]),
            (n_x[triangles[:, 2]] - n_x[triangles[:, 0]]))

        num_triangles = self.number_of_shells()
        proj = triangle_proj[:num_triangles] + triangle_proj[num_triangles:]
        proj = proj * invert

        if normalize:
            proj = normalize_vectors(proj)

        return proj

    def calc_shell_areas(self, vertices=None):
        if vertices is None:
            normals = self.calc_shell_normals()
        else:
            normals = self.calc_shell_normals(n_x=vertices)

        return np.linalg.norm(normals, axis=1)

    def calc_shell_centers(self, n_x=None, shape='quadrilaterals'):
        if n_x is None:
            n_x = self.nodes

        if shape == 'triangels':
            triangles_from_shells = np.append(
                self.shells[:, :3], self.shells[:, [0, 2, 3]], axis=0)
            cell_centers = np.sum(n_x[triangles_from_shells], axis=1) / 3.0

        if shape == 'quadrilaterals':
            cell_centers = np.sum(n_x[self.shells], axis=1) / 4.0

        return cell_centers

    def calc_node_orientation(self, unchanged='z', clean_unchanged=True,
                              return_normal=False, flip_normal=1):
        valid_unchanged = ['x', 'y', 'z']
        if unchanged not in valid_unchanged:
            raise ValueError('unchanged must be one of {}.'.format(
                unchanged))

        cell_normals = self.calc_shell_normals() * flip_normal

        if clean_unchanged:
            cell_normals[:, valid_unchanged.index(unchanged)] = 0

        cell_normals = normalize_vectors(cell_normals)

        node_normals = np.zeros((self.number_of_nodes(), 3))

        for i, vertices in enumerate(self.shells):
            node_normals[vertices, :] += cell_normals[i, :]

        node_normals = normalize_vectors(node_normals)

        # TODO: which coordinate points outwards
        global_frame = np.zeros((2, 3))
        global_frame[0, 1] = 1
        global_frame[1, 2] = 1

        local_frame = global_frame.copy()

        orientation = np.zeros((self.number_of_nodes(), 3))

        for i, normal in enumerate(node_normals):
            local_frame[0, :] = normal
            rotation = R.align_vectors(global_frame, local_frame)
            orientation[i, :] = rotation[0].as_euler('xyz')

        self.node_orientations = orientation

        if return_normal:
            return orientation, node_normals
        return orientation

    def match_names(self, names, edgenumbers, shellnumbers):
        self.edge_names = np.empty(len(edgenumbers), dtype='U')
        self.shell_names = np.empty(len(shellnumbers), dtype='U')

        print('Named regions found:')
        print('{0:<15} {1:<10} {2:<10}'.format('Name', 'Edges', 'Shells'))

        for name in names:
            self.edge_names = np.where(edgenumbers == name[1],
                                       name[2], self.edge_names)
            self.shell_names = np.where(shellnumbers == name[1],
                                        name[2], self.shell_names)
            n_edgetype = np.count_nonzero(edgenumbers == name[1])
            n_shelltype = np.count_nonzero(shellnumbers == name[1])
            print('{:<15} {:<10} {:<10}'.format(name[2], n_edgetype,
                                                n_shelltype))

        print('{:-^36}'.format(''))
        print('{:<15} {:<10} {:<10}'.format('Sum', len(edgenumbers),
                                            len(shellnumbers)))

    #TODO: fix if, test it
    def set_clamp_constraint(self, fixed_nodes, dead_z=False):
        assert(isinstance(fixed_nodes, (slice, list, int)))
        if not self.node_constraints.any():
            self.node_constraints = np.full(
                (self.number_of_nodes(), 6),  False, dtype='?')
        self.node_constraints[fixed_nodes, :3] = True
        if dead_z:
            self.node_constraints[:, 2] = True

def normalize_vectors(vectors):
    length = np.linalg.norm(vectors, axis=1)
    return np.divide(vectors.transpose(), length).transpose()


if __name__ == "__main__":
    pass
