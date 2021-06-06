#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import configparser
import numpy as np
import meshio

from mbdynprecice.helper import Mesh


class MBDynPrep:
    ''' A simple preprocessor that provides the data
    for the precice-mbdyn adapter. '''

    def __init__(self, case_name=None, in_mm=False):
        self.name = ''
        self.mesh = Mesh()
        self.problem_dict = dict()
        self.control_dict = dict()
        self.nodes_dict = dict()
        self.material_dict = dict()

        # TODO: case name input instead of mesh name
        if case_name:
            if os.path.dirname(case_name):
                directory = os.path.dirname(case_name)
                config_name = os.path.join(directory, 'mbdyn-config')
            else:
                config_name = 'mbdyn-config'
            self.read_config(config_name)
            self.read_mesh(case_name, self.nodes_dict.get('fixed nodes'),
                           mm_to_m=in_mm)

    def read_config(self, file_name='mbdyn-config'):
        parser = configparser.ConfigParser()
        parser.optionxform = str
        parser.read(file_name)
        self.control_dict = dict(dict(parser)['control'])
        self.problem_dict = dict(dict(parser)['problem'])
        self.nodes_dict = dict(dict(parser)['nodes'])
        self.material_dict = dict(dict(parser)['material'])
        if 'fixed nodes' in self.nodes_dict.keys():
            if ':' in self.nodes_dict['fixed nodes']:
                splitter = self.nodes_dict['fixed nodes'].split(':')
                self.nodes_dict['fixed nodes'] = slice(int(splitter[0]),
                                                       int(splitter[1]), 1)
            elif ',' in self.nodes_dict['fixed nodes']:
                splitter = self.nodes_dict['fixed nodes'].split(',')
                self.nodes_dict['fixed nodes'] = list()
                for node in splitter:
                    self.nodes_dict['fixed nodes'].append(int(node) - 1)
            else:
                self.nodes_dict['fixed nodes'] = int(
                    self.nodes_dict['fixed nodes'])

    def read_mesh(self, file_name: str, fixed_nodes=None, mm_to_m=False):
        mesh_type = meshio.extension_to_filetype.get(
            os.path.splitext(file_name)[1])
        assert not isinstance(mesh_type, type(None)), \
            "Mesh file format of '{}' uncompatible!".format(file_name)

        print("Reading mesh file '{}' of type '{}'".format(
            os.path.basename(file_name), mesh_type))

        mfile = meshio.read(file_name)
        self.name = os.path.splitext(os.path.basename(file_name))[0]
        self.mesh.name = file_name

        assert len(mfile.points) > 0, \
            "Mesh does not contain any nodes!"

        self.mesh.nodes = mfile.points
        print("Number of nodes: {}".format(len(self.mesh.nodes)))
        if mm_to_m:
            self.mesh.nodes *= 0.001

        assert 'quad' in mfile.cells_dict, \
            "Mesh does not contain any quadrilateral elements!"

        self.mesh.shells = mfile.get_cells_type('quad')
        print("Number of quads: {}".format(len(self.mesh.shells)))

        if 'line' in mfile.cells_dict:
            self.mesh.edges = mfile.get_cells_type('line')
            print("Number of lines: {}".format(len(self.mesh.edges)))

        if mesh_type == 'gmsh' and mfile.field_data and \
                all(tp in mfile.cells_dict for tp in ('line', 'quad')):
            print("Reading physical names if 'gmsh' file:")

            edge_type = mfile.get_cell_data('gmsh:physical', 'line')
            shell_type = mfile.get_cell_data('gmsh:physical', 'quad')
            physical_names = mfile.field_data

            tab_str = '{0:<15} {1:<10} {2:<10}'

            print("Named elements found:")
            print(tab_str.format('Name', 'Lines', 'Quads'))

            self.mesh.edge_names = np.empty(edge_type.shape, dtype='<U8')
            self.mesh.shell_names = np.empty(shell_type.shape, dtype='<U8')

            for name, val in physical_names.items():
                element_number = val[0]
                element_type = val[1]
                if element_type == 1:
                    matching = (edge_type == element_number)
                    self.mesh.edge_names[matching] = name
                    print(tab_str.format(name, np.count_nonzero(matching), 0))
                elif element_type == 2:
                    matching = (shell_type == element_number)
                    self.mesh.shell_names[matching] = name
                    print(tab_str.format(name, 0, np.count_nonzero(matching)))

            print("{:-^36}".format(''))
            print(tab_str.format('Sum', np.count_nonzero(edge_type),
                                 np.count_nonzero(shell_type)))

            self.mesh.constraints_from_edge_names()

        if fixed_nodes is not None:
            self.mesh.set_clamp_constraint(fixed_nodes)


if __name__ == '__main__':
    prep = MBDynPrep('test-files/kite-v5-1500cells.msh')
