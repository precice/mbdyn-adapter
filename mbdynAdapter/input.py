#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class MBDynInput:
    '''This is a quick summary line used as a description of the object.'''

    def __init__(self):
        self._indexing = 0
        self._file = 'case'
        self.data = DataBlock()
        self.problem = ProblemBlock()
        self.control = ControlBlock()
        self.nodes = NodesBlock()
        # drivers block not implemented
        self.elements = ElementsBlock()
        # parallel block not implemented

    def create_from_prep(self, prep):
        self._file = prep.name

        self.set_data_block()
        self.set_problem_block(prep.problem_dict)
        self.set_nodes_block(prep.nodes_dict, prep.material_dict, prep.mesh)
        self.set_elements_block(prep.material_dict, prep.mesh)
        self.set_control_block(prep.control_dict, prep.material_dict)

    def set_data_block(self, new_problem='initial value'):
        self.data.entries['problem'] = new_problem

    def set_problem_block(self, problem_dict):
        self.problem.name = self.data.entries['problem']
        for key, value in problem_dict.items():
            self.problem.entries[key] = value

    def set_control_block(self, control_dict, material_dict):
        self.count_cards_by_type()
        for key, value in control_dict.items():
            self.control.entries[key] = value
        for key, value in material_dict.items():
            self.control.variables[key] = float(value)

    def set_nodes_block(self, nodes_dict, material_dict, mesh):
        self._indexing = 10**(int(np.ceil(np.log10(mesh.number_of_nodes()))))

        self.nodes.nodes_from_mesh(
            mesh, set_orientation=int(nodes_dict['orientation']))

        # check if damping is set
        if 'C' in material_dict.keys():
            if 'damping node' in nodes_dict.keys():
                self.nodes.create_node(
                    mesh.number_of_nodes(), 'static',
                    self.nodes.vector_to_string(
                        np.fromstring(nodes_dict['damping node'], sep=',')))
            else:
                self.nodes.create_node(
                    mesh.number_of_nodes(), 'static',
                    np.fromstring('0, 0, 0'))

    def set_elements_block(self, material_dict, mesh):
        assert self._indexing > 0, 'NodesBlock needs to be set first.'

        # check if damping is set
        if 'C' in material_dict.keys():
            self.elements.damping_from_mesh(mesh, self._indexing,
                                            mesh.number_of_nodes())

        self.elements.bodies_from_mesh(mesh, self._indexing)
        self.elements.constraints_from_mesh(mesh, self._indexing)
        self.elements.plates_from_mesh(mesh, self._indexing)
        self.elements.create_force_coupling(mesh.number_of_nodes(),
                                            self._file, self._indexing)

    def count_cards_by_type(self):
        count = dict()
        for key, value in self.nodes.entries.items():
            count[key] = len(value)

        for key, value in self.elements.entries.items():
            count[key] = len(value)

        count_keys = count.keys()
        if 'membrane4eas' in count_keys:
            self.control.entries['plates'] += count.pop('membrane4eas')

        if 'shell4easans' in count_keys:
            self.control.entries['plates'] += count.pop('shell4easans')

        if 'force' in count_keys:
            self.control.entries['forces'] += count.pop('force')

        if 'joint' in count_keys:
            self.control.entries['joints'] += count.pop('joint')

        if 'body' in count_keys:
            self.control.entries['rigid bodies'] += count.pop('body')

        if 'structural' in count_keys:
            self.control.entries['structural nodes'] += count.pop('structural')

        if 'gravity' in count_keys:
            self.control.entries['gravity'] = None

        if len(count) > 0:
            raise ValueError('unknown type found during count')

    def update_time_step(self, time_step):
        self.problem.entries['time step'] = time_step

    def write_input_file(self, file_name=None):
        if file_name:
            self._file = file_name
        with open(self._file + '.mbd', 'w') as input_file:
            input_file.write(self.data.get_block_str())
            input_file.write(self.problem.get_block_str())
            input_file.write(self.control.get_block_str())
            input_file.write(self.nodes.get_block_str())
            input_file.write(self.elements.get_block_str())


class Block:
    indent_str = '    '
    new_line_str = '\n'
    line_end_str = ';{}'.format(new_line_str)

    def __init__(self):
        self.name = ''
        self.entries = dict()

    def __len__(self):
        return len(self.entries)

    def _open_str(self):
        return 'begin: {name}{lend}'.format(name=self.name,
                                            lend=self.line_end_str)

    def _close_str(self):
        return 'end: {name}{lend}'.format(name=self.name,
                                          lend=self.line_end_str)

    def _block_sep_str(self):
        return '{nline}#{line}{nline}# [{bname} BLOCK]{nline}{nline}'.format(
            line='-'*77, nline=self.new_line_str, bname=self.name.upper())

    def get_block_str(self):
        block_str = self._block_sep_str() + self._open_str()
        for key, value in self.entries.items():
            if value is None or value == '':
                block_str += '{idt}{card}{lend}'.format(
                    idt=self.indent_str, card=key, lend=self.line_end_str)
            elif isinstance(value, list):
                for item in value:
                    block_str += '{idt}{card}: {val}{lend}'.format(
                        idt=self.indent_str, card=key, val=item,
                        lend=self.line_end_str)
            else:
                block_str += '{idt}{card}: {val}{lend}'.format(
                    idt=self.indent_str, card=key, val=value,
                    lend=self.line_end_str)
        block_str += self._close_str()
        return block_str

    def vector_to_string(self, vector):
        out = np.array2string(vector, separator=',', threshold=len(vector),
                              formatter={'float_kind': '{: }'.format},
                              max_line_width=80)
        return self.shave(out)

    def matrix_to_string(self, matrix):
        array_str = np.array2string(matrix.ravel(), separator=',',
                                    formatter={'float_kind': '{: }'.format})
        return 'matr,{}'.format(self.shave(array_str))

    def euler_to_string(self, euler, order='123'):
        euler_str = np.array2string(euler.ravel(), separator=',',
                                    formatter={'float_kind': '{: }'.format})
        return 'euler{},{}'.format(order, self.shave(euler_str))

    def shave(self, input_string):
        remove_str = '''["']''' + self.new_line_str
        table = str.maketrans('', '', remove_str)
        return input_string.translate(table)


class DataBlock(Block):
    def __init__(self):
        super().__init__()
        self.name = 'data'
        self.entries['problem'] = 'initial value'


class ProblemBlock(Block):
    def __init__(self):
        super().__init__()
        self.name = 'initial value'
        self.entries['initial time'] = 0
        self.entries['final time'] = 'forever'
        self.entries['time step'] = 0
        self.entries['method'] = 'ms, 0.6'
        self.entries['tolerance'] = 1e-4
        self.entries['max iterations'] = 300
        self.entries['linear solver'] = 'umfpack'
        self.entries['output'] = 'iterations'
        self.entries['threads'] = 'auto'
        self.entries['modify residual test'] = None


class ControlBlock(Block):
    def __init__(self):
        super().__init__()
        self.name = 'control data'
        self.entries['beams'] = 0
        self.entries['forces'] = 0
        self.entries['joints'] = 0
        self.entries['plates'] = 0
        self.entries['rigid bodies'] = 0
        self.entries['structural nodes'] = 0
        self.entries['default output'] = 'none, plates'
        self.entries['output frequency'] = 10
        self.variables = dict()

    def get_block_str(self):
        block_str = super().get_block_str()
        if len(self.variables) > 0:
            block_str += self.new_line_str
            for var_label, var_value in self.variables.items():
                var_type = 'real' if isinstance(var_value, float) else 'int'
                block_str += 'set: {nset} {lbl}{spc} = {val}{lend}'.format(
                    nset=var_type, lbl=var_label,
                    spc=(11-len(var_label)) * ' ', val=var_value,
                    lend=self.line_end_str)
        return block_str


class NodesBlock(Block):
    __dynamic_node_str = '{nindex}, {ntype},{npos}, {ndir}, {v0}, {o0}'
    __static_node_str = '{nindex}, {ntype},{npos}, {ndir}, {v0}, {o0}, output, no'
    __displacement_node_str = '{nindex}, {ntype},{npos}, {v0}'

    def __init__(self):
        super().__init__()
        self.name = 'nodes'
        self.entries['structural'] = list()

    def __len__(self):
        return len(self.entries['structural'])

    def create_node(self, index, node_type, coordinates,
                    orientation='eye', v_0='null', omega_0='null'):
        if node_type == 'dynamic':
            node_str = self.__dynamic_node_str.format(
                nindex=index, ntype=node_type, npos=coordinates,
                ndir=orientation, v0=v_0, o0=omega_0)
        elif node_type == 'static':
            node_str = self.__static_node_str.format(
                nindex=index, ntype=node_type, npos=coordinates,
                ndir=orientation, v0=v_0, o0=omega_0)
        elif node_type == 'dynamic displacement':
            node_str = self.__displacement_node_str.format(
                nindex=index, ntype=node_type, npos=coordinates, v0=v_0)
        self.entries['structural'].append(node_str)

    def remove_node(self, node_index):
        for i, node_str in self.entries['structural'].items():
            if node_str.split(',')[0] == node_index:
                self.entries['structural'].pop(i)

    # TODO: constraint free membrane node as dynamic displacement
    def nodes_from_mesh(self, input_mesh, node_type='dynamic',
                        set_orientation=0):
        if set_orientation:
            input_mesh.calc_node_orientation(flip_normal=set_orientation)
            for i, xyz in enumerate(input_mesh.nodes):
                self.create_node(i, node_type,
                                 self.vector_to_string(xyz),
                                 self.euler_to_string(
                                     input_mesh.node_orientations[i, :]))
        else:
            for i, xyz in enumerate(input_mesh.nodes):
                self.create_node(i, node_type, self.vector_to_string(xyz))


class ElementsBlock(Block):
    __body_str = '{blabel}, {nlabel}, {bmass}, {coffset}, {binertia}'
    __gravity_str = '{vec}, const, {val}'
    # orientation matrix
    __force_coupling_str = '{flabel}, external structural, socket, create, yes, path, "{sname}.sock", coupling, tight, labels, no, orientation, none, accelerations, no, {num},{nodes}'
    __plate_str = '{mlabel},{mcorners}, {constlaw}{mstress}'
    __total_pin_joint_str = '''{jlabel}, total pin joint, {jnode}, position, reference, node, null, position orientation, reference, node, eye, rotation orientation, reference, node, eye, position, {jpos}, position orientation, {jorientation}, rotation orientation, {jorientation}, position constraint,{pstraint}, null, orientation constraint,{ostraint}, null'''
    __deformable_displacement_joint_str = '{jlabel}, deformable displacement joint, {node1}, null, {node2}, null, {constlaw}'
    __clamp_str = '{clabel}, clamp, {cnode}, node, node'

    def __init__(self):
        super().__init__()
        self.name = 'elements'
        self.entries['body'] = list()
        self.entries['force'] = list()
        self.entries['membrane4eas'] = list()
        self.entries['shell4easans'] = list()
        self.entries['joint'] = list()

    def create_body(self, label, node, mass, offset, inertia):
        body_str = self.__body_str.format(blabel=label, nlabel=node,
                                          bmass=mass, coffset=offset,
                                          binertia=inertia)
        self.entries['body'].append(body_str)

    # TODO: add inertia to body option, custom variables
    def bodies_from_mesh(self, input_mesh, indexing):
        node_areas = np.zeros(input_mesh.number_of_nodes())
        for i, shell_area in enumerate(input_mesh.calc_shell_areas()):
            node_areas[input_mesh.shells[i, :]] += shell_area / 4

        index_offset = 6 * indexing
        for i, node_area in enumerate(node_areas):
            self.create_body(index_offset+i, i, '{}*rho*t'.format(node_area),
                             'null', 'diag, 0.0, 0.0, 0.0')

    def create_gravity(self, direction=np.array([0, -1, 0]), value=9.81):
        gravity_str = self.__gravity_str.format(
            vec=self.vector_to_string(direction), val=value)
        self.entries['gravity'] = gravity_str

    def create_force_coupling(self, num_nodes, case_name, indexing):
        index_offset = 5 * indexing
        nodes_str = self.vector_to_string(np.arange(num_nodes))
        force_str = self.__force_coupling_str.format(
            flabel=index_offset, sname=case_name, num=num_nodes,
            nodes=self.new_line_str+nodes_str)
        self.entries['force'].append(force_str)

    def create_plate(self, label, vertices, constitutive_law,
                     prestress='', key='shell4easans'):
        if prestress:
            prestress = ', prestress, ' + prestress

        plate_str = self.__plate_str.format(
            mlabel=label, mcorners=vertices, constlaw=constitutive_law,
            mstress=prestress)
        self.entries[key].append(plate_str)

    # TODO: custom variables, add prestress
    def plates_from_mesh(self, input_mesh, indexing):
        index_offset = 4 * indexing
        isotropic_law = 'isotropic, E, Et, nu, nut, thickness, t'
        for i, vertices in enumerate(input_mesh.shells):
            # TODO quick and dirty fix sofar
            try:
                plate_type = input_mesh.shell_names[i].casefold()
            except IndexError:
                plate_type = 'membrane'
            if 'shell' in plate_type:
                self.create_plate(index_offset+i,
                                  self.vector_to_string(vertices),
                                  isotropic_law)
            elif 'membrane' in plate_type:
                self.create_plate(index_offset+i,
                                  self.vector_to_string(vertices),
                                  isotropic_law, key='membrane4eas')

    def create_total_pin_joint(self, label, node, joint_position,
                               position_constraint, orientation_constraint,
                               node_orientation='eye'):
        joint_str = self.__total_pin_joint_str.format(
            jlabel=label, jnode=node, jpos=joint_position,
            pstraint=position_constraint, ostraint=orientation_constraint,
            jorientation=node_orientation)
        self.entries['joint'].append(joint_str)

    # TODO: constraints for dynamic displacement nodes if possible
    def constraints_from_mesh(self, input_mesh, indexing):
        index_offset = 3 * indexing

        constraint_str = input_mesh.node_constraints.astype(np.str_)
        constraint_str = np.char.replace(constraint_str, 'True', 'active')
        constraint_str = np.char.replace(constraint_str, 'False', 'inactive')

        for i, xyz in enumerate(input_mesh.nodes):
            self.create_total_pin_joint(
                index_offset+i, i, self.vector_to_string(xyz),
                self.vector_to_string(constraint_str[i, :3]),
                self.vector_to_string(constraint_str[i, 3:]))

    def create_deformable_displacement_joint(self, label, first_node,
                                             second_node, constitutive_law):
        joint_str = self.__deformable_displacement_joint_str.format(
            jlabel=label, node1=first_node, node2=second_node,
            constlaw=constitutive_law)
        self.entries['joint'].append(joint_str)

    def create_clamp(self, label, node):
        clamp_str = self.__clamp_str.format(clabel=label, cnode=node)
        self.entries['joint'].append(clamp_str)

    # TODO: set which constitutive law option
    def damping_from_mesh(self, input_mesh, indexing, fixed_point_label):
        index_offset = 2 * indexing + 1
        self.create_clamp(index_offset-1, fixed_point_label)
        cl_linviseliso_str = 'linear viscoelastic isotropic, K, C'
        cl_linvis_str = 'linear viscous, C'
        for i in range(input_mesh.number_of_nodes()):
            self.create_deformable_displacement_joint(
                label=index_offset+i, first_node=fixed_point_label,
                second_node=i, constitutive_law=cl_linvis_str)


if __name__ == "__main__":
    pass
