#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation as R

from mbc_py_interface import mbcNodal


def visualize_mesh_3d(nodes, normals, scaling=1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)

    # plot nodes
    ax.plot(nodes[:, 0], nodes[:, 1], nodes[:, 2], 'o',
            markersize=10, color='green', alpha=0.2)

    # plot normals
    for i, vector in enumerate(normals):
        ax.quiver(
            nodes[i, 0], nodes[i, 1], nodes[i, 2],
            vector[0] * scaling, vector[1] * scaling, vector[2] * scaling,
            color='red', alpha=.8, lw=3,
        )

    # labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # settings
    ax.view_init(90, 0)

    # title
    plt.title('Mesh normals')

    # diplay
    plt.draw()
    plt.show()


def visualize_mesh_3d_reduced(nodes, normals, scaling=1):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot nodes
    ax.plot(nodes[:, 0], nodes[:, 1], 'o',
            markersize=1, color='green', alpha=0.2)

    # plot normals
    for i, vector in enumerate(normals):
        ax.quiver(
            nodes[i, 0], nodes[i, 1],
            vector[0] * scaling, vector[1] * scaling,
            color='red', alpha=.8, width=0.001,
        )

    # labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    # settings
    ax.set_aspect('equal')

    padding = 1.2
    bottom, top = ax.get_ylim()
    left, right = ax.get_xlim()
    middle = [(right + left) / 2, (top + bottom) / 2]
    height = (top - bottom) * padding
    width = (right - left) * padding / 2

    ax.set_xlim(middle[0] - width, middle[0] + width)
    ax.set_ylim(middle[1] - height, middle[1] + height)

    # title
    plt.title('Mesh normals')

    # diplay
    plt.draw()
    plt.show()


def animate_nodes_and_normals(nodes, normals, num_dt, num_nodes):

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(0.2, 0.65), ylim=(0, 0.41))
    dots, = ax.plot([], [], lw=2, ls='', marker='.')
    arrows = ax.quiver(nodes[0, :, 0], nodes[0, :, 1],
                       normals[0, :, 0], normals[0, :, 1],
                       color='r', units='xy', width=0.0025)

    double_nodes = np.zeros((num_dt, 2*num_nodes, 3))
    double_normals = np.zeros((num_dt, 2*num_nodes, 3))

    double_nodes[:, :num_nodes, :] = nodes
    double_nodes[:, num_nodes:, :] = nodes
    double_normals[:, :num_nodes, :] = normals * 0.01
    double_normals[:, num_nodes:, :] = -normals * 0.01

    # initialization function: plot the background of each frame
    def init():
        dots.set_data([], [])
        arrows.set_offsets([])
        arrows.set_UVC([], [])
        return dots, arrows

    # animation function.  This is called sequentially
    def animate(i):
        dots.set_data(nodes[i, :, 0], nodes[i, :, 1])
        arrows.set_offsets(double_nodes[i, :, 0:2])
        arrows.set_UVC(double_normals[i, :, 0], double_normals[i, :, 1])
        return dots, arrows

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=num_dt, interval=30, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save('basic_animation.mp4', fps=30,
              extra_args=['-vcodec', 'libx264'], dpi=200)

    plt.show()


def extract_nodes_and_normals(file_name):
    data = np.genfromtxt(file_name, delimiter=' ')

    num_nodes = int(np.max(data[:, 0])) + 1
    num_time_steps = int(len(data[:, 0]) / num_nodes)

    nodes_pos = np.zeros((num_time_steps, num_nodes, 3))
    for i in range(num_time_steps):
        nodes_pos[i, :, :] = data[num_nodes*i:num_nodes*(i+1), 1:4]

    nodes_euler = np.zeros((num_time_steps, num_nodes, 3))
    for i in range(num_time_steps):
        nodes_euler[i, :, :] = data[num_nodes*i:num_nodes*(i+1), 4:7]

    nodes_euler[:, :, 0] = nodes_euler[:, :, 0] * -1

    nodes_normals = np.zeros((num_time_steps, num_nodes, 3))
    normal = np.array([0, 1, 0])
    for i in range(num_time_steps):
        tmp_rotation = R.from_euler('xyz', nodes_euler[i, :, :], degrees=1)
        nodes_normals[i, :, :] = tmp_rotation.apply(normal)

    animate_nodes_and_normals(nodes_pos, nodes_normals,
                              num_time_steps, num_nodes)

    return nodes_normals, tmp_rotation.as_matrix()
