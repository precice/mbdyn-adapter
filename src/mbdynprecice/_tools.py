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

# %%


def n_triangles(points):
    v1, v2, v3 = points[1:, :] - points[0, :]
    normal = 0.5 * (np.cross(v1, v2) + np.cross(v2, v3))

    area = np.linalg.norm(normal)

    return normal / area


def n_lsquare_height(points):
    # Note: only for 3D planes
    average = np.average(points, axis=0)
    reduced = points-average

    matrix_a = reduced[:, :2]
    vector_b = reduced[:, 2]

    inverse = np.linalg.inv(np.dot(matrix_a.T, matrix_a))
    a_tilde, b_tilde = inverse.dot(np.dot(matrix_a.T, vector_b))

    normal = np.array([-a_tilde, -b_tilde, 1])

    return normal / np.linalg.norm(normal)


def n_lsquare_orthogonal(points):
    # Note: points.shape = (m_samples, n_coordinates)
    m, n = points.shape

    average = np.average(points, axis=0)
    reduced = np.asmatrix(points-average)

    c_matrix = np.zeros((n, n))
    for i in range(m):
        c_matrix += np.matmul(reduced[i, :].T, reduced[i, :])
    
    test = np.cov(reduced.T)

    eigen_val, eigen_vec = np.linalg.eig(c_matrix)
    smallest = np.argmin(eigen_val)

    normal = eigen_vec[:, smallest]

    return normal / np.linalg.norm(normal)

# %%


def plane_with_noise_data(n, delta_x=0.25, delta_y=0.5, offset=5, length=5):
    noise = 0.5

    data = np.zeros((n, 3))
    data[:, 0] = [np.random.uniform(2*length)-length for i in range(n)]
    data[:, 1] = [np.random.uniform(2*length)-length for i in range(n)]
    for i in range(n):
        data[i, 2] = data[i, 0] * delta_x + data[i, 1] * delta_y \
            + offset + np.random.normal(scale=noise)
    return data

import meshio
from functools import lru_cache
from timeit import default_timer

def fsi_data(f_path):
    mesh = meshio.read(f_path)
    if 'quad' in mesh.cells_dict:
        for quad in mesh.cells_dict['quad']:
            yield mesh.points[quad]

def normal_algorithms(org_vtk, def_vtk):
    vtk_locs = [org_vtk, def_vtk]
    func = {'triangles':    n_triangles,
            'height':       n_lsquare_height,
            'orthogonal':   n_lsquare_orthogonal}
    normal = {'triangles':  np.empty((0,3)),
              'height':     np.empty((0,3)),
              'orthogonal': np.empty((0,3))}
    for method in func:
        print("Using method {}!".format(method))
        for vtk in vtk_locs:
            quads = fsi_data(vtk)
            start = default_timer()
            for quad in quads:
                try:
                    n = func[method](quad)
                except:
                    print('prop singulary')
                    n = np.array([0, 0, 0])
                normal[method] = np.append(
                    normal[method], n.reshape((1,-1)), axis=0)
            end = default_timer()
            print("Elapsed time: {:.3f}".format((end-start)*1000))
    return normal

def compare_algorithms():
    vtk_0 = 'test-files/kite-v5-1500cells_00001.vtk'
    vtk_1 = 'test-files/kite-v5-1500cells_01209.vtk'
    
    normals = normal_algorithms(vtk_0, vtk_1)
    alphas = np.empty((0))
    
    tri, zdir, ortho = normals.keys()
    for ls in [zdir, ortho]:
        for i, normal in enumerate(normals[ls]):
            scalar = np.dot(normals[tri][i,:], normal)
            if abs(scalar) > 1:
                scalar /= abs(scalar)
            angle = np.arccos(scalar) * 180 / np.pi
            alphas = np.append(alphas, angle)
    
    indices = alphas > 90
    alphas -= 180 * indices
    alphas = np.absolute(alphas)
    alphas = np.split(alphas, 4)
    indices = np.split(indices, 4)
    for a, i in zip(alphas, indices):
        print("average change: {:.3f}".format(np.average(a)))
        print("max change: {:.3f}".format(np.max(a)))
        print("wrong orientation: {}".format(np.sum(i)))
        
    # fig, ax = plt.subplots()
    # for a in alphas:
    #     ax.plot(a)
    # plt.show()
    
    return alphas
    
erg = compare_algorithms()

#%%
def test_orthogonal_regression(data=None):

    if isinstance(data, type(None)):
        data = plane_with_noise_data(100)

    # data = np.array([[1.12000000e+00, 3.76613819e-03, 1.45161682e-09],
    #     [1.08130613e+00, 9.10306656e-03, 1.31469972e-09],
    #     [1.08130613e+00, 9.10306656e-03, 5.00000013e-02],
    #     [1.12000000e+00, 3.76613819e-03, 5.00000013e-02]])

    # data[2:,1] += 0.2

    # data = np.array([[ 5.2765869e-02,  1.0268845e-01, -4.9462944e-02],
    #     [ 7.1312808e-02,  6.3278370e-02, -4.7897752e-02],
    #     [ 8.3707243e-02,  7.2916396e-02, -2.6163558e-05],
    #     [ 5.2508812e-02,  1.0379118e-01, -4.4938151e-05]])

    data = np.array([[4.3525200e-02,  5.8813106e-02, -5.0000295e-02],
                     [1.1218552e-02,  3.0126929e-02, -4.9997341e-02],
                     [1.1219671e-02,  3.0127788e-02,  2.2274909e-09],
                     [4.3526553e-02,  5.8813721e-02,  2.8943647e-09]])

    # data = np.array([[ 1,  0, 0.2],
    #     [ 0.5,  0, 1.5],
    #     [-0.5,  0, -0.1],
    #     [0.1,  0,  -1]])

    ave = np.average(data, axis=0)

    data_dash = np.asmatrix(data - ave)

    cov_data = np.zeros((3, 3))
    for sample in range(data.shape[0]):
        cov_data += np.matmul(data_dash[sample, :].T,  data_dash[sample, :])

    eigen_val, eigen_vec = np.linalg.eig(cov_data)
    # print(eigen_val, eigen_vec)

    normal_index = np.argmin(eigen_val)
    # normal_index = 1

    normal = eigen_vec[:, normal_index]
    # print(normal)

    import timeit

    t = timeit.Timer(lambda: n_triangles(data))
    print(t.timeit(10000))
    t = timeit.Timer(lambda: n_lsquare_height(data))
    print(t.timeit(10000))
    t = timeit.Timer(lambda: n_lsquare_orthogonal(data))
    print(t.timeit(10000))

    normal_tri = n_triangles(data)
    normal_hei = n_lsquare_height(data)
    normal_or = n_lsquare_orthogonal(data)

    print(normal_tri)
    print(normal_hei, np.arccos(np.dot(normal_tri, normal_hei))*180/np.pi)
    print(normal_or, np.arccos(np.dot(normal_tri, normal_or))*180/np.pi)

    normal = normal_or

    # normal = eigen_vec[:3, normal_index] / np.linalg.norm(eigen_vec[:3, normal_index])

    # normal_tmp = normal
    # normal_tmp[0] = normal[1]
    # normal_tmp[1] = -normal[0]

    # normal = normal_tmp
    # d = -ave.dot(normal)

    # # create x,y
    # xx, yy = np.meshgrid(np.arange(0, 1, 0.1), np.arange(0, 1, 0.1))

    # # calculate corresponding z
    # z = -(-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]

    # plt3d.plot_surface(xx, yy, z)

    # plt3d.quiver(ave[0], ave[1], ave[2],
    #              normal[0], normal[1], normal[2])

    # fig, ax = plt.subplots()

    # ax.quiver(ave[0], ave[1], normal[0]  , normal[1], units='xy')
    # ax.plot(ave[0], ave[1], 'o')
    # ax.plot(ave[0] + normal[0], ave[1] + normal[1], 'o')
    # ax.plot(data[:,0], data[:,1], 'o')
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_aspect('equal')

    fig = plt.figure()
    ax2 = plt.subplot(111, projection='3d')

    # print(np.dot(eigen_vec.T[0], eigen_vec.T[1]))
    # print(np.dot(eigen_vec.T[1], eigen_vec.T[2]))
    # print(np.dot(eigen_vec.T[2], eigen_vec.T[0]))

    for vec in eigen_vec.T:
        ax2.quiver(ave[0], ave[1], ave[2],
                   vec[0]*10, vec[1]*10, vec[2]*10)

    x_min = np.min(data[:, 0])
    x_max = np.max(data[:, 0])

    y_min = np.min(data[:, 1])
    y_max = np.max(data[:, 1])

    z_min = np.min(data[:, 2])
    z_max = np.max(data[:, 2])

    center = [x_min+0.5*(x_max-x_min), y_min+0.5 *
              (y_max-y_min), z_min+0.5*(z_max-z_min)]

    length = 1.2 * np.max([x_max-x_min, y_max-y_min, z_max-z_min])

    ax2.set_xlim(center[0]-length/2, center[0]+length/2)
    ax2.set_ylim(center[1]-length/2, center[1]+length/2)
    ax2.set_zlim(center[2]-length/2, center[2]+length/2)

    # ax2.quiver(ave[0], ave[1], ave[2],
    #            normal[0]*10, normal[1]*10, normal[2]*10)

    ax2.scatter(ave[0], ave[1], ave[2], color='red')
    ax2.scatter(data[:, 0], data[:, 1], data[:, 2], color='blue')

    d = -ave.dot(normal)
    xx, yy = np.meshgrid(np.arange(x_min, x_max, length/10),
                         np.arange(y_min, y_max, length/10))
    # calculate corresponding z
    zz = (-normal[0] * xx - normal[1] * yy - d) * 1. / normal[2]

    ax2.plot_surface(xx, yy, zz)

    # plt.show()


test_orthogonal_regression()
