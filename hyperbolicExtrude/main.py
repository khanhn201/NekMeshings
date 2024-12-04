import numpy as np


from mesh_surface import mesh_surface
from plot_utils import plot_slices_sidebyside, plot_slices_3d
from slice_surface import slice_surface
from extrude_slice import extrude_slice
from extrude_slice_2d import extrude_slice_2d
from extrude_slice_gordonhall import extrude_slice_gordonhall


print("Meshing surface")
# elements = mesh_surface('./nrelVI.step', 10)
# np.save("elements.npy", elements)
# elements = np.load("elements.npy")
# min_vals, max_vals = np.min(elements.reshape(-1, 3), axis=0), np.max(elements.reshape(-1, 3), axis=0)
# print("Min values (x, y, z):", min_vals)
# print("Max values (x, y, z):", max_vals)
# x_min, x_max = min_vals[0], max_vals[0]
#
# print("Slicing surface along x axis")
# # slice_positions = np.linspace(x_min+0.1, x_max-0.1, n_slices)
#
# slice_positions = np.concatenate([
#     np.linspace(-5000, -1250, 15),
#     np.linspace(-1250, -500, 8),
#     np.linspace(-500, 500, 6),
#     np.linspace(500, 1250, 8),
#     np.linspace(1250, 5000, 15)
# ])
# slices = slice_surface(elements, slice_positions)
#
# # plot_slices_sidebyside(slices)
# plot_slices_3d(slices)
# # np.save("slice5.npy", slices[5])
# # np.save("slice6.npy", slices[6])
# np.save("slice9.npy", slices[9])
#
#
print("Extruding slices")
# slice = np.load("slice7.npy")
# slice = np.load("slice6.npy")
# slice = np.load("slice5.npy")
# extrude_slice_2d(slice)
# extrude_slice(slice)
# extrude_slice_gordonhall(slice)
# slice = np.load("slice9.npy")
# extruded = extrude_slice_gordonhall(slice)

# extruded_volume = []
# for i in range(len(slice_positions)):
#     extruded = extrude_slice_gordonhall(slices[i])
#     extruded_volume.append(extruded)
# extruded_volume = np.array(extruded_volume)
# np.save("extruded_volume.npy", extruded_volume)


def check_lefthand(element):
    corner_list = [
        [0, 1, 3, 4],
        [1, 2, 0, 5],
        [2, 3, 1, 6],
        [3, 0, 2, 7],
        [4, 7, 5, 0],
        [5, 4, 6, 1],
        [6, 5, 7, 2],
        [7, 6, 4, 3]
    ]
    for face in corner_list:
        node1, node2, node4, node5 = element[face[0]], element[face[1]], element[face[2]], element[face[3]]
        vec12 = node2 - node1
        vec14 = node4 - node1
        vec15 = node5 - node1
        cross_vec = np.cross(vec12, vec14)
        dot_prod = np.dot(cross_vec, vec15)
        if dot_prod < 0:
            print('Left-handed element detected!')
            return False



extruded_volume = np.load("extruded_volume.npy")
print(extruded_volume.shape)
# print(extruded_volume)
elements = []
for i in range(extruded_volume.shape[0]-1):
    for j in range(extruded_volume.shape[1]-1):
        for k in range(extruded_volume.shape[2]):
            k_next = k + 1
            if k == extruded_volume.shape[2]-1:
                k_next = 0
            node1 = extruded_volume[i, j, k]
            node2 = extruded_volume[i, j, k_next]
            node3 = extruded_volume[i, j+1, k_next]
            node4 = extruded_volume[i, j+1, k]
            node5 = extruded_volume[i+1, j, k]
            node6 = extruded_volume[i+1, j, k_next]
            node7 = extruded_volume[i+1, j+1, k_next]
            node8 = extruded_volume[i+1, j+1, k]
            element = np.array([
                node1, node2, node3, node4, node5, node6, node7, node8
            ])
            check_lefthand(element)
            elements.append(element)

            


