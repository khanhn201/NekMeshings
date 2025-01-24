import numpy as np


# from mesh_surface import mesh_surface
from plot_utils import plot_slices_sidebyside, plot_slices_3d
from slice_surface import slice_surface
from extrude_slice import extrude_slice
from extrude_slice_2d import extrude_slice_2d
from extrude_slice_gordonhall import extrude_slice_gordonhall

from plot_element import plot_hex_elements
from export_rea import export_rea


print("Meshing surface")
# elements = mesh_surface('./nrelVI.step', 10)  # Rect mesh of the surface
# np.save("elements.npy", elements)
elements = np.load("elements.npy")
min_vals, max_vals = np.min(elements.reshape(-1, 3), axis=0), np.max(elements.reshape(-1, 3), axis=0)
print("Min values (x, y, z):", min_vals)
print("Max values (x, y, z):", max_vals)
x_min, x_max = min_vals[0], max_vals[0]

print("Slicing surface along x axis")
slice_count = 208
# slice_count = 32
slice_positions = np.concatenate([
    np.linspace(-5000, -1250, slice_count//4, endpoint=False),
    np.linspace(-1250, -500, (slice_count//16)*3, endpoint=False),
    np.linspace(-500, 500, slice_count//8, endpoint=False),
    np.linspace(500, 1250, (slice_count//16)*3, endpoint=False),
    np.linspace(1250, 5000, slice_count//4)
])
slices = slice_surface(elements, slice_positions)
print(slices[0][0])
print(slices[-1][0])
plot_slices_3d(slices)

for i in range(len(slice_positions)):
    print(slices[i].shape)
    x = slices[i][0]
    output_file = f'./slices/slice{i:04}.txt'

    with open(output_file, 'w') as f:
        f.write('x y z\n')
        for row in slices[i]:
            f.write(f"{row[0]} {row[1]} {row[2]}\n")

# print("Extruding slices")
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
        if dot_prod <= 0:
            print('Left-handed/ill-shaped element detected!')
            return False



extruded_volume = np.load("extruded_volume.npy")
print(extruded_volume.shape)
# print(extruded_volume)
elements = []
boundaries = []
count_elements = 0
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
            if j == 0:
                boundaries.append((count_elements, 1-1, 1)) # Wall
            if i == 0:
                boundaries.append((count_elements, 5-1, 4)) # Outflow
            if i == extruded_volume.shape[0] - 2:
                boundaries.append((count_elements, 6-1, 4))
            if j == extruded_volume.shape[1] - 2:
                boundaries.append((count_elements, 3-1, 4))
            count_elements += 1
elements = np.array(elements)

plot_hex_elements(elements, boundaries, 1) # Outflow
plot_hex_elements(elements, boundaries, 4) # Wall
export_rea('turb_hyp.rea', elements, boundaries)

