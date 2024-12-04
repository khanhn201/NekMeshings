import csv

from read_msh import read_msh
from tet2hex import tet2hex
from export_rea import export_rea
from plot_element import plot_hex_elements, plot_tet_elements, plot_element, plot_faces

msh_file = "./nreltet.msh"
elements, boundaries = read_msh(msh_file)
# plot_tet_elements(elements, boundaries)
elements, boundaries = tet2hex(elements, boundaries)
# plot_hex_elements(elements, boundaries)


def export_faces_to_csv(faces_with_tag, filename):
    header = [
        "index",
        "x1", "y1", "z1",
        "x2", "y2", "z2",
        "x3", "y3", "z3",
        "x4", "y4", "z4"
    ]
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for i, face in enumerate(faces_with_tag):
            row = [i]
            for vertex in face:
                row.extend(vertex)
            writer.writerow(row)
faces = [
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 4, 8, 7],
    [4, 1, 5, 8],
    [1, 4, 3, 2],
    [5, 6, 7, 8]
]
faces = [[node - 1 for node in face] for face in faces]
def get_faces_with_tag(hex_elements, hex_boundaries, tag):
    faces_with_tag1 = []
    for element_id, face_index, face_tag in hex_boundaries:
        if face_tag == tag:
            face_nodes = faces[face_index]
            face_coords = [hex_elements[element_id][node] for node in face_nodes]
            faces_with_tag1.append(face_coords)
    return faces_with_tag1

faces_with_tag = get_faces_with_tag(elements, boundaries, tag=4)
# plot_faces(faces_with_tag)
print(len(faces_with_tag))

# export_faces_to_csv(faces_with_tag, "turbineBladeWall")

# export_rea('turbine.rea', elements, boundaries)
