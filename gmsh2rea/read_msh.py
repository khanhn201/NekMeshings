import gmsh
import numpy as np

from export_rea import export_rea
# from plot_element import plot_element

faces = [
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 4, 8, 7],
    [4, 1, 5, 8],
    [1, 4, 3, 2],
    [5, 6, 7, 8]
]
faces = [[node - 1 for node in face] for face in faces]

def is_left_handed(node_coords):
    v1, v2, v3, v4, v5, v6, v7, v8 = [np.array(node) for node in node_coords]
    
    corners = [
        (v1, v2, v4, v5),
        (v2, v3, v1, v6),
        (v3, v4, v2, v7),
        (v4, v1, v3, v8),
        (v5, v8, v6, v1),
        (v6, v5, v7, v2),
        (v7, v6, v8, v3),
        (v8, v7, v5, v4),
    ]
    for corner in corners:
        vv1, vv2, vv3, vv4 = corner
        vector1 = vv2 - vv1
        vector2 = vv3 - vv1
        vector3 = vv4 - vv1
        cross_product = np.cross(vector1, vector2)

        if np.dot(cross_product, vector3) <= 0:
            return True  # Left-handed
    return False



msh_file = "./turbine.msh"

gmsh.initialize()
gmsh.open(msh_file)

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_dict = {tag: node_coords[3*i:3*(i+1)] for i, tag in enumerate(node_tags)}

_, element_tags, node_tags = gmsh.model.mesh.getElements(dim=3)
element_tags = element_tags[0]
node_tags = node_tags[0]

print(gmsh.model.getPhysicalGroups())

_, surface_element_tags, surface_node_tags = gmsh.model.mesh.getElements(dim=2)
surface_element_tags = surface_element_tags[0]
surface_node_tags = surface_node_tags[0]

print(element_tags, len(element_tags))
print(node_tags, len(node_tags))
print(surface_element_tags, len(surface_element_tags))
print(surface_node_tags, len(surface_node_tags))

        
n_element = 0
count_left_hand = 0
elements = []
surface_faces = []
for i, tags in enumerate(element_tags):
    n_element += 1
    nodes = node_tags[20*i:20*i+8]
    node_coords = [node_dict[nodeid] for nodeid in nodes]
    if is_left_handed(node_coords):
        count_left_hand += 1
        # plotElement(nodeCoords)
        
    elements.append(node_coords)

    for face_num, face in enumerate(faces):
        face_nodes = [nodes[node] for node in face]
        if all(node in surface_node_tags for node in face_nodes):
            surface_faces.append((n_element, face_num))

elements = np.array(elements)
print(elements.shape)
print(len(surface_faces))
print(f'found {count_left_hand} left-handed elements')
gmsh.finalize()
export_rea('output.rea', elements)
