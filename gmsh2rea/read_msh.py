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

physical_groups = gmsh.model.getPhysicalGroups()
physical_3d_groups = [group for group in physical_groups if group[0] == 3]
assert len(physical_3d_groups) == 1, f'expected one 3D physical group, got {len(physical_3d_groups)}'
physical_2d_groups = [group for group in physical_groups if group[0] == 2]
print(f'found {len(physical_2d_groups)} surface groups')




# Start
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
node_dict = {tag: node_coords[3*i:3*(i+1)] for i, tag in enumerate(node_tags)}

_, element_tags, node_tags = gmsh.model.mesh.getElements(dim=3)
element_tags = element_tags[0]
node_tags = node_tags[0]


face_nodes_to_physical_2d = {}
for group in physical_2d_groups:
    dim, tag = group
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    for entity in entities:
        _, entity_element_tags, entity_element_nodes = gmsh.model.mesh.getElements(dim, entity)
        entity_element_tags = entity_element_tags[0]
        entity_element_nodes = entity_element_nodes[0]
        for i in range(len(entity_element_tags)):
            nodes_in_faces = entity_element_nodes[8*i:8*i+4]
            nodes_in_faces.sort()
            face_nodes_to_physical_2d[tuple(nodes_in_faces)] = tag


n_element = 0
count_left_hand = 0
elements = []
boundary_conditions = []

for i, tags in enumerate(element_tags):
    n_element += 1
    nodes = node_tags[20*i:20*i+8]
    node_coords = [node_dict[nodeid] for nodeid in nodes]
    if is_left_handed(node_coords):
        count_left_hand += 1
        # plotElement(nodeCoords)
    elements.append(node_coords)
    
    # Check if face is boundary condition
    for face_num, face in enumerate(faces):
        face_nodes = [nodes[node] for node in face]
        face_nodes.sort()
        tag = face_nodes_to_physical_2d.get(tuple(face_nodes))
        if tag:
            boundary_conditions.append((n_element, face_num, tag))


elements = np.array(elements)
print(elements.shape)
print(f'found {count_left_hand} left-handed elements')

print(f'found {len(boundary_conditions)} surface boundary conditions')
boundaries = np.zeros((n_element, 6))
for element, face, tag in boundary_conditions:
    boundaries[element, face] = tag
gmsh.finalize()
export_rea('output.rea', elements, boundaries)
