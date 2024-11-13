import gmsh
import numpy as np

faces = [
    [1, 3, 2],
    [1, 4, 3],
    [1, 2, 4],
    [2, 3, 4]
]
faces = [[node - 1 for node in face] for face in faces]


def is_left_handed(node_coords):
    v1, v2, v3, v4 = [np.array(node) for node in node_coords]
    
    corners = [
        (v1, v2, v3, v4),
        (v2, v3, v1, v4),
        (v3, v1, v2, v4),
        (v4, v1, v3, v2),
    ]
    for corner in corners:
        vv1, vv2, vv3, vv4 = corner
        vector1 = vv2 - vv1
        vector2 = vv3 - vv1
        vector3 = vv4 - vv1
        cross_product = np.cross(vector1, vector2)

        if np.dot(cross_product, vector3) <= 0:
            return True
    return False

def read_msh(msh_file):
    gmsh.initialize()
    gmsh.open(msh_file)

    physical_groups = gmsh.model.getPhysicalGroups()
    physical_3d_groups = [group for group in physical_groups if group[0] == 3]
    assert len(physical_3d_groups) == 1, f'expected one 3D physical group, got {len(physical_3d_groups)}'
    physical_2d_groups = [group for group in physical_groups if group[0] == 2]
    print(f'found {len(physical_2d_groups)} surface groups')

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_dict = {tag: node_coords[3*i:3*(i+1)] for i, tag in enumerate(node_tags)}

    element_types, element_tags, node_tags = gmsh.model.mesh.getElements(dim=3)
    assert element_types[0] == 4, f'only support tet4, found a different 3d element type'
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
                nodes_in_faces = entity_element_nodes[3*i:3*i+3]
                nodes_in_faces.sort()
                face_nodes_to_physical_2d[tuple(nodes_in_faces)] = tag


    count_elements = 0
    count_left_hand = 0
    elements = []
    boundaries = []

    for i, tags in enumerate(element_tags):
        nodes = node_tags[4*i:4*i+4]
        node_coords = [node_dict[nodeid] for nodeid in nodes]
        if is_left_handed(node_coords):
            count_left_hand += 1
        elements.append(node_coords)
        
        for face_num, face in enumerate(faces):
            face_nodes = [nodes[node] for node in face]
            face_nodes.sort()
            tag = face_nodes_to_physical_2d.get(tuple(face_nodes))
            if tag:
                boundaries.append((count_elements, face_num, tag))
        count_elements += 1


    elements = np.array(elements)
    print(f'found {elements.shape[0]} 3d elements')
    assert count_left_hand == 0, f'found {count_left_hand} left-handed elements, cannot fix left-handed elements yet'

    print(f'found {len(boundaries)} 2d surface boundary conditions')
    gmsh.finalize()
    return elements, boundaries

