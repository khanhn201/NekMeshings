import gmsh
import numpy as np

def mesh_surface(file, mesh_size):
    gmsh.initialize()

    v = gmsh.model.occ.importShapes(file)

    gmsh.model.occ.remove([(v[0][0], v[0][1])], True) # Remove external cylinder
    gmsh.model.occ.rotate([(v[1][0], v[1][1])], 0, 0, 0, 0, 0, 1, -np.pi / 4) # Rotate blade from pi/2 to horizontal

    gmsh.model.occ.synchronize()
    print(gmsh.model.getEntities(dim=3))
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.model.mesh.generate(2)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_dict = {tag: node_coords[3*i:3*(i+1)] for i, tag in enumerate(node_tags)}

    element_types, element_tags, node_tags = gmsh.model.mesh.getElements(dim=2)
    element_types = element_types[0]
    element_tags = element_tags[0]
    node_tags = node_tags[0]
    print(node_tags, len(node_tags))
    print(element_tags, len(element_tags))

    elements = []
    for i in range(len(element_tags)):
        nodes = node_tags[3*i:3*i+3]
        node_coords = [node_dict[nodeid] for nodeid in nodes]
        elements.append(node_coords)
    elements = np.array(elements)
    print(elements.shape)
    gmsh.finalize()
    return elements

