import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection



def plot_element(node_coords):
    faces = [
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 4, 8, 7],
        [4, 1, 5, 8],
        [1, 4, 3, 2],
        [5, 6, 7, 8]
    ]
    faces = [[node - 1 for node in face] for face in faces]

    
    center = [sum(coords) / len(coords) for coords in zip(*node_coords)]
    centered_coords = [[x - center[0], y - center[1], z - center[2]] for x, y, z in node_coords]
    # centered_coords = nodeCoords

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    poly3d = [[centered_coords[vertex] for vertex in face] for face in faces]
    poly3d_collection = Poly3DCollection(poly3d, alpha=0.5, edgecolor='k')
    ax.add_collection3d(poly3d_collection)

    x, y, z = zip(*centered_coords)
    ax.scatter(x, y, z, s=50)
    
    for i, (xi, yi, zi) in enumerate(centered_coords):
        ax.text(xi, yi, zi, f"{i+1}", color="k", fontsize=12, ha='center')
    plt.show()

def plot_tet_elements(tet_elements, tet_boundaries):
    faces = [
        [1, 3, 2],
        [1, 4, 3],
        [1, 2, 4],
        [2, 3, 4]
    ]
    faces = [[node - 1 for node in face] for face in faces]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Define some colors for different boundary tags
    boundary_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
    
    # for hex_element in hex_elements:
    #     vertices = [hex_element[i] for i in range(8)]
    #     element_faces = [
    #         [vertices[i] for i in face] for face in faces
    #     ]
    #     poly3d = [[list(vertex) for vertex in face] for face in element_faces]
    #     ax.add_collection3d(Poly3DCollection(poly3d, facecolors='grey', linewidths=0.1, edgecolors='black', alpha=0.3))
    
    # Highlight boundaries with specific colors
    for (element_idx, face, tag) in tet_boundaries:
        color = boundary_colors[tag % len(boundary_colors)]
        tet_element = tet_elements[element_idx]
        vertices = [tet_element[i] for i in range(4)]
        
        # Choose the face based on the index
        boundary_faces = [
            [vertices[i] for i in boundary] for boundary in faces
        ]
        boundary_face = boundary_faces[face]
        
        # Add boundary face to plot
        poly = Poly3DCollection([boundary_face], facecolors=color, edgecolors='k', alpha=0.8)
        poly.set_alpha(0.2)
        ax.add_collection3d(poly)

        
    all_vertices = np.concatenate(tet_elements)
    x_min, y_min, z_min = np.min(all_vertices, axis=0)
    x_max, y_max, z_max = np.max(all_vertices, axis=0)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])
    # Show the plot
    plt.show()


def plot_hex_elements(hex_elements, hex_boundaries):
    faces = [
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 4, 8, 7],
        [4, 1, 5, 8],
        [1, 4, 3, 2],
        [5, 6, 7, 8]
    ]
    faces = [[node - 1 for node in face] for face in faces]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Define some colors for different boundary tags
    face_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
    boundary_colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
    
    
    # for element_idx, hex_element in enumerate(hex_elements):
    #     vertices = [hex_element[i] for i in range(8)]
    #     element_faces = [
    #         [vertices[i] for i in face] for face in faces
    #     ]
    #     poly3d = [[list(vertex) for vertex in face] for face in element_faces]
    #     # ax.add_collection3d(Poly3DCollection(poly3d, facecolors='grey', linewidths=0.1, edgecolors='black', alpha=0.3))
    #
    #     centroid = np.mean(vertices, axis=0)
    #     ax.text(centroid[0], centroid[1], centroid[2], f"E{element_idx}", color="blue", fontsize=10)
    #
    #     # Label each vertex
    #     for vertex_idx, vertex in enumerate(vertices):
    #         ax.text(vertex[0], vertex[1], vertex[2], f"{element_idx}V{vertex_idx}", color="black", fontsize=8)
    #
    #     for i in range(6):
    #         color = face_colors[i % len(face_colors)]
    #         hex_element = hex_elements[element_idx]
    #         vertices = [hex_element[i] for i in range(8)]
    #
    #         # Use the appropriate face from the face list
    #         boundary_face_vertices = [vertices[i] for i in faces[i]]
    #         poly = Poly3DCollection([boundary_face_vertices], facecolors=color, edgecolors='k')
    #         poly.set_alpha(0.8)
    #         ax.add_collection3d(poly)
    
    # Highlight boundaries with specific colors
    for (element_idx, face, tag) in hex_boundaries:
        if tag != 4:
            continue
        color = boundary_colors[tag % len(boundary_colors)]
        hex_element = hex_elements[element_idx]
        vertices = [hex_element[i] for i in range(8)]

        # Choose the face based on the index
        boundary_faces = [
            [vertices[i] for i in boundary] for boundary in faces
        ]
        boundary_face = boundary_faces[face]

        # Add boundary face to plot
        poly = Poly3DCollection([boundary_face], facecolors=color, edgecolors='k', alpha=0.8)
        poly.set_alpha(0.2)
        ax.add_collection3d(poly)

        
    all_vertices = np.concatenate(hex_elements)
    x_min, y_min, z_min = np.min(all_vertices, axis=0)
    x_max, y_max, z_max = np.max(all_vertices, axis=0)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])
    # Show the plot
    plt.show()
    
def plot_faces(faces_with_tag1):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    coords = []
    # Loop through each face
    for face in faces_with_tag1:
        # Convert to a numpy array for easier manipulation
        face_coords = np.array(face)
        
        # Add a polygon for the face
        poly = Poly3DCollection([face_coords], alpha=0.7, edgecolor="black")
        poly.set_alpha(0.5)
        ax.add_collection3d(poly)
        coords.append(face_coords)
        
    # Set labels and aspect
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
    all_vertices = np.concatenate(faces_with_tag1)
    x_min, y_min, z_min = np.min(all_vertices, axis=0)
    x_max, y_max, z_max = np.max(all_vertices, axis=0)
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])


    
    plt.show()

