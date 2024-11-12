import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

faces = [
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 4, 8, 7],
    [4, 1, 5, 8],
    [1, 4, 3, 2],
    [5, 6, 7, 8]
]
faces = [[node - 1 for node in face] for face in faces]



def plot_element(node_coords):
    
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
