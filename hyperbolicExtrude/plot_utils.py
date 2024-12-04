import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_slices_3d(slices):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for slice_points in slices:
        if len(slice_points)>0:  # Only plot non-empty slices
            slice_points = slice_points
            ax.plot(slice_points[:, 0], slice_points[:, 1], slice_points[:, 2])

    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")

    def set_aspect_equal_3d(ax):
        """Set 3D plot aspect ratio to 1:1:1."""
        extents = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
        centers = np.mean(extents, axis=1)
        ranges = np.max(extents, axis=1) - np.min(extents, axis=1)
        max_range = max(ranges)
        for ctr, rng in zip(centers, ranges):
            ax.set_xlim3d([centers[0] - max_range/2, centers[0] + max_range/2])
            ax.set_ylim3d([centers[1] - max_range/2, centers[1] + max_range/2])
            ax.set_zlim3d([centers[2] - max_range/2, centers[2] + max_range/2])

    set_aspect_equal_3d(ax)
    plt.show()
def plot_slices_sidebyside(slices):
    fig, axes = plt.subplots(8, 2, figsize=(20, 30))
    y_values = []
    for slice_list in slices:
        if len(slice_list) > 0:
            slice_points = slice_list
            y_values.extend(slice_points[:, 1])
    y_min = min(y_values)
    y_max = max(y_values)

    for i in range(len(slices)):
        row = i // 2
        col = i % 2
        ax = axes[row, col]

        if len(slices[i]) > 0:
            slice_points = slices[i]
            ax.plot(slice_points[:, 1], slice_points[:, 2])
            ax.set_title(f"x = {slice_points[0, 0]:.2f}")
        ax.set_xlabel("Y")
        ax.set_ylabel("Z")
        ax.set_aspect('equal')
        ax.grid(True)
        
        ax.set_xlim([y_min, y_max])

    plt.show()


