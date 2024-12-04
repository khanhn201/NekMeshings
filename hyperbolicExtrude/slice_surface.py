import numpy as np


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_intersect_points(element, x_position):
    slice_points = []
    edges = [(element[0], element[1]), (element[1], element[2]), (element[2], element[0])];
    
    for p_start, p_end in edges:
        if (p_start[0] - x_position) * (p_end[0] - x_position) < 0:
            t = (x_position - p_start[0]) / (p_end[0] - p_start[0])
            y = p_start[1] + t * (p_end[1] - p_start[1])
            z = p_start[2] + t * (p_end[2] - p_start[2])
            slice_points.append([x_position, y, z])
    
    return slice_points

def check_orientation_2d(points):
    area = 0
    for i in range(len(points)):
        p1 = points[i]
        p2 = points[(i + 1) % len(points)]  # Wrap around to the first point
        area += (p2[1] - p1[1]) * (p2[2] + p1[2])  # Cross product in 2D

    return area > 0


def ordered_points_from_edges(slice_points, tolerance=1e-4):
    if len(slice_points) == 0:
        return np.array([])

    def find_or_create_key(point, point_map):
        for key in point_map.keys():
            if np.linalg.norm(np.array(point) - np.array(key)) <= tolerance:
                return key
        point_map[tuple(point)] = True
        return tuple(point)

    connectivity = {}
    point_map = {}

    for points in slice_points:
        if len(points) == 2:
            p1 = find_or_create_key(points[0], point_map)
            p2 = find_or_create_key(points[1], point_map)
            connectivity.setdefault(p1, []).append(p2)
            connectivity.setdefault(p2, []).append(p1)

    ordered = [list(connectivity.keys())[0]]
    neighbors = connectivity[list(connectivity.keys())[0]]
    current = [p for p in neighbors if p not in ordered][0]
    ordered.append(current)

    # Traverse through the edges
    while True:
        neighbors = connectivity[current]
        neighbors_next = [p for p in neighbors if p not in ordered]
        if len(neighbors_next) < 1: # closed loop
            break
        next_point = neighbors_next[0]
        ordered.append(next_point)
        current = next_point
    
    return np.array(ordered)

def slice_surface(elements, slice_positions):
    slices = [[] for _ in slice_positions]

    for element in elements:
        for i in range(len(slice_positions)):
            intersect_points = get_intersect_points(element, slice_positions[i])
            if len(intersect_points) == 2:
                slices[i].append(intersect_points)
    np_slices = []
    for i in range(len(slice_positions)):
        edges = slices[i]
        ordered_points = ordered_points_from_edges(edges)
        if not check_orientation_2d(ordered_points):
            ordered_points = ordered_points[::-1]
        np_slices.append(np.array(ordered_points))

    return np_slices

