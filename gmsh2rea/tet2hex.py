import numpy as np
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
            return True
    return False


def tet2hex(elements, boundaries):
    count_left_hand = 0
    hex_elements = []
    hex_boundaries = []
    n_elements = 0
    for i, element in enumerate(elements):
        v1, v2, v3, v4 = element

        m12 = (v1 + v2) / 2
        m13 = (v1 + v3) / 2
        m14 = (v1 + v4) / 2
        m23 = (v2 + v3) / 2
        m24 = (v2 + v4) / 2
        m34 = (v3 + v4) / 2
        m123 = (v1 + v2 + v3) / 3
        m124 = (v1 + v2 + v4) / 3
        m134 = (v1 + v3 + v4) / 3
        m234 = (v2 + v3 + v4) / 3
        m1234 = (v1 + v2 + v3 + v4) / 4

        hex1 = [v1, m12, m123, m13, m14, m124, m1234, m134]
        hex2 = [m12, v2, m23, m123, m124, m24, m234, m1234]
        hex3 = [m23, v3, m13, m123, m234, m34, m134, m1234]
        hex4 = [m124, m24, m234, m1234, m14, v4, m34, m134]

        if is_left_handed(hex1):
            count_left_hand += 1
        if is_left_handed(hex2):
            count_left_hand += 1
        if is_left_handed(hex3):
            count_left_hand += 1
        if is_left_handed(hex4):
            count_left_hand += 1
        bc = [(item[1], item[2]) for item in boundaries if item[0] == i]
        for face, tag in bc:
            if face == 0:
                hex_boundaries.append((n_elements + 0, 5-1, tag))
                hex_boundaries.append((n_elements + 1, 5-1, tag))
                hex_boundaries.append((n_elements + 2, 5-1, tag))
            if face == 1:
                hex_boundaries.append((n_elements + 0, 4-1, tag))
                hex_boundaries.append((n_elements + 2, 2-1, tag))
                hex_boundaries.append((n_elements + 3, 6-1, tag))
            if face == 2:
                hex_boundaries.append((n_elements + 0, 1-1, tag))
                hex_boundaries.append((n_elements + 1, 1-1, tag))
                hex_boundaries.append((n_elements + 3, 1-1, tag))
            if face == 3:
                hex_boundaries.append((n_elements + 1, 2-1, tag))
                hex_boundaries.append((n_elements + 2, 1-1, tag))
                hex_boundaries.append((n_elements + 3, 2-1, tag))
                

        hex_elements.extend([hex1, hex2, hex3, hex4])
        n_elements += 4


    assert count_left_hand == 0, f'found {count_left_hand} left-handed elements, cannot fix left-handed elements yet'
    hex_elements = np.array(hex_elements)
    print(hex_elements.shape)
    print(len(hex_boundaries))
    return hex_elements, hex_boundaries

