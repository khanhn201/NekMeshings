import numpy as np
import numpy.linalg as la
import scipy as sp
import matplotlib.pyplot as plt

from scipy.interpolate import splprep, splev

    
def angular_force(s1, s2, s3, k_theta=1.0):
    s1_norm = np.linalg.norm(s1)
    s2_norm = np.linalg.norm(s2)
    s3_norm = np.linalg.norm(s3)
    
    cos_theta_1 = np.dot(s3, s1) / (s3_norm * s1_norm)
    cos_theta_2 = np.dot(s3, s2) / (s3_norm * s2_norm)
    
    cos_theta_1 = np.clip(cos_theta_1, -1.0, 1.0)
    cos_theta_2 = np.clip(cos_theta_2, -1.0, 1.0)
    
    theta_1 = np.arccos(cos_theta_1)
    theta_2 = np.arccos(cos_theta_2)
    
    force_spring_1 = k_theta * theta_1 * (s1 / s1_norm)  # Force on s3 due to s1
    force_spring_2 = k_theta * theta_2 * (s2 / s2_norm)  # Force on s3 due to s2
    
    net_force = force_spring_1 + force_spring_2
    return net_force


def splinefit(points):
    min_index = np.argmin(points[:, 1])
    points = np.roll(points, -min_index, axis=0)
    points = np.vstack([points, points[0]])
    max_index = np.argmax(points[:, 1])
    tck, u = splprep([points[:, 1], points[:, 2]], s=0, per=True)
    return tck, u[max_index]

def extrude_slice_gordonhall(slice):
    """
        slice: numpy array of shape (N, 3)
        - ordered series of points along the slice
        - not closed
    """
    n_segments = 20
    k_max = 15
    x = slice[0,0]
    tck, param_at_max = splinefit(slice)
    s_fine = np.concatenate([
        np.linspace(0, param_at_max, n_segments // 2, endpoint=False),
        np.linspace(param_at_max, 1, n_segments // 2, endpoint=False),
    ])  # y_min at 0, y_max at param_at_max
    y_fine, z_fine = splev(s_fine, tck)
    y_at_max, z_at_max = splev(param_at_max, tck)
    
    R_max = 2000

    all_layers_y = [y_fine]
    all_layers_z = [z_fine]
    angle_offset = np.arctan2(z_fine[0], -y_fine[0])
    angle_offset_max = np.arctan2(z_at_max, -y_at_max)
    if angle_offset_max < 0:
        angle_offset_max += 2*np.pi
    farfield_points = np.zeros((len(s_fine),2))
    for t in range(n_segments):
        if t < n_segments//2:
            farfield_points[t] = np.array([
                -R_max*np.cos((angle_offset_max-angle_offset)*t/(n_segments//2) + angle_offset),
                R_max*np.sin((angle_offset_max-angle_offset)*t/(n_segments//2) + angle_offset)
            ])
        else:
            farfield_points[t] = np.array([
                -R_max*np.cos((2*np.pi-angle_offset_max+angle_offset)*(t-(n_segments//2))/(n_segments//2) + angle_offset_max),
                R_max*np.sin((2*np.pi-angle_offset_max+angle_offset)*(t-(n_segments//2))/(n_segments//2) + angle_offset_max)
            ])
    for k in range(k_max):
        s = (k+1)/k_max
        y_fine_next = []
        z_fine_next = []
        for t in range(n_segments):
            U = (1-s)*np.array([
                y_fine[t], z_fine[t]
            ]) + s*farfield_points[t]
            y_fine_next.append(U[0])
            z_fine_next.append(U[1])
        all_layers_y.append(y_fine_next)
        all_layers_z.append(z_fine_next)


    n_smoothing = 100
    for i in range(n_smoothing):
        for k in range(1, k_max - 1):
            for t in range(n_segments):
                if t == 0 or t == n_segments//2:
                    continue
                next_t = (t + 1)%n_segments
                prev_t = t - 1
                s1 = np.array([
                    all_layers_y[k-1][prev_t] - all_layers_y[k-1][t],
                    all_layers_z[k-1][prev_t] - all_layers_z[k-1][t],
                ])
                s2 = np.array([
                    all_layers_y[k-1][next_t] - all_layers_y[k-1][t],
                    all_layers_z[k-1][next_t] - all_layers_z[k-1][t],
                ])
                s3 = np.array([
                    all_layers_y[k][t] - all_layers_y[k-1][t],
                    all_layers_z[k][t] - all_layers_z[k-1][t],
                ])
                force_on_s3 = angular_force(s1, s2, s3)
                torque = np.cross(s3, force_on_s3)
                delta_theta = torque * 0.00005
                rotation_matrix = np.array([
                    [np.cos(delta_theta), -np.sin(delta_theta)],
                    [np.sin(delta_theta), np.cos(delta_theta)]
                ])
                s3_rotated = np.dot(rotation_matrix, s3)
                s3_new = np.dot(s3, s3_rotated) / np.dot(s3_rotated, s3_rotated) * s3_rotated
                all_layers_y[k][t] = all_layers_y[k-1][t] + s3_new[0]
                all_layers_z[k][t] = all_layers_z[k-1][t] + s3_new[1]

    # Plot result
    # for layer_y, layer_z in zip(all_layers_y, all_layers_z):
    #
    #     layer_y_closed = np.append(layer_y, layer_y[0])
    #     layer_z_closed = np.append(layer_z, layer_z[0])
    #     plt.plot(layer_y_closed, layer_z_closed, 'k-', linewidth=0.5)
    # for i in range(len(all_layers_y) - 1):
    #     layer_y = all_layers_y[i]
    #     layer_z = all_layers_z[i]
    #     next_layer_y = all_layers_y[i + 1]
    #     next_layer_z = all_layers_z[i + 1]
    #
    #     for j in range(len(layer_y)):
    #         plt.plot(
    #             [layer_y[j], next_layer_y[j]],
    #             [layer_z[j], next_layer_z[j]],
    #             'k-',
    #             linewidth=0.5
    #         )
    # plt.title("Multiple Layer Extrusion")
    # plt.xlabel("Y")
    # plt.ylabel("Z")
    # plt.axis('equal')
    # plt.show()

    # Assemble the final extruded array
    stacked_array = np.zeros((k_max, n_segments, 3))
    for k in range(k_max):
        y_layer = np.array(all_layers_y[k])
        z_layer = np.array(all_layers_z[k])
        x_layer = np.full_like(y_layer, x) 
        stacked_array[k, :, 0] = x_layer
        stacked_array[k, :, 1] = y_layer
        stacked_array[k, :, 2] = z_layer
    return stacked_array



