"""
Generation of Body-Fitted Coordinates Using Hyperbolic Partial Differential Equations
Joseph L. Steger and Denny S. Chaussee
https://doi.org/10.1137/0901031

Using cell volume scheme
Differences in notation:
- eta parameterize the surface
- zeta is the march direction perpendicular to the surface
"""



import numpy as np
import numpy.linalg as la
import scipy as sp
import matplotlib.pyplot as plt

from scipy.interpolate import splprep, splev

def splinefit(points):
    min_index = np.argmin(points[:, 1])
    points = np.roll(points, -min_index, axis=0)
    points = np.vstack([points, points[0]])
    max_index = np.argmax(points[:, 1])
    tck, u = splprep([points[:, 1], points[:, 2]], s=0, per=True)
    return tck, u[max_index]


def extrude_slice_2d(slice):
    """
        slice: numpy array of shape (N, 3)
        - ordered series of points along the slice
        - not closed
    """
    n_segments = 100
    k_max = 20
    x = slice[0,0]
    tck, param_at_max = splinefit(slice)
    s_fine = np.concatenate([
        np.linspace(0, param_at_max, n_segments // 2, endpoint=False),
        np.linspace(param_at_max, 1, n_segments // 2, endpoint=False),
    ])  # y_min at 0, y_max at param_at_max
    y_fine, z_fine = splev(s_fine, tck)

    all_layers_y = [y_fine]
    all_layers_z = [z_fine]
    prevV = np.zeros((n_segments,))
    for k in range(k_max):
        y_fine = all_layers_y[-1]
        z_fine = all_layers_z[-1]
        big_mat = np.zeros((n_segments*2, n_segments*2))
        big_g = np.zeros((n_segments*2,))
        for j in range(n_segments):
            zetap1y, zetap1z = y_fine[0], z_fine[0]
            if j < n_segments - 1:
                zetap1y, zetap1z = y_fine[j+1], z_fine[j+1]
            zetam1y, zetam1z = y_fine[-1], z_fine[-1]
            if j > 0:
                zetam1y, zetam1z = y_fine[j-1], z_fine[j-1]
            xij = np.array([y_fine[j], z_fine[j]])
            xijp1 = np.array([zetap1y, zetap1z])
            xijm1 = np.array([zetam1y, zetam1z])
            
            length1 = np.linalg.norm(xijp1 - xij)
            length2 = np.linalg.norm(xij - xijm1)
            V = (length1 + length2) / 2*10
            x_eta = (xijp1-xijm1)/2
            x_zeta = -np.array([
                x_eta[1]*prevV[j]/la.norm(x_eta)**2,
                x_eta[0]*prevV[j]/la.norm(x_eta)**2
            ])
            A = np.array([
                x_zeta,
                [x_zeta[1], -x_zeta[0]]
            ])
            Bi = la.inv(np.array([
                x_eta,
                [-x_eta[1], x_eta[0]]
            ]))
            if not sp.linalg.issymmetric(Bi@A, rtol=1E-6):
                print(Bi@A)
            f = np.array([
                0, 
                V + prevV[j]
            ])
            if k == 0:
                f = np.array([
                    0, 
                    V + V
                ])
            C = Bi@A

            prevV[j] = V

            big_mat[
                j*2:j*2 + 2,
                j*2:j*2 + 2
            ] = np.eye(2)
            if j > 0:
                big_mat[
                    j*2:j*2 + 2,
                    (j-1)*2:(j-1)*2 + 2
                ] = -C/2
            else:
                big_mat[
                    j*2:j*2 + 2,
                    (n_segments-1)*2:(n_segments-1)*2 + 2
                ] = -C/2
            if j < n_segments - 1:
                big_mat[
                    j*2:j*2 + 2,
                    (j+1)*2:(j+1)*2 + 2
                ] = C/2
            else:
                big_mat[
                    j*2:j*2 + 2,
                    0*2:0*2 + 2
                ] = C/2
            big_g[j*2: j*2 + 2] = Bi@f + xij

        # if k == 1:
        #     print(big_mat)
        #     print(list(big_mat[0]))
        #     print(big_g)
        y_fine_next = np.zeros_like(y_fine)
        z_fine_next = np.zeros_like(z_fine)
        delta = la.solve(big_mat, big_g)
        for j in range(n_segments):
            delta_seg = delta[j*2  : j*2  + 2]
            y_fine_next[j] = delta_seg[0]
            z_fine_next[j] = delta_seg[1]
        all_layers_y.append(y_fine_next)
        all_layers_z.append(z_fine_next)



    # Plot result
    for layer_y, layer_z in zip(all_layers_y, all_layers_z):
        plt.plot(layer_y, layer_z, 'k-', linewidth=0.5)
    for i in range(len(all_layers_y) - 1):
        layer_y = all_layers_y[i]
        layer_z = all_layers_z[i]
        next_layer_y = all_layers_y[i + 1]
        next_layer_z = all_layers_z[i + 1]
        
        for j in range(len(layer_y)):
            plt.plot(
                [layer_y[j], next_layer_y[j]],
                [layer_z[j], next_layer_z[j]],
                'k-',
                linewidth=0.5
            )
    plt.title("Multiple Layer Extrusion")
    plt.xlabel("Y")
    plt.ylabel("Z")
    plt.axis('equal')
    plt.show()



