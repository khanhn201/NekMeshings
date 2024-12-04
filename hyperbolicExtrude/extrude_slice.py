"""
Efficient Mesh Generation and Deformation for Aerodynamic Shape Optimization
Ney R. Secco, Gaetan K. W. Kenway, Ping He, Charles Mader and Joaquim R. R. A. Martins
https://doi.org/10.2514/1.J059491

Enhancements of a three-dimensional hyperbolic grid generation scheme
William M. Chan and Joseph L. Steger
https://doi.org/10.1016/0096-3003(92)90073-A
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

def calc_d_xi_frac(all_layers_y, all_layers_z, j, k, l):
    rjp1 = np.array([j+1, all_layers_y[l][k], all_layers_z[l][k]]) 
    rj = np.array([j, all_layers_y[l][k], all_layers_z[l][k]]) 
    rjm1 = np.array([j-1, all_layers_y[l][k], all_layers_z[l][k]]) 
    return (la.norm(rjp1 - rj) + la.norm(rjm1 - rj))

def calc_d_xi(all_layers_y, all_layers_z, j, k, l):
    return calc_d_xi_frac(all_layers_y, all_layers_z, j, k, l-1)\
            / calc_d_xi_frac(all_layers_y, all_layers_z, j, k, l)
def calc_d_eta_frac(all_layers_y, all_layers_z, j, k, l):
    kp1 = (k + 1) % len(all_layers_y[0])
    km1 = k - 1
    if km1 == -1:
        km1 += len(all_layers_y[0])
    rkp1 = np.array([j, all_layers_y[l][kp1], all_layers_z[l][kp1]]) 
    rk = np.array([j, all_layers_y[l][k], all_layers_z[l][k]]) 
    rkm1 = np.array([j, all_layers_y[l][km1], all_layers_z[l][km1]]) 
    return (la.norm(rkp1 - rk) + la.norm(rkm1 - rk))

def calc_d_eta(all_layers_y, all_layers_z, j, k, l):
    return calc_d_eta_frac(all_layers_y, all_layers_z, j, k, l-1)\
            / calc_d_eta_frac(all_layers_y, all_layers_z, j, k, l)

def calc_axi_aeta(all_layers_y, all_layers_z, j, k, l):
    rjp1 = np.array([j+1, all_layers_y[l][k], all_layers_z[l][k]]) 
    rj = np.array([j, all_layers_y[l][k], all_layers_z[l][k]]) 
    rjm1 = np.array([j-1, all_layers_y[l][k], all_layers_z[l][k]]) 
    kp1 = (k + 1) % len(all_layers_y[0])
    km1 = k - 1
    if km1 == -1:
        km1 += len(all_layers_y[0])
    rkp1 = np.array([j, all_layers_y[l][kp1], all_layers_z[l][kp1]]) 
    rk = np.array([j, all_layers_y[l][k], all_layers_z[l][k]]) 
    rkm1 = np.array([j, all_layers_y[l][km1], all_layers_z[l][km1]]) 

    rpj = rjp1 - rj
    rmj = rjm1 - rj
    rpk = rkp1 - rk
    rmk = rkm1 - rk
    rpj = rpj/la.norm(rpj)
    rmj = rmj/la.norm(rmj)
    rpk = rpk/la.norm(rpk)
    rmk = rmk/la.norm(rmk)

    nhat = np.cross(rpj - rmj, rpk - rmk)
    nhat = nhat/la.norm(nhat)
    alpha = np.arccos(np.dot(nhat, rpj))
    beta = np.arccos(np.dot(nhat, rpk))
    a_xi = 1
    if alpha < np.pi / 2:
        a_xi = 1 / (1 - np.cos(alpha)**2)
    a_eta = 1
    if beta < np.pi / 2:
        a_eta = 1 / (1 - np.cos(beta)**2)
    return a_xi, a_eta

def calc_rxixi_retaeta(all_layers_y, all_layers_z, j, k, l):
    rjp1 = np.array([j+1, all_layers_y[l][k], all_layers_z[l][k]]) 
    rj = np.array([j, all_layers_y[l][k], all_layers_z[l][k]]) 
    rjm1 = np.array([j-1, all_layers_y[l][k], all_layers_z[l][k]]) 
    kp1 = (k + 1) % len(all_layers_y[0])
    km1 = k - 1
    if km1 == -1:
        km1 += len(all_layers_y[0])
    # km1 = (k - 1) % len(all_layers_y)
    rkp1 = np.array([j, all_layers_y[l][kp1], all_layers_z[l][kp1]]) 
    rk = np.array([j, all_layers_y[l][k], all_layers_z[l][k]]) 
    rkm1 = np.array([j, all_layers_y[l][km1], all_layers_z[l][km1]]) 
    return 0.5*(rjp1 - 2*rj + rjm1), 0.5*(rkp1 - 2*rk + rkm1)
    

    

def extrude_slice(slice):
    """
        slice: numpy array of shape (N, 3)
        - ordered series of points along the slice
        - not closed
    """
    n_segments = 100
    x_count = 3
    theta = 0
    eps_i = 0.0
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
    prevV = np.zeros((x_count, n_segments))
    for k in range(k_max):
        y_fine = all_layers_y[-1]
        z_fine = all_layers_z[-1]
        big_mat = np.zeros((n_segments*x_count*3, n_segments*x_count*3))
        big_g = np.zeros((n_segments*x_count*3,))
        for i in range(x_count):
            for j in range(n_segments):
                zetap1y, zetap1z = y_fine[0], z_fine[0]
                if j < n_segments - 1:
                    zetap1y, zetap1z = y_fine[j+1], z_fine[j+1]
                zetam1y, zetam1z = y_fine[-1], z_fine[-1]
                if j > 0:
                    zetam1y, zetam1z = y_fine[j-1], z_fine[j-1]
                xij = np.array([i, y_fine[j], z_fine[j]])
                xijp1 = np.array([i, zetap1y, zetap1z])
                xijm1 = np.array([i, zetam1y, zetam1z])
                
                length1 = np.linalg.norm(xijp1 - xij)
                length2 = np.linalg.norm(xij - xijm1)
                V = (length1 + length2) / 2*10
                x_xi = np.array([1, 0, 0])
                x_eta = (xijp1-xijm1)/2
                Ci = la.inv(np.array([
                    x_xi,
                    x_eta,
                    np.cross(x_xi,x_eta)
                ]))
                # print(x_eta)
                x_zeta = Ci@np.array([0, 0, prevV[i, j]])
                if k == 0:
                    x_zeta = Ci@np.array([0, 0, V])
                prevV[i, j] = V
                A = np.array([
                    x_zeta,
                    [0, 0, 0],
                    np.cross(x_eta,x_zeta)
                ])
                B = np.array([
                    [0, 0, 0],
                    x_zeta,
                    np.cross(x_zeta,x_xi)
                ])
                # print(np.array([
                #     x_xi,
                #     x_eta,
                #     np.cross(x_xi,x_eta)
                # ]))
                # print('A', A)
                # print('B', B)
                # print('C', np.array([
                #     x_xi,
                #     x_eta,
                #     np.cross(x_xi,x_eta)
                # ]))
                P = (1+theta)/2*Ci@A - eps_i*np.eye(3)
                Q = -(1+theta)/2*Ci@A - eps_i*np.eye(3)
                R = (1+theta)/2*Ci@B - eps_i*np.eye(3)
                S = -(1+theta)/2*Ci@B - eps_i*np.eye(3)
                T = (1 - 2*eps_i - 2*eps_i)*np.eye(3)
                assert sp.linalg.issymmetric(P, rtol=1e-6)
                assert sp.linalg.issymmetric(Q, rtol=1e-6)
                assert sp.linalg.issymmetric(R, rtol=1e-6)
                assert sp.linalg.issymmetric(S, rtol=1e-6)
                assert sp.linalg.issymmetric(T, rtol=1e-6)
                g = Ci@np.array([0, 0, V])
                if k > 0:
                    eps_e = eps_i / 2
                    N_xi = la.norm(Ci@A)
                    N_eta = la.norm(Ci@B)
                    d_eta = calc_d_eta(all_layers_y, all_layers_z, i, j, k)
                    d_xi = calc_d_xi(all_layers_y, all_layers_z, i, j, k)
                    Sl = np.sqrt(np.min([k, k_max//4*3])/k_max)
                    a_xi, a_eta = calc_axi_aeta(all_layers_y, all_layers_z, i, j, k)
                    R_xi = Sl * np.max([d_xi**(2/Sl), 0.1])* a_xi
                    R_eta = Sl * np.max([d_eta**(2/Sl), 0.1]) * a_eta
                    eps_e_xi = eps_e*R_xi*N_xi
                    eps_e_eta = eps_e*R_eta*N_eta
                    r_xi_xi, r_eta_eta = calc_rxixi_retaeta(all_layers_y, all_layers_z, i, j, k)
                    # g += -eps_e_xi*r_xi_xi - eps_e_eta*r_eta_eta

                big_mat[j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3] = T
                if j > 0:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        (j-1)*x_count*3 + i*3 : (j-1)*x_count*3 + i*3 + 3
                    ] = S
                else:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        (n_segments-1)*x_count*3 + i*3 : (n_segments-1)*x_count*3 + i*3 + 3
                    ] = S

                if j < n_segments-1:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        (j+1)*x_count*3 + i*3 : (j+1)*x_count*3 + i*3 + 3
                    ] = R
                else:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        0*x_count*3 + i*3 : 0*x_count*3 + i*3 + 3
                    ] = R
                if i > 0:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        j*x_count*3 + (i-1)*3 : j*x_count*3 + (i-1)*3 + 3
                    ] = Q
                else:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        j*x_count*3 + (x_count-1)*3 : j*x_count*3 + (x_count-1)*3 + 3
                    ] = Q
                if i < x_count-1:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        j*x_count*3 + (i+1)*3 : j*x_count*3 + (i+1)*3 + 3
                    ] = P
                else:
                    big_mat[
                        j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3,
                        j*x_count*3 + 0*3 : j*x_count*3 + 0*3 + 3
                    ] = P
                big_g[j*x_count*3 + i*3 : j*x_count*3 + i*3 + 3] = g
        # if k == 1:
        #     print(big_mat)
        #     print(list(big_mat[0]))
        #     print(big_g)
        y_fine_next = np.zeros_like(y_fine)
        z_fine_next = np.zeros_like(z_fine)
        delta = la.solve(big_mat, big_g)
        for j in range(n_segments):
            delta_seg = delta[j*x_count*3  + (x_count//2)*3 : j*x_count*3 + (x_count//2)*3 + 3]
            delta_seg = delta[j*x_count*3  : j*x_count*3  + 3]
            y_fine_next[j] = y_fine[j] + delta_seg[1]
            z_fine_next[j] = z_fine[j] + delta_seg[2]
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



