function u_elem = solveLaplace(elements, p_coord, value)
    % Input:
    %   elements: (N, 4, 3) quad elements, with coordinates
    %   p_coord: (P, 3) boundary point coordinates
    %   value:   (P, 1) Dirichlet values for boundary points
    % Output:
    %   u: (total_pts x 1) solution at each mesh point

    N = size(elements, 1);

    % Step 1: Build global point list from elements
    all_pts = reshape(elements, [], 3); % (4N x 3)
    [unique_pts, ~, idx_map] = uniquetol(all_pts,'ByRows',true);  % (P, 3)
    total_pts = size(unique_pts, 1);

    % Step 2: Build connectivity: elements as indices into unique_pts
    elem_conn = reshape(idx_map, [N, 4]);  % (N x 4)

    % Step 3: Map boundary coordinates to global indices
    P = size(p_coord, 1);
    boundary_idx = zeros(P,1);
    for i = 1:P
        pt = p_coord(i, :);
    
        tol =1e-5;
        dists = sqrt(sum((unique_pts - pt).^2, 2));
        
        idx = find(dists < tol, 1);
        if isempty(idx)
            error('Boundary point not found within tolerance.');
        end
        boundary_idx(i) = idx;
    end

    % Step 4: Assemble stiffness matrix using Q1 elements
    A = sparse(total_pts, total_pts);
    b = zeros(total_pts,1);

    % Gauss points and weights
    gp = [-1/sqrt(3), 1/sqrt(3)];
    w = [1, 1];

    for e = 1:N
        nodes = elem_conn(e,:);
        coords = unique_pts(nodes, 1:2); % Only use x, y

        Ke = zeros(4,4);
        for i = 1:2
            for j = 1:2
                xi = gp(i); eta = gp(j);
                dN_dxi  = 0.25 * [-1+eta,  1-eta,  1+eta, -1-eta];
                dN_deta = 0.25 * [-1+xi,  -1-xi,  1+xi,   1-xi];

                J = [dN_dxi; dN_deta] * coords;
                detJ = det(J);
                invJ = inv(J);
                dN_dxy = invJ' * [dN_dxi; dN_deta];

                Ke = Ke + (dN_dxy' * dN_dxy) * detJ * w(i)*w(j);
            end
        end

        A(nodes,nodes) = A(nodes,nodes) + Ke;
    end

    % Step 5: Apply Dirichlet conditions
    u = zeros(total_pts,1);
    u(boundary_idx) = value;
    fixed = boundary_idx;
    free = setdiff(1:total_pts, fixed);

    b = b - A * u;
    u(free) = A(free, free) \ b(free);
    
    u_elem = reshape(u(idx_map), [N, 4]);
end
