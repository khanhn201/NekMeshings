function relaxed_elements = relaxQuadMesh(elements, boundary_coords, max_iter)
    if nargin < 3
        max_iter = 20; % number of smoothing iterations
    end

    N = size(elements, 1);

    % Step 1: Flatten and get unique points
    all_pts = reshape(elements, [], 3);  % (4*N, 3)
    [unique_pts, ~, idx_map] = uniquetol(all_pts,'ByRows',true);  % (P, 3)
    P = size(unique_pts, 1);

    % Step 2: Build element-to-node mapping
    elem_conn = reshape(idx_map, [N, 4]);  % (N, 4)

    % Step 3: Identify boundary node indices with tolerance
    tol = 1e-8;
    is_fixed = false(P,1);
    for i = 1:size(boundary_coords,1)
        dists = sqrt(sum((unique_pts - boundary_coords(i,:)).^2, 2));
        idx = find(dists < tol, 1);
        if ~isempty(idx)
            is_fixed(idx) = true;
        else
            warning('Boundary node not found for point %d', i);
        end
    end

    % Step 4: Build adjacency list
    adj = cell(P,1);
    for e = 1:N
        nodes = elem_conn(e,:);
        for i = 1:4
            ni = nodes(i);
            nj = nodes(mod(i,4)+1);  % cyclic
            adj{ni}(end+1) = nj;
            adj{nj}(end+1) = ni;
        end
    end

    % Step 5: Laplacian smoothing
    relaxed_pts = unique_pts;

    for iter = 1:max_iter
        new_pts = relaxed_pts;
        for i = 1:P
            if is_fixed(i), continue; end  % skip boundary nodes
            neighbors = unique(adj{i});
            new_pts(i,:) = mean(relaxed_pts(neighbors,:), 1);
        end
        relaxed_pts = new_pts;
    end

    % Step 6: Map smoothed points back to elements
    relaxed_elements = elements;
    for k = 1:numel(idx_map)
        relaxed_elements(k) = relaxed_pts(idx_map(k));
    end
    relaxed_elements = reshape(relaxed_elements, size(elements));
end
