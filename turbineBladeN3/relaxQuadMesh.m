function relaxed_elements = relaxQuadMesh(elements, boundary_coords, max_iter)
    if nargin < 3
        max_iter = 20; % number of smoothing iterations
    end

    N = size(elements, 1);

    % Step 1: Flatten and get unique points
    all_pts = reshape(elements, [], 3);  % (4*N, 3)
    [unique_pts, ~, idx_map] = uniquetol(all_pts, 1e-8,'ByRows',true);  % (P, 3)
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
        % for i = 1:4
        % for j = 1:3
        %     ni = nodes(i);
        %     nj = nodes(mod(i+j-1,4)+1);  % cyclic
        %     adj{ni}(end+1) = nj;
        %     adj{nj}(end+1) = ni;
        % end
        % end
        for i = 1:4
            ni = nodes(i);
            nj = nodes(mod(i,4)+1);   % next neighbor
            nk = nodes(mod(i-2,4)+1); % previous neighbor
            adj{ni}(end+1) = nj;
            adj{ni}(end+1) = nk;
        end

    end

    % Step 5: Laplacian smoothing
    relaxed_pts = unique_pts;

    for iter = 1:max_iter
        new_pts = relaxed_pts;
        for i = 1:P
            if is_fixed(i), continue; end  % skip boundary nodes
            neighbors = unique(adj{i});
            new_pts(i,:) = mean(new_pts(neighbors,:), 1);
        end
        % relaxed_pts = 0.5*relaxed_pts + 0.5*new_pts;
        relaxed_pts = new_pts;
        % Plot mesh at current iteration
        % figure(1); clf; hold on; axis equal;
        % for e = 1:N
        %     pts = relaxed_pts(elem_conn(e,:), :);  % 4x3
        %     fill(pts([1:4 1],1), pts([1:4 1],2), [0.8 0.9 1], 'EdgeColor', 'k');
        % end
        % title(['Relaxation Iteration ', num2str(iter)]);
        % xlabel('X'); ylabel('Y');
        % drawnow;
    end

    % Step 6: Map smoothed points back to elements
    
    relaxed_elements = reshape(relaxed_pts(idx_map, :), size(elements));
end
