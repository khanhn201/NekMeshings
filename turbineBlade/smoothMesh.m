function [smoothedElements] = smoothMesh(elements, boundaries)
    sf = 0.99;
    boundary_elements = boundaries(:,1);
    config;
    tol = 1e-6;
    n = size(elements, 1);
    X = reshape(elements, n*4, 3);
    [Xnew,ia,ic] = uniquetol(X,tol,'ByRows',true);
    Hexes = reshape(ic,n,4);
    X = Xnew;
    for i = 1:n_smooth
        X_new = X;
        for j = 1:size(X,1)
            [row_indices, ~] = find(Hexes == j);
            element_indices = unique(row_indices);

            is_boundary = ~isempty(intersect(element_indices, boundary_elements));
            if ~is_boundary
                neighbors_all = Hexes(element_indices, :)(:);
                neighbors_indices = unique(neighbors_all);
                neighbors_indices(neighbors_indices == j) = [];

                new_positions = zeros(length(element_indices), size(X,2));
                for k = 1:length(element_indices) % Shrink each element
                    element = Hexes(element_indices(k), :);
                    element_coords = X(element, :);
                    centroid = mean(element_coords, 1);
                    shrunk_coords = sf * element_coords + (1-sf)* centroid;
                    idx = find(element == j, 1);
                    new_positions(k, :) = shrunk_coords(idx, :);
                end

                X_new(j, :) = mean(new_positions, 1);
            end
        end
        X = X_new;
    end
    Hexes = Hexes(:);
    smoothedElements = X(Hexes, :);
    smoothedElements = reshape(smoothedElements, n, 4, 3);
end

