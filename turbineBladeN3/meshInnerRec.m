function [elements,boundaries] = meshInnerRec(pp, arc_length)
    config
    tail_count = k_inner+4;




    N_total = n_top + 2*k_inner+1;
    i = 0:(N_total-1);
    cheb_nodes = cos(i*pi / (N_total-1));
    cheb_nodes_mapped = (cheb_nodes + 1) / 2;
    s_fine = [linspace(0, arc_length, n_top+n_bottom+4*k_inner + 1)(1:end-1)];


    all_points = ppval(pp, s_fine)';
    points_top = ppval(pp, s_fine(1:n_top + 2*k_inner+1))';
    points_bottom = [ppval(pp, s_fine(n_top + 2*k_inner + 1:end))'; points_top(1,:)];
    n = n_top + 2*k_inner;


    point_at_min = points_top(1,:)';
    point_at_max = points_bottom(1,:)';
    z = point_at_min(3);

    elements = [];
    boundaries = [];

    % Left tip
    p14 = findBisectNode(point_at_min', points_top(2, :), points_bottom(end-1, :));
    element = [];
    element(1,:) = point_at_min';
    element(2,:) = points_bottom(end-1, :);
    element(3,:) = p14;
    element(4,:) = points_top(2, :);
    elements(end+1, :, :) = element;

    all_points_sub = [];
    for i=tail_count:2*N_total+1-tail_count
        p1 = all_points(i, :);
        p2 = all_points(i-1, :);
        p3 = all_points(i+1, :);
        p4 = p1 - 0.04*findBisect(p1, p2, p3);
        p5 = p2 - 0.04*findBisect(p2, all_points(i-2, :), p1);
        all_points_sub(end+1, :) = p5;
        if i==tail_count
            k1 = p5;
        end
        if i==N_total
            k4 = p4;
        end
        if i==2*N_total+1-tail_count
            k2 = p4;
        end

        element = [];
        element(1,:) = p1;
        element(2,:) = p2;
        element(3,:) = p5;
        element(4,:) = p4;
        elements(end+1, :, :) = element;
    end
    all_points_sub(end+1, :) = p4;

    k3 = (k1 + k2)/2;
    sub_count = size(all_points_sub, 1);
    for i=1:(sub_count-3)/2
        rat = (i+1)/((sub_count-1)/2);
        ratp = (i)/((sub_count-1)/2);

        p1 = all_points_sub(i+1, :);
        p2 = all_points_sub(i, :);
        p3 = ratp*k4 + (1-ratp)*k3;
        p4 = rat*k4 + (1-rat)*k3;
        if i==1
            k6 = p3;
        end

        element = [];
        element(1,:) = p1;
        element(2,:) = p2;
        element(3,:) = p3;
        element(4,:) = p4;
        elements(end+1, :, :) = element;

        p1 = rat*k4 + (1-rat)*k3;
        p2 = ratp*k4 + (1-ratp)*k3;
        p3 = all_points_sub(end-i+1, :);
        p4 = all_points_sub(end-i, :);

        element = [];
        element(1,:) = p1;
        element(2,:) = p2;
        element(3,:) = p3;
        element(4,:) = p4;
        elements(end+1, :, :) = element;
    end

    for i=3:tail_count-2
        rat = (i-2)/(tail_count-3);
        ratp = (i-3)/(tail_count-3);
        p1 = all_points(i, :);
        p2 = all_points(i-1, :);
        p3 = ratp*k3 + (1-ratp)*p14;
        p4 = rat*k3 + (1-rat)*p14;

        element = [];
        element(1,:) = p1;
        element(2,:) = p2;
        element(3,:) = p3;
        element(4,:) = p4;
        elements(end+1, :, :) = element;

        p1 = all_points(end-i+3, :);
        p2 = all_points(end-i+2, :);
        p3 = rat*k3 + (1-rat)*p14;
        p4 = ratp*k3 + (1-ratp)*p14;

        element = [];
        element(1,:) = p1;
        element(2,:) = p2;
        element(3,:) = p3;
        element(4,:) = p4;
        elements(end+1, :, :) = element;
    end
    element = [];
    element(1,:) = k1;
    element(2,:) = points_top(tail_count-1, :);
    element(3,:) = points_top(tail_count-2, :);
    element(4,:) = p3;
    elements(end+1, :, :) = element;
    element = [];
    element(1,:) = p3;
    element(2,:) = points_bottom(end-tail_count+3, :);
    element(3,:) = points_bottom(end-tail_count+2, :);
    element(4,:) = k2;
    elements(end+1, :, :) = element;
    element = [];
    element(1,:) = k1;
    element(2,:) = p3;
    element(3,:) = k2;
    element(4,:) = k6;
    elements(end+1, :, :) = element;


    checkCounterClockwise(elements)
    boundaries_coords = ppval(pp, s_fine)';
    % size(boundaries_coords)
    elements = relaxQuadMesh(elements, boundaries_coords, 500);
end

function p4 = findBisectNode(p1, p2, p3)
    % Vectors defining the plane
    p4 = p1 + 2.25/4*(p2 - 2*p1 + p3);
end
% function p4 = findBisectNode(p1, p2, p3)
%     % Vectors defining the plane
%     P1 = p1(2:3);
%     P2 = p2(2:3);
%     P3 = p3(2:3);
%     v12 = P2 - P1; % Vector from P1 to P2
%     v13 = P3 - P1; % Vector from P1 to P3
%
%     A = [v12(1), v12(2);
%          v13(1), v13(2)];
%     b = [dot(v12, P2);
%          dot(v13, P3)];
%     P4 = (A \ b)';
%     p4= [];
%     p4(1) = p1(1);
%     p4(2:3) = P4;
% end

function U=bilinearInterp(corners, s, t, sp, tp)
    P0_t = (1-t) * corners(1, :) + t * corners(2, :);
    P1_t = (1-tp) * corners(4, :) + tp * corners(3, :);
    P2_s = (1-s) * corners(1, :) + s * corners(4, :);
    P3_s = (1-sp) * corners(2, :) + sp * corners(3, :);
    P00 = corners(1, :);
    P10 = corners(4, :);
    P01 = corners(2, :);
    P11 = corners(3, :);
    
    U = (1-s) * P0_t + s * P1_t + (1-t) * P2_s + t * P3_s ...
        - (1-s)*(1-t)*P00 - s*(1-t)*P10 - (1-s)*t*P01 - s*t*P11;

end


function checkCounterClockwise(elements)
    for k = 1:size(elements, 1)
        element = squeeze(elements(k, :, :));
        for i = 1:4
            v1 = element(i, :);
            v2 = element(mod(i, 4) + 1, :);
            v3 = element(mod(i + 1, 4) + 1, :);
            edge1 = v2 - v1;
            edge2 = v3 - v2;
            
            cross_prod = cross(edge1, edge2);
            
            if cross_prod(3) <= 0
                fprintf('Inner Element %d is not counterclockwise at corner %d.\n', k, i);
            end
        end
    end
end


function bisector = findBisect(p1, p2, p3)
    vec1 = p2 - p1;
    vec2 = p3 - p1;
    vec1 = vec1 / norm(vec1);
    vec2 = vec2/ norm(vec2);
    
    bisector = vec1 + vec2;
    bisector = bisector' / norm(bisector);
    
    cross_prod = cross([vec1], [vec2]);
    
    if cross_prod(3) > 0
        bisector = -bisector;
    end
    bisector = bisector';
end
