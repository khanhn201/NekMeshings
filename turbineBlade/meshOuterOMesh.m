function [elements,boundaries, pp_coarse] = meshOuterOMesh(pp, arc_length, arc_length_at_max_y, flipped)
    config


    s_fine = [linspace(0, arc_length_at_max_y, n_top + 2*k_inner +1)(1:end-1), ...
          linspace(arc_length_at_max_y, arc_length, n_bottom + 2*k_inner + 1)(1:end-1)];

    if flipped == true 
        s_fine(s_fine >= arc_length_at_max_y) = s_fine(s_fine >= arc_length_at_max_y) - arc_length;
    end
    points_top = ppval(pp, s_fine(1:n_top + 2*k_inner+1))';
    % points_leading = ppval(pp, s_fine(n_top + 1:n_top + 2*n_leading + 1))';
    points_bottom = [ppval(pp, s_fine(n_top + 2*k_inner + 1:end))'; points_top(1,:)];

    point_at_min = ppval(pp, 0);
    point_at_max = ppval(pp, arc_length_at_max_y);
    angle_offset = atan2(point_at_min(3), -point_at_min(2)) + atan2(-point_at_max(3), point_at_max(2));
    angle_offset = angle_offset/2;
    x = point_at_min(1);

    elements = [];
    boundaries = [];

    % pp_coarse = splinefit([s_fine'; arc_length], all_points', [s_fine'; arc_length], "periodic", true);
    all_points = ppval(pp, s_fine)';
    if flipped == true
        s_fine1 = circshift(s_fine, -n_top - 2*k_inner);
        all_points1 = circshift(all_points, -n_top - 2*k_inner);
        s_fine1 = [s_fine1'; arc_length_at_max_y];
        all_points1 = [all_points1; all_points1(1, :)];
    else
        s_fine1 = [s_fine'; arc_length];
        all_points1 = [all_points; points_top(1,:)];
    end
    pp_coarse = spline(s_fine1, all_points1');

    layer_prev = all_points;
    for k = 1:2
        layer_next = [];
        for t = 1:length(all_points)
            p1 = all_points(t, :);
            if t == 1
                p2 = all_points(end, :);
            else
                p2 = all_points(t - 1, :);
            end
            if t == length(all_points)
                p3 = all_points(1, :);
            else
                p3 = all_points(t + 1, :);
            end
            bisector = findBisect(p1, p2, p3);
            U = p1 + first_layer_thickness*(mult^k)*bisector;
            layer_next = [layer_next; U];
            if (t > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(t-1, :);
                element(4,:) = layer_prev(t, :);
                elements(end+1, :, :) = element;
                if k == 1
                    boundaries(end+1, :) = [size(elements, 1); 1;];
                end
            end
        end
        element = [];
        element(1,:) = layer_next(1, :);
        element(2,:) = layer_next(end, :);
        element(3,:) = layer_prev(end, :);
        element(4,:) = layer_prev(1, :);
        elements(end+1, :, :) = element;
        if k == 1
            boundaries(end+1, :) = [size(elements, 1); 1;];
        end

        layer_prev = layer_next;
        all_points = layer_next;
    end

    % Quad top
    for k = 1:k_outer+1
        layer_next = [];
        s = (mult^k-mult)/(mult^(k_outer+1)-mult);
        Rz = [1, 0, 0;
              0, cos(angle_offset*(1-s)), sin(angle_offset*(1-s));
              0, -sin(angle_offset*(1-s)),  cos(angle_offset*(1-s));
              ];
        for t = 1:n_top+1
            r = (t - 1)/(n_top);
            U = (1 - s) * all_points(t + k_inner,:) + s * (Rz*[x; -R_a + r*2*R_a; R_t])';
            layer_next = [layer_next; U];
            if (k > 1) && (t > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(t-1, :);
                element(4,:) = layer_prev(t, :);
                elements(end+1, :, :) = element;
                if k == k_outer+1
                    boundaries(end+1, :) = [size(elements, 1); 2;];
                end
            end
        end        
        layer_prev = layer_next;
    end
    % Quad bottom
    for k = 1:k_outer+1
        layer_next = [];
        s = (mult^k-mult)/(mult^(k_outer+1)-mult);
        Rz = [1, 0, 0;
              0, cos(angle_offset*(1-s)), sin(angle_offset*(1-s));
              0, -sin(angle_offset*(1-s)),  cos(angle_offset*(1-s));
              ];
        for t = 1:n_bottom + 1
            r = (t - 1)/(n_bottom);
            U = (1 - s) * all_points(t + n_top  + 3* k_inner,:) + s * (Rz *[x; R_a - r*2*R_a; R_b])';
            layer_next = [layer_next; U];
            if (k > 1) && (t > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(t-1, :);
                element(4,:) = layer_prev(t, :);
                elements(end+1, :, :) = element;
                if k == k_outer+1
                    boundaries(end+1, :) = [size(elements, 1); 3;];
                end
            end
        end
        layer_prev = layer_next;
    end
    % Quad left
    for k = 1:k_outer+1
        layer_next = [];
        s = (mult^k-mult)/(mult^(k_outer+1)-mult);
        Rz = [1, 0, 0;
              0, cos(angle_offset*(1-s)), sin(angle_offset*(1-s));
              0, -sin(angle_offset*(1-s)),  cos(angle_offset*(1-s));
              ];
        for t = 1:2*k_inner + 1
            if t > k_inner
                point = all_points(t-k_inner,:);
            else
                point = all_points(end - k_inner + t,:);
            end
            r = (t - 1)/(2*k_inner);
            U = (1 - s) * point + s * (Rz *[x; -R_a; R_b + r*(R_t-R_b)])';
            layer_next = [layer_next; U];
            if (k > 1) && (t > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(t-1, :);
                element(4,:) = layer_prev(t, :);
                elements(end+1, :, :) = element;
                if k == k_outer+1
                    boundaries(end+1, :) = [size(elements, 1); 8;];
                end
            end
        end
        layer_prev = layer_next;
    end
    % Quad right
    for k = 1:k_outer+1
        layer_next = [];
        s = (mult^k-mult)/(mult^(k_outer+1)-mult);
        Rz = [1, 0, 0;
              0, cos(angle_offset*(1-s)), sin(angle_offset*(1-s));
              0, -sin(angle_offset*(1-s)),  cos(angle_offset*(1-s));
              ];
        for t = 1:2*k_inner + 1
            point = all_points(n_top + k_inner + t,:);
            r = (t - 1)/(2*k_inner);
            U = (1 - s) * point + s * (Rz*[x; R_a; R_t - r*(R_t-R_b)])';
            layer_next = [layer_next; U];
            if (k > 1) && (t > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(t-1, :);
                element(4,:) = layer_prev(t, :);
                elements(end+1, :, :) = element;
                if k == k_outer+1
                    boundaries(end+1, :) = [size(elements, 1); 8;];
                end
            end
        end
        layer_prev = layer_next;
    end

    checkCounterClockwise(elements)
end

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
            
            if cross_prod(1) <= 0
                fprintf('Element %d is not counterclockwise at corner %d.\n', k, i);
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
    
    if cross_prod(1) > 0
        bisector = -bisector;
    end
    bisector = bisector';
end
