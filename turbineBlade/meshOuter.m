function [elements,boundaries, pp_coarse] = meshOuter(pp, arc_length, arc_length_at_max_y)
    config

    s_fine = [linspace(0, arc_length_at_max_y, n_top+1)(1:end-1), ...
          linspace(arc_length_at_max_y, arc_length, n_bottom+1)(1:end-1)];
    points_top = ppval(pp, s_fine(1:n_top+1))';
    % points_leading = ppval(pp, s_fine(n_top + 1:n_top + 2*n_leading + 1))';
    points_bottom = [ppval(pp, s_fine(n_top + 1:end))'; points_top(1,:)];

    point_at_min = ppval(pp, 0);
    point_at_max = ppval(pp, arc_length_at_max_y);
    x = point_at_min(1);

    elements = [];
    boundaries = [];

    all_points = [ppval(pp, s_fine)'; points_top(1,:)];
    pp_coarse = splinefit([s_fine'; arc_length], all_points', [s_fine'; arc_length], "periodic", true);
    %

    % Quad top
    for k = 1:k_inner+1
        layer_next = [];
        s = (mult^k-mult)/(mult^(k_inner+1)-mult)*delta_inner;
        for t = 1:n_top+1
            r = (t - 1)/(n_top);
            U = (1 - s) * points_top(t,:) + s * [x; -R_a + r*2*R_a; R_t]';
            layer_next = [layer_next; U];
            if (k > 1) && (t > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(t-1, :);
                element(4,:) = layer_prev(t, :);
                elements(end+1, :, :) = element;
                if k == 2
                    boundaries(end+1, :) = [size(elements, 1); 1;];
                end
            end
        end        
        layer_prev = layer_next;
    end
    for k = 1:k_outer+1
        layer_next = [];
        s = delta_inner + (1-delta_inner)*(k-1)/(k_outer);
        for t = 1:n_top+1
            r = (t - 1)/(n_top);
            U = (1 - s) * points_top(t,:) + s * [x; -R_a + r*2*R_a; R_t]';
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
    for k = 1:k_inner+1
        layer_next = [];
        s = (mult^k-mult)/(mult^(k_inner+1)-mult)*delta_inner;
        for t = 1:n_bottom + 1
            r = (t - 1)/(n_bottom);
            U = (1 - s) * points_bottom(t,:) + s * [x; R_a - r*2*R_a; R_b]';
            layer_next = [layer_next; U];
            if (k > 1) && (t > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(t-1, :);
                element(4,:) = layer_prev(t, :);
                elements(end+1, :, :) = element;
                if k == 2
                    boundaries(end+1, :) = [size(elements, 1); 1;];
                end
            end
        end
        layer_prev = layer_next;
    end
    for k = 1:k_outer+1
        layer_next = [];
        s = delta_inner + (1-delta_inner)*(k-1)/(k_outer);
        for t = 1:n_bottom + 1
            r = (t - 1)/(n_bottom);
            U = (1 - s) * points_bottom(t,:) + s * [x; R_a - r*2*R_a; R_b]';
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
    % Quad right
    % for k = 1:k_inner+1
    %     layer_next = [];
    %     s = (mult^k-mult)/(mult^(k_inner+1)-mult)*delta_inner;
    %     for t = 1:2*n_leading + 1
    %         r = (t - 1)/(2*n_leading);
    %         U = (1 - s) * points_leading(t,:) + s * [x; R_a; R_t - r*(R_t-R_b)]';
    %         layer_next = [layer_next; U];
    %         if (k > 1) && (t > 1)
    %             element = [];
    %             element(1,:) = layer_next(end, :);
    %             element(2,:) = layer_next(end-1, :);
    %             element(3,:) = layer_prev(t-1, :);
    %             element(4,:) = layer_prev(t, :);
    %             elements(end+1, :, :) = element;
    %             if k == 2
    %                 boundaries(end+1, :) = [size(elements, 1); 1;];
    %             end
    %         end
    %     end
    %     layer_prev = layer_next;
    % end
    % for k = 1:k_outer+1
    %     layer_next = [];
    %     s = delta_inner + (1-delta_inner)*(k-1)/(k_outer);
    %     for t = 1:2*n_leading + 1
    %         r = (t - 1)/(2*n_leading);
    %         U = (1 - s) * points_leading(t,:) + s * [x; R_a; R_t - r*(R_t-R_b)]';
    %         layer_next = [layer_next; U];
    %         if (k > 1) && (t > 1)
    %             element = [];
    %             element(1,:) = layer_next(end, :);
    %             element(2,:) = layer_next(end-1, :);
    %             element(3,:) = layer_prev(t-1, :);
    %             element(4,:) = layer_prev(t, :);
    %             elements(end+1, :, :) = element;
    %             if k == k_outer+1
    %                 boundaries(end+1, :) = [size(elements, 1); 4;];
    %             end
    %         end
    %     end
    %     layer_prev = layer_next;
    % end

    right_diamond_points = [
        point_at_max';
        (1-delta_inner)*point_at_max' + delta_inner*[x; R_a; R_b]';
        (1-diamond_mult*delta_inner)*point_at_max' + diamond_mult*delta_inner*[x; R_a; (R_t + R_b)/2]';
        (1-delta_inner)*point_at_max' + delta_inner*[x; R_a; R_t]';
    ];
 
    for k = 1:k_inner+1
        layer_next = [];

        s = (mult^k-mult)/(mult^(k_inner+1)-mult);
        sp = (k-1)/k_inner;
        for tt = 1: k_inner+1
            t = (mult^tt-mult)/(mult^(k_inner+1)-mult);
            tp = (tt- 1)/k_inner;
            U = bilinearInterp(right_diamond_points, s, t, sp, tp);
            layer_next = [layer_next; U];
            if (k > 1) && (tt > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(tt-1, :);
                element(4,:) = layer_prev(tt, :);
                elements(end+1, :, :) = element;
            end
        end
        layer_prev = layer_next;
    end

    % Quad right outer
    right_diamond_top_points = [
        (1-delta_inner)*point_at_max' + delta_inner*[x; R_a; R_t]';
        (1-diamond_mult*delta_inner)*point_at_max' + diamond_mult*delta_inner*[x; R_a; (R_t + R_b)/2]';
        [x; R_a; (R_t + R_b)/2]';
        [x; R_a; R_t]';
    ];
    for k = 1:k_outer+1
        layer_next = [];
        sp = (k-1)/k_outer;
        for tt = 1: k_inner+1
            tp = (tt- 1)/k_inner;
            U = bilinearInterp(right_diamond_top_points, sp, tp, sp, tp);
            layer_next = [layer_next; U];
            if (k > 1) && (tt > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(tt-1, :);
                element(4,:) = layer_prev(tt, :);
                elements(end+1, :, :) = element;
                if k == k_outer+1
                    boundaries(end+1, :) = [size(elements, 1); 5;];
                end
            end
        end
        layer_prev = layer_next;
    end
    right_diamond_bottom_points = [
        (1-diamond_mult*delta_inner)*point_at_max' + diamond_mult*delta_inner*[x; R_a; (R_t + R_b)/2]';
        (1-delta_inner)*point_at_max' + delta_inner*[x; R_a; R_b]';
        [x; R_a; R_b]';
        [x; R_a; (R_t + R_b)/2]';
    ];
    for k = 1:k_outer+1
        layer_next = [];
        sp = (k-1)/k_outer;
        for tt = 1: k_inner+1
            tp = (tt- 1)/k_inner;
            U = bilinearInterp(right_diamond_bottom_points, sp, tp, sp, tp);
            layer_next = [layer_next; U];
            if (k > 1) && (tt > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(tt-1, :);
                element(4,:) = layer_prev(tt, :);
                elements(end+1, :, :) = element;
                if k == k_outer+1
                    boundaries(end+1, :) = [size(elements, 1); 6;];
                end
            end
        end
        layer_prev = layer_next;
    end

    % Quad left
    left_diamond_points = [
        point_at_min';
        (1-delta_inner)*point_at_min' + delta_inner*[x; -R_a; R_t]';
        (1-diamond_mult*delta_inner)*point_at_min' + diamond_mult*delta_inner*[x; -R_a; (R_t + R_b)/2]';
        (1-delta_inner)*point_at_min' + delta_inner*[x; -R_a; R_b]';
    ];
 
    for k = 1:k_inner+1
        layer_next = [];

        s = (mult^k-mult)/(mult^(k_inner+1)-mult);
        sp = (k-1)/k_inner;
        for tt = 1: k_inner+1
            t = (mult^tt-mult)/(mult^(k_inner+1)-mult);
            tp = (tt- 1)/k_inner;
            U = bilinearInterp(left_diamond_points, s, t, sp, tp);
            layer_next = [layer_next; U];
            if (k > 1) && (tt > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(tt-1, :);
                element(4,:) = layer_prev(tt, :);
                elements(end+1, :, :) = element;
            end
        end
        layer_prev = layer_next;
    end

    % Quad left outer
    left_diamond_top_points = [
        (1-diamond_mult*delta_inner)*point_at_min' + diamond_mult*delta_inner*[x; -R_a; (R_t + R_b)/2]';
        (1-delta_inner)*point_at_min' + delta_inner*[x; -R_a; R_t]';
        [x; -R_a; R_t]';
        [x; -R_a; (R_t + R_b)/2]';
    ];
    for k = 1:k_outer+1
        layer_next = [];
        sp = (k-1)/k_outer;
        for tt = 1: k_inner+1
            tp = (tt- 1)/k_inner;
            U = bilinearInterp(left_diamond_top_points, sp, tp, sp, tp);
            layer_next = [layer_next; U];
            if (k > 1) && (tt > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(tt-1, :);
                element(4,:) = layer_prev(tt, :);
                elements(end+1, :, :) = element;
                if k == k_outer+1
                    boundaries(end+1, :) = [size(elements, 1); 7;];
                end
            end
        end
        layer_prev = layer_next;
    end
    left_diamond_bottom_points = [
        (1-delta_inner)*point_at_min' + delta_inner*[x; -R_a; R_b]';
        (1-diamond_mult*delta_inner)*point_at_min' + diamond_mult*delta_inner*[x; -R_a; (R_t + R_b)/2]';
        [x; -R_a; (R_t + R_b)/2]';
        [x; -R_a; R_b]';
    ];
    for k = 1:k_outer+1
        layer_next = [];
        sp = (k-1)/k_outer;
        for tt = 1: k_inner+1
            tp = (tt- 1)/k_inner;
            U = bilinearInterp(left_diamond_bottom_points, sp, tp, sp, tp);
            layer_next = [layer_next; U];
            if (k > 1) && (tt > 1)
                element = [];
                element(1,:) = layer_next(end, :);
                element(2,:) = layer_next(end-1, :);
                element(3,:) = layer_prev(tt-1, :);
                element(4,:) = layer_prev(tt, :);
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


