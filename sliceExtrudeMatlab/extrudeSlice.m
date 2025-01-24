function stacked_array = extrudeSlice(pp, arc_length, arc_length_at_max_y)
    n_segments = 80;
    k_max = 10;
    layer_size = 20;
    R_max = 1000;
    
    s_fine = [linspace(0, arc_length_at_max_y, n_segments / 2 + 1)(1:end-1), ...
          linspace(arc_length_at_max_y, arc_length, n_segments / 2 + 1)(1:end-1)];
    points = ppval(pp, s_fine)';

    % plotSlice(points)
    point_at_max = ppval(pp, arc_length_at_max_y);
    % ppval(pp, 0)
    z_fine = points(:, 3);
    y_fine = points(:, 2);

    angle_offset = atan2(points(1, 3), -points(1, 2));
    angle_offset_max = atan2(point_at_max(3), -point_at_max(2))
    if angle_offset_max < 0
        angle_offset_max = angle_offset_max + 2 * pi;
    end

    % Compute farfield points
    farfield_angles = [linspace(angle_offset, angle_offset_max, n_segments / 2 + 1)(1:end-1), ...
          linspace(-(pi*2-angle_offset_max), angle_offset , n_segments / 2 + 1)(1:end-1)];
    farfield_points = zeros(length(s_fine), 2);
    for t = 1:n_segments
        farfield_points(t, :) = -R_max * [cos(farfield_angles(t)),...
                                          -sin(farfield_angles(t))];
    end

    % Initialize layers
    all_layers_y = cell(k_max, 1);
    all_layers_z = cell(k_max, 1);
    all_layers_y{1} = points(:,2);
    all_layers_z{1} = points(:,3);

    left_bisector = calculateBisect([all_layers_y{1}, all_layers_z{1}], 1);
    right_bisector = calculateBisect([all_layers_y{1}, all_layers_z{1}], n_segments/2 + 1);
    [left_trisect1, left_trisect2] = calculateTrisect([all_layers_y{1}, all_layers_z{1}], 1);
    [right_trisect1, right_trisect2] = calculateTrisect([all_layers_y{1}, all_layers_z{1}], n_segments/2 + 1);

    % First layer with diamonds
    for k = 1:1
        y_fine_next = [];
        z_fine_next = [];
        for t = 1:n_segments
            area = calculateArea([all_layers_y{k}, all_layers_z{k}], t);
            area = 1;
            if t == 1
                bisector = left_trisect1;
                U = [all_layers_y{k}(t); all_layers_z{k}(t)] + bisector*layer_size/area;
                y_fine_next(end+1) = U(1);
                z_fine_next(end+1) = U(2);
            end
            if t == n_segments/2 + 1
                bisector = right_trisect1;
                U = [all_layers_y{k}(t); all_layers_z{k}(t)] + bisector*layer_size/area;
                y_fine_next(end+1) = U(1);
                z_fine_next(end+1) = U(2);
            end
            bisector = calculateBisect([all_layers_y{k}, all_layers_z{k}], t);
            U = [all_layers_y{k}(t); all_layers_z{k}(t)] + bisector*layer_size/area;
            y_fine_next(end+1) = U(1);
            z_fine_next(end+1) = U(2);

            if t == 1
                bisector = left_trisect2;
                U = [all_layers_y{k}(t); all_layers_z{k}(t)] + bisector*layer_size/area;
                y_fine_next(end+1) = U(1);
                z_fine_next(end+1) = U(2);
            end
            if t == n_segments/2 + 1
                bisector = right_trisect2;
                U = [all_layers_y{k}(t); all_layers_z{k}(t)] + bisector*layer_size/area;
                y_fine_next(end+1) = U(1);
                z_fine_next(end+1) = U(2);
            end
        end
        all_layers_y{k + 1} = y_fine_next';
        all_layers_z{k + 1} = z_fine_next';
    end

    % More layers
    for k = 2:k_max - 1
        y_fine_next = zeros(size(all_layers_y{k},1),1);
        z_fine_next = zeros(size(all_layers_y{k},1),1);
        for t = 1:size(all_layers_y{k}, 1)/2
            bisector = calculateBisect([all_layers_y{k}, all_layers_z{k}], t);
            area = calculateArea([all_layers_y{k}, all_layers_z{k}], t);
            area = 1;
            U = [all_layers_y{k}(t); all_layers_z{k}(t)] + bisector*layer_size/area;
            y_fine_next(t) = U(1);
            z_fine_next(t) = U(2);
        end
        for t = size(all_layers_y{k}, 1)/2 + 1:size(all_layers_y{k}, 1)
            bisector = calculateBisect([all_layers_y{k}, all_layers_z{k}], t);
            area = calculateArea([all_layers_y{k}, all_layers_z{k}], t);
            area = 1;
            U = [all_layers_y{k}(t); all_layers_z{k}(t)] + bisector*layer_size/area;
            y_fine_next(t) = U(1);
            z_fine_next(t) = U(2);
        end
        all_layers_y{k + 1} = y_fine_next;
        all_layers_z{k + 1} = z_fine_next;
    end

    % % Compute intermediate layers
    % for k = 1:k_max - 1
    %     s = (k + 1) / k_max;
    %     y_fine_next = zeros(size(y_fine));
    %     z_fine_next = zeros(size(z_fine));
    %     for t = 1:n_segments
    %         U = (1 - s) * [y_fine(t); z_fine(t)] + s * farfield_points(t, :)';
    %         y_fine_next(t) = U(1);
    %         z_fine_next(t) = U(2);
    %     end
    %     all_layers_y{k + 1} = y_fine_next;
    %     all_layers_z{k + 1} = z_fine_next;
    % end
    %
    %

    % Assemble the final extruded array
    stacked_array = zeros(k_max, n_segments + 4, 3);
    for k = 2:k_max
        stacked_array(k-1, :, 2) = all_layers_y{k};
        stacked_array(k-1, :, 3) = all_layers_z{k};
    end
    plotMesh([all_layers_y{1}, all_layers_z{1}], stacked_array);
end

function [trivec1, trivec2]=calculateTrisect(points, i)
    N = size(points, 1);
    prev_point = points(mod(i-2, N) + 1, :);
    curr_point = points(i, :);
    next_point = points(mod(i, N) + 1, :);
    
    vec1 = next_point - curr_point;
    vec2 = prev_point - curr_point;
    
    vec1 = vec1 / norm(vec1);
    vec2 = vec2 / norm(vec2);
    
    angle = atan2d(vec2(2), vec2(1)) - atan2d(vec1(2), vec1(1));
    if angle < 0
        angle = angle + 360; 
    end
    
    angle2 = angle / 3;
    angle1 = 2 * angle / 3;
    
    trivec1 = rotateVector(vec1, angle1);
    trivec2 = rotateVector(vec1, angle2);
end

function rotated_vec = rotateVector(vec, angle)
    rotation_matrix = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
    rotated_vec = rotation_matrix * vec(:);
end

function area = calculateArea(points,i)
    N = size(points, 1);

    prev_point = points(mod(i-2, N) + 1, :);
    curr_point = points(i, :);
    next_point = points(mod(i, N) + 1, :);
    
    vec1 = prev_point - curr_point;
    vec2 = next_point - curr_point;
    area = norm(vec1)/2 + norm(vec2)/2;
end

function bisector = calculateBisect(points,i)
    N = size(points, 1);

    prev_point = points(mod(i-2, N) + 1, :);
    curr_point = points(i, :);
    next_point = points(mod(i, N) + 1, :);
    
    vec1 = prev_point - curr_point;
    vec2 = next_point - curr_point;
    vec1 = vec1 / norm(vec1);
    vec2 = vec2 / norm(vec2);
    
    bisector = vec1 + vec2;
    bisector = bisector' / norm(bisector);
    
    cross_prod = cross([vec1, 0], [vec2, 0]);
    
    if cross_prod(3) > 0
        bisector = -bisector;
    end
end


function plotMesh(first_layer, stacked_array)
    figure;
    hold on;

    layer_y = first_layer(:, 1);
    layer_z = first_layer(:, 2);

    layer_y_closed = [layer_y; layer_y(1)];
    layer_z_closed = [layer_z; layer_z(1)];

    plot(layer_y_closed, layer_z_closed, 'k-', 'LineWidth', 0.5);

    next_layer_y = stacked_array(1, :, 2)';
    next_layer_z = stacked_array(1, :, 3)';

    for j = 1:size(first_layer,1)/2
        k = j + 2;
        if j == 1
            plot(
                [layer_y(j), next_layer_y(j)], ...
                [layer_z(j), next_layer_z(j)], ...
                'k-', 'LineWidth', 0.5
            );
        end
        plot(
            [layer_y(j), next_layer_y(k)], ...
            [layer_z(j), next_layer_z(k)], ...
            'k-', 'LineWidth', 0.5
        );
    end
    for j = size(first_layer,1)/2+1:size(first_layer,1)
        k = j + 4;
        if j == size(first_layer,1)/2+1
            plot(
                [layer_y(j), next_layer_y(j+2)], ...
                [layer_z(j), next_layer_z(j+2)], ...
                'k-', 'LineWidth', 0.5
            );
        end
        plot(
            [layer_y(j), next_layer_y(k)], ...
            [layer_z(j), next_layer_z(k)], ...
            'k-', 'LineWidth', 0.5
        );
    end

    for k = 1:size(stacked_array, 1)
        layer_y = stacked_array(k, :, 2)';
        layer_z = stacked_array(k, :, 3)';

        layer_y_closed = [layer_y; layer_y(1)];
        layer_z_closed = [layer_z; layer_z(1)];

        plot(layer_y_closed, layer_z_closed, 'k-', 'LineWidth', 0.5);
    end

    % Plot the connections between consecutive layers
    for i = 1:(size(stacked_array, 1)-2)
        layer_y = stacked_array(i, :, 2)';
        layer_z = stacked_array(i, :, 3)';
        next_layer_y = stacked_array(i+1, :, 2)';
        next_layer_z = stacked_array(i+1, :, 3)';

        for j = 1:length(layer_y)
            plot(
                [layer_y(j), next_layer_y(j)], ...
                [layer_z(j), next_layer_z(j)], ...
                'k-', 'LineWidth', 0.5
            );
        end
    end

    title('Multiple Layer Extrusion');
    xlabel('Y');
    ylabel('Z');
    axis equal;
    hold off;
end


