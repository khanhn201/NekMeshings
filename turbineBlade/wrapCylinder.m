function [elements, boundaries] = wrapCylinder(xs)
    config
    elements = [];
    boundaries = [];

    r_square = xs(2);
    r_cir = xs(end);
    function point=bottom_curve(u)
        if u <= k_cyl
            point = [
                -r_square;
                u/k_cyl*r_square;
                0;
            ]';
        end
        if u > k_cyl && u <= 3*k_cyl
            point = [
                (u-2*k_cyl)/k_cyl*r_square;
                r_square;
                0;
            ]';
        end
        if u > 3*k_cyl
            point = [
                r_square;
                (4*k_cyl-u)/k_cyl*r_square;
                0;
            ]';
        end
    end
    function point=top_curve(u)
        point = [
            -r_cir*cos(u/4/k_cyl*pi);
            r_cir*sin(u/4/k_cyl*pi);
            0;
        ]';
    end
    function point=curve_interp(u, v)
        point = v*top_curve(u) + (1-v)*bottom_curve(u);
    end

    p00 = [0;0;0]';
    p11 = [-r_square;r_square;0]';
    p12 = [r_square;r_square;0]';

    element = zeros(4,3);
    element(1, 1:3) = p11;
    element(2, 1:3) = [
        -r_square;
        (k_cyl-1)/k_cyl*r_square;
        0;
    ]';
    element(3, 1:3) = (k_cyl-1)/k_cyl*p11 + 1/k_cyl*p00;
    element(4, 1:3) = [
        -(k_cyl-1)/k_cyl*r_square;
        r_square;
        0;
    ]';
    elements(end+1, :, :) = element;

    element = zeros(4,3);
    element(1, 1:3) = p12;
    element(2, 1:3) = [
        (k_cyl-1)/k_cyl*r_square;
        r_square;
        0;
    ]';
    element(3, 1:3) = (k_cyl-1)/k_cyl*p12 + 1/k_cyl*p00;
    element(4, 1:3) = [
        r_square;
        (k_cyl-1)/k_cyl*r_square;
        0;
    ]';
    elements(end+1, :, :) = element;

    for i = 1:k_cyl-1
        element = zeros(4,3);
        element(1, 1:3) = [
            -r_square;
            i/k_cyl*r_square;
            0;
        ]';
        element(2, 1:3) = [
            -r_square;
            (i-1)/k_cyl*r_square;
            0;
        ]';
        element(3, 1:3) = (i-1)/k_cyl*p11 + (k_cyl-i+1)/k_cyl*p00;
        element(4, 1:3) = (i)/k_cyl*p11 + (k_cyl-i)/k_cyl*p00;
        elements(end+1, :, :) = element;

        element = zeros(4,3);
        element(1, 1:3) = [
            (i+1-k_cyl)/k_cyl*r_square;
            r_square;
            0;
        ]';
        element(2, 1:3) = [
            (i-k_cyl)/k_cyl*r_square;
            r_square;
            0;
        ]';
        element(3, 1:3) = (i)/k_cyl*p00 + (k_cyl-i)/k_cyl*p11;
        element(4, 1:3) = (i+1)/k_cyl*p00 + (k_cyl-i-1)/k_cyl*p11;
        elements(end+1, :, :) = element;

        element = zeros(4,3);
        element(1, 1:3) = [
            (i)/k_cyl*r_square;
            r_square;
            0;
        ]';
        element(2, 1:3) = [
            (i-1)/k_cyl*r_square;
            r_square;
            0;
        ]';
        element(3, 1:3) = (i-1)/k_cyl*p12 + (k_cyl-i+1)/k_cyl*p00;
        element(4, 1:3) = (i)/k_cyl*p12 + (k_cyl-i)/k_cyl*p00;
        elements(end+1, :, :) = element;

        element = zeros(4,3);
        element(1, 1:3) = [
            r_square;
            (k_cyl-i-1)/k_cyl*r_square;
            0;
        ]';
        element(2, 1:3) = [
            r_square;
            (k_cyl-i)/k_cyl*r_square;
            0;
        ]';
        element(3, 1:3) = (i)/k_cyl*p00 + (k_cyl-i)/k_cyl*p12;
        element(4, 1:3) = (i-1)/k_cyl*p00 + (k_cyl-i-1)/k_cyl*p12;
        elements(end+1, :, :) = element;
    end


    for i = 3:size(xs,1)
        r_prev = xs(i-1);
        r = xs(i);
        v = (r-r_square)/(r_cir-r_square);
        v_prev = (r_prev-r_square)/(r_cir-r_square);

        for j = 1:4*k_cyl
            element = zeros(4,3);
            element(1, 1:3) = curve_interp(j,v);
            element(2, 1:3) = curve_interp(j-1,v);
            element(3, 1:3) = curve_interp(j-1,v_prev);
            element(4, 1:3) = curve_interp(j,v_prev);

            elements(end+1, :, :) = element;
            if i == size(xs,1)
                boundaries(end+1, :) = [size(elements, 1); 9;];
            end
        end

    end
    elements(:, :, 2) = elements(:, :, 2)/xs(end)*(xs(end)-R_a);
    elements(:, :, 2) = elements(:, :, 2) + R_a;
    checkCounterClockwise(elements);
    % plotElements(elements, [])
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
                fprintf('Element %d is not counterclockwise at corner %d.\n', k, i);
            end
        end
    end
end


function plotElements(elements, boundaries)
    figure;
    hold on;
    for k = 1:size(elements, 1)
        element = elements(k, :, :);
        element = squeeze(element);
        
        x = element(:, 1);
        y = element(:, 2);
        z = element(:, 3);
        
        x = [x; x(1)];
        y = [y; y(1)];
        z = [z; z(1)];
        % if ismember(k, boundaries(:, 1))
        %     plot(y, z, 'r-', 'LineWidth', 0.5);
        % else
        %     plot(y, z, 'k-', 'LineWidth', 0.5);
        % end
        centroid_y = mean(y(1:end-1));  % exclude the repeated first point
        centroid_x = mean(x(1:end-1));
        
        % Plot element number
        text(centroid_x, centroid_y, num2str(k), 'FontSize', 14, 'Color', 'b');
        plot(x, y, 'k-', 'LineWidth', 0.5);

    end
    axis equal;
    hold off;
end
