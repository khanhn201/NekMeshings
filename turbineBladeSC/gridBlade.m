function [elements, boundaries] = gridBlade(zs)
    config
    elements = [];
    boundaries = [];

    r_square = zs(2);
    r_cir = zs(end);
    pd1 = [-R_a; 0; z_shift]';
    pd2 = [R_a; 0; z_shift]';
    pd3 = [-R_a; 0; zs(end)]';
    pd4 = [R_a; 0; zs(end)]';

    k = n_top + 2;
    function point=bottom_curve(u)
        point = pd2*(u)/k + pd1*(k-u)/k;
    end
    function point=top_curve(u)
        point = pd4*(u)/k + pd3*(k-u)/k;
    end
    function point=curve_interp(u, v)
        point = v*top_curve(u) + (1-v)*bottom_curve(u);
    end


    pd5 = [-R_a; 0; R_a/sqrt(3)]';
    pd6 = [R_a; 0; R_a/sqrt(3)]';

    % element = zeros(4,3);
    % element(1, 1:3) = pd2;
    % element(2, 1:3) = bottom_curve(n_top-1);
    % element(3, 1:3) = (n_top-1)/n_top*pd2;
    % element(4, 1:3) = bottom_curve(n_top+1);
    % elements(end+1, :, :) = element;
    for j = 1:k
        element = zeros(4,3);
        element(1, 1:3) = bottom_curve(j);
        element(2, 1:3) = bottom_curve(j-1);
        element(3, 1:3) = pd5*(k/2-j+1)/k*2;
        element(4, 1:3) = pd5*(k/2-j)/k*2;
        if j > k/2
            element(3, 1:3) = pd6*(j-k/2-1)/k*2;
            element(4, 1:3) = pd6*(j-k/2)/k*2;
        end

        elements(end+1, :, :) = element;
        if i == size(zs,1)
            boundaries(end+1, :) = [size(elements, 1); 9;];
        end
    end

    for i = 3:size(zs,1)
        r_prev = zs(i-1);
        r = zs(i);
        v = (r-r_square)/(r_cir-r_square);
        v_prev = (r_prev-r_square)/(r_cir-r_square);

        for j = 1:k
            element = zeros(4,3);
            element(1, 1:3) = curve_interp(j,v);
            element(2, 1:3) = curve_interp(j-1,v);
            element(3, 1:3) = curve_interp(j-1,v_prev);
            element(4, 1:3) = curve_interp(j,v_prev);

            elements(end+1, :, :) = element;
            if i == size(zs,1)
                boundaries(end+1, :) = [size(elements, 1); 9;];
            end
        end

    end
    checkCounterClockwise(elements);

    % plotElements(elements, []);
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
            
            if cross_prod(2) >= 0
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
        centroid_y = mean(z(1:end-1));  % exclude the repeated first point
        centroid_x = mean(x(1:end-1));
        %
        % % Plot element number
        text(centroid_x, centroid_y, num2str(k), 'FontSize', 14, 'Color', 'b');
        plot(x, z, 'k-', 'LineWidth', 0.5);

    end
    axis equal;
    hold off;
end
