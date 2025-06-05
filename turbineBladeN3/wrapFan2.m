function [elements, boundaries] = wrapFanDiamond(zs)
    config
    elements = [];
    boundaries = [];

    displacement = R_a/sqrt(3);
    zs(:) = zs(:) - displacement;

    k = size(zs,1);
    kh = k_cyl;

    r_square = zs(kh);
    r_max = zs(end);

    pd1 = [-r_square*cos(- pi/6); 0; r_square*sin(- pi/6)]';
    pd2 = [-r_square*cos(pi/3 - pi/6); 0; r_square*sin(pi/3 - pi/6)]';
    pd3 = [-r_square*cos(2*pi/3 - pi/6); 0; r_square*sin(2*pi/3 - pi/6)]';


    km = 0:2*kh-2;  % 2*kh-1 points
    theta = cos(pi * km / (2*kh-2));
    theta = 1-(theta + 1)/2;
    fv = theta;
    theta = (-pi/6) + (2*pi/3) * theta;


    u1 = [-zs(kh:end)' * cos(-pi/6); 
           zeros(1, k-kh+1); 
          zs(kh:end)' * sin(-pi/6)]';
    u2 = [-zs(kh:end)' * cos(pi/2); 
           zeros(1, k-kh+1); 
          zs(kh:end)' * sin(pi/2)]';
    v1 = [linspace(pd1, pd2, kh)'(1:end-1,:); linspace(pd2, pd3, kh)'];
    v2 = [-r_max * cos(theta); 
           zeros(1, 2*kh-1); 
          r_max * sin(theta)]';
    [elements1, boundaries1] = meshQuad(u1, u2, v1, v2, 3, 0, fv);

    for elem=1:size(elements1,1)
        elements(end+1, :, :) = elements1(elem, :, :);
        [isBoundary, idx] = ismember(elem, boundaries1(:, 1));
        if isBoundary
            boundaries(end+1, :) = [size(elements, 1); 9;];
        end
    end

    zss = [0.0; zs(2:kh)];

    u1 = [-zss' * cos(-pi/6); 
           zeros(1, kh); 
          zss' * sin(-pi/6)]';
    v1 = [-zss' * cos(pi/2); 
           zeros(1, kh); 
          zss' * sin(pi/2)]';
    v2 = linspace(pd1, pd2, kh)';
    u2 = linspace(pd3, pd2, kh)';
    [elements1, boundaries1] = meshQuad(u1, u2, v1, v2, 3, 0, 0);
    for elem=1:size(elements1,1)
        elements(end+1, :, :) = elements1(elem, :, :);
        [isBoundary, idx] = ismember(elem, boundaries1(:, 1));
    end


    checkCounterClockwise(elements);
    % plotElements(elements, [], []);

    % Displace
    elements(:, :, 1) -= R_a;
    elements(:, :, 3) += displacement;


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


function plotElements(elements, boundaries, points)
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
    size(points)
    plot(points(:, 1), points(:, 3), 'k-', 'LineWidth', 0.5);
    axis equal;
    hold off;
end
