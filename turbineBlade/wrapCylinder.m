function [elements, boundaries] = wrapCylinder(xs)
    config
    elements = [];
    boundaries = [];

    % Inner 
    r = xs(2);
    rt =(r + xs(3))/2;
    p0 = [0;0;0;]';
    p1 = [-rt/sqrt(2); rt/sqrt(2); 0;]';
    p2 = [-r; 0; 0;]';
    p3 = [0; r; 0;]';
    p4 = (p1+p2+p3)/4;

    for i=1:2
        element = [];
        element(1, :) = p0;
        element(2, :) = p3;
        element(3, :) = (p3+p1)/2;
        element(4, :) = p4;

        theta = -90 *(i-1)* pi / 180;
        R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 0];
        element = (R * element')';
        elements(end+1, :, :) = element;

        element = [];
        element(1, :) = p0;
        element(2, :) = p4;
        element(3, :) = (p2+p1)/2;
        element(4, :) = p2;

        theta = -90 *(i-1)* pi / 180;
        R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 0];
        element = (R * element')';
        elements(end+1, :, :) = element;

        element = [];
        element(1, :) = p4;
        element(2, :) = (p1+p3)/2;
        element(3, :) = p1;
        element(4, :) = (p1+p2)/2;

        theta = -90 *(i-1)* pi / 180;
        R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 0];
        element = (R * element')';
        elements(end+1, :, :) = element;
    end

    for i = 1:2
        r_prev = xs(2);
        r = xs(3);

        for j = 1:4
            element = [];
            rm = r_prev;
            if j == 1
                element(1, 1:3) = p2;
                element(2, 1:3) = (p1 + p2)/2;
            end
            if j == 2
                element(1, 1:3) = (p1 + p2)/2;
                element(2, 1:3) = p1;
            end
            if j == 3
                element(1, 1:3) = p1;
                element(2, 1:3) = (p1 + p3)/2;
            end
            if j == 4
                element(1, 1:3) = (p1 + p3)/2;
                element(2, 1:3) = p3;
            end
            element(4, 1) = -r*cosd((j-1)*22.5);
            element(4, 2) = r*sind((j-1)*22.5);
            element(3, 1) = -r*cosd((j)*22.5);
            element(3, 2) = r*sind((j)*22.5);
            element(1:4,3) = 0; % z coord

            theta = -90 * (i-1)*pi / 180;
            R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 0];
            element = (R * element')';
            elements(end+1, :, :) = element;
        end

    end

    for i = 4:size(xs,1)
        r_prev = xs(i-1);
        r = xs(i);

        for j = 1:8
            element = [];
            element(1, 1) = -r_prev;
            element(1, 2) = 0;
            element(2, 1) = -r_prev*cosd(22.5);
            element(2, 2) = r_prev*sind(22.5);
            element(4, 1) = -r;
            element(4, 2) = 0;
            element(3, 1) = -r*cosd(22.5);
            element(3, 2) = r*sind(22.5);
            element(1:4,3) = 0; % z coord
            theta = -22.5 * (j-1)*pi / 180;
            R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 0];
            element = (R * element')';

            elements(end+1, :, :) = element;
            if i == size(xs,1)
                boundaries(end+1, :) = [size(elements, 1); 9;];
            end
        end

    end
    elements(:, :, 2) = elements(:, :, 2) + R_a;
    checkCounterClockwise(elements);
end


function [elements, boundaries] = wrapCylinder4Angles(xs)
    config
    elements = [];
    boundaries = [];

    % Inner 
    r = xs(2);
    element = [];
    element(1, 1) = 0;
    element(1, 2) = 0;
    element(2, 1) = r;
    element(2, 2) = 0;
    element(3, 1) = r/sqrt(2);
    element(3, 2) = r/sqrt(2);
    element(4, 1) = 0;
    element(4, 2) = r;
    element(1:4,3) = 0; % z coord
    elements(end+1, :, :) = element;
    element = [];
    element(1, 1) = 0;
    element(1, 2) = 0;
    element(2, 1) = 0;
    element(2, 2) = r;
    element(3, 1) = -r/sqrt(2);
    element(3, 2) = r/sqrt(2);
    element(4, 1) = -r;
    element(4, 2) = 0;
    element(1:4,3) = 0; % z coord
    elements(end+1, :, :) = element;

    for i = 3:size(xs,1)
        r_prev = xs(i-1);
        r = xs(i);

        for j = 1:4
            element = [];
            element(1, 1) = -r_prev;
            element(1, 2) = 0;
            element(2, 1) = -r_prev/sqrt(2);
            element(2, 2) = r_prev/sqrt(2);
            element(4, 1) = -r;
            element(4, 2) = 0;
            element(3, 1) = -r/sqrt(2);
            element(3, 2) = r/sqrt(2);
            element(1:4,3) = 0; % z coord
            theta = -45 * (j-1)*pi / 180;
            R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 0];
            element = (R * element')';

            elements(end+1, :, :) = element;
            if i == size(xs,1)
                boundaries(end+1, :) = [size(elements, 1); 9;];
            end
        end

    end
    elements(:, :, 2) = elements(:, :, 2) + R_a;
    checkCounterClockwise(elements);
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
        plot(x, y, 'k-', 'LineWidth', 0.5);

    end
    axis equal;
    hold off;
end
