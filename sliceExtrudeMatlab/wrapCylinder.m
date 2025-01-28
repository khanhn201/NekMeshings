function [elements, boundaries] = wrapCylinder(xs)
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


