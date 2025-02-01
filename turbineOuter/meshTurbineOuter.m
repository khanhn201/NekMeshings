function [elements,boundaries] = meshTurbineOuter()
    config
    elements = [];
    boundaries = [];

    R = R_x;

    % Square
    for i = 1:n_segs
        for j = 1:n_segs
            p1 = [0, (i-1)/n_segs * 2*Ra - Ra, (j-1)/n_segs * 2*Ra - Ra];
            p2 = [0, i/n_segs * 2*Ra - Ra, (j-1)/n_segs * 2*Ra - Ra];
            p3 = [0, i/n_segs * 2*Ra - Ra, j/n_segs * 2*Ra - Ra];
            p4 = [0, (i-1)/n_segs * 2*Ra - Ra, j/n_segs * 2*Ra - Ra];
            element(1,:) = p1;
            element(2,:) = p2;
            element(3,:) = p3;
            element(4,:) = p4;

            elements(end+1, :, :) = element;
            boundaries(end+1, :) = [size(elements, 1); 1;];
        end
    end

    % Inner Circle
    for k = 1:k_inner
        for i = 1:n_segs
            p1 = [0, (i-1)/n_segs * 2*Ra - Ra, Ra];
            p2 = [0, (i)/n_segs * 2*Ra - Ra, Ra];
            p3 = [0, -R*cosd(45 + (i)/n_segs*90), R*sind(45 + (i)/n_segs*90)];
            p4 = [0, -R*cosd(45 + (i-1)/n_segs*90), R*sind(45 + (i-1)/n_segs*90)];
            element(1,:) = (k_inner-k+1)*p1/k_inner + (k-1)*p4/k_inner;
            element(2,:) = (k_inner-k+1)*p2/k_inner + (k-1)*p3/k_inner;
            element(3,:) = (k_inner-k)*p2/k_inner + (k)*p3/k_inner;
            element(4,:) = (k_inner-k)*p1/k_inner + (k)*p4/k_inner;

            elements(end+1, :, :) = element;
            boundaries(end+1, :) = [size(elements, 1); 1;];
        end
    end
    for k = 1:k_inner
        for i = 1:n_segs
            p1 = [0, Ra, (n_segs-i+1)/n_segs * 2*Ra - Ra];
            p2 = [0, Ra, (n_segs-i)/n_segs * 2*Ra - Ra];
            p3 = [0, -R*cosd(135 + (i)/n_segs*90), R*sind(135 + (i)/n_segs*90)];
            p4 = [0, -R*cosd(135 + (i-1)/n_segs*90), R*sind(135 + (i-1)/n_segs*90)];
            element(1,:) = (k_inner-k+1)*p1/k_inner + (k-1)*p4/k_inner;
            element(2,:) = (k_inner-k+1)*p2/k_inner + (k-1)*p3/k_inner;
            element(3,:) = (k_inner-k)*p2/k_inner + (k)*p3/k_inner;
            element(4,:) = (k_inner-k)*p1/k_inner + (k)*p4/k_inner;

            elements(end+1, :, :) = element;
            boundaries(end+1, :) = [size(elements, 1); 1;];
        end
    end
    for k = 1:k_inner
        for i = 1:n_segs
            p1 = [0, (i-1)/n_segs * 2*Ra - Ra, -Ra];
            p2 = [0, (i)/n_segs * 2*Ra - Ra, -Ra];
            p3 = [0, -R*cosd(45 + (i)/n_segs*90), -R*sind(45 + (i)/n_segs*90)];
            p4 = [0, -R*cosd(45 + (i-1)/n_segs*90), -R*sind(45 + (i-1)/n_segs*90)];
            element(1,:) = (k_inner-k+1)*p1/k_inner + (k-1)*p4/k_inner;
            element(2,:) = (k_inner-k)*p1/k_inner + (k)*p4/k_inner;
            element(3,:) = (k_inner-k)*p2/k_inner + (k)*p3/k_inner;
            element(4,:) = (k_inner-k+1)*p2/k_inner + (k-1)*p3/k_inner;

            elements(end+1, :, :) = element;
            boundaries(end+1, :) = [size(elements, 1); 1;];
        end
    end
    for k = 1:k_inner
        for i = 1:n_segs
            p1 = [0, -Ra, (i-1)/n_segs * 2*Ra - Ra];
            p2 = [0, -Ra, (i)/n_segs * 2*Ra - Ra];
            p3 = [0, -R*cosd(-45 + (i)/n_segs*90), R*sind(-45 + (i)/n_segs*90)];
            p4 = [0, -R*cosd(-45 + (i-1)/n_segs*90), R*sind(-45 + (i-1)/n_segs*90)];
            element(1,:) = (k_inner-k+1)*p1/k_inner + (k-1)*p4/k_inner;
            element(2,:) = (k_inner-k+1)*p2/k_inner + (k-1)*p3/k_inner;
            element(3,:) = (k_inner-k)*p2/k_inner + (k)*p3/k_inner;
            element(4,:) = (k_inner-k)*p1/k_inner + (k)*p4/k_inner;

            elements(end+1, :, :) = element;
            boundaries(end+1, :) = [size(elements, 1); 1;];
        end
    end


    % Outer
    for k = 1:k_outer
        for i = 1:n_segs
            p1 = [0, -R*cosd(45 + (i)/n_segs*90), R*sind(45 + (i)/n_segs*90)];
            p2 = [0, -R*cosd(45 + (i-1)/n_segs*90), R*sind(45 + (i-1)/n_segs*90)];
            p3 = [0, (i-1)/n_segs * 2*R_far_horizontal - R_far_horizontal, R_far_top];
            p4 = [0, (i)/n_segs * 2*R_far_horizontal - R_far_horizontal, R_far_top];
            element(1,:) = (k_outer-k+1)*p1/k_outer + (k-1)*p4/k_outer;
            element(2,:) = (k_outer-k)*p1/k_outer + (k)*p4/k_outer;
            element(3,:) = (k_outer-k)*p2/k_outer + (k)*p3/k_outer;
            element(4,:) = (k_outer-k+1)*p2/k_outer + (k-1)*p3/k_outer;

            elements(end+1, :, :) = element;
            if k == 1
                boundaries(end+1, :) = [size(elements, 1); 2;];
            end
            if k == k_outer
                boundaries(end+1, :) = [size(elements, 1); 3;];
            end
        end
    end
    for k = 1:k_outer
        for i = 1:n_segs
            p1 = [0, -R*cosd(135 + (i)/n_segs*90), R*sind(135 + (i)/n_segs*90)];
            p2 = [0, -R*cosd(135 + (i-1)/n_segs*90), R*sind(135 + (i-1)/n_segs*90)];
            p3 = [0, R_far_horizontal, (n_segs-i+1)/n_segs * (R_far_top-R_far_bottom) + R_far_bottom];
            p4 = [0, R_far_horizontal, (n_segs-i)/n_segs * (R_far_top-R_far_bottom) + R_far_bottom];
            element(1,:) = (k_outer-k+1)*p1/k_outer + (k-1)*p4/k_outer;
            element(2,:) = (k_outer-k)*p1/k_outer + (k)*p4/k_outer;
            element(3,:) = (k_outer-k)*p2/k_outer + (k)*p3/k_outer;
            element(4,:) = (k_outer-k+1)*p2/k_outer + (k-1)*p3/k_outer;

            elements(end+1, :, :) = element;
            if k == 1
                boundaries(end+1, :) = [size(elements, 1); 4;];
            end
            if k == k_outer
                boundaries(end+1, :) = [size(elements, 1); 5;];
            end
        end
    end
    for k = 1:k_outer
        for i = 1:n_segs
            p1 = [0, -R*cosd(45 + (i)/n_segs*90), -R*sind(45 + (i)/n_segs*90)];
            p2 = [0, -R*cosd(45 + (i-1)/n_segs*90), -R*sind(45 + (i-1)/n_segs*90)];
            p3 = [0, (i-1)/n_segs * 2*R_far_horizontal - R_far_horizontal, R_far_bottom];
            p4 = [0, (i)/n_segs * 2*R_far_horizontal - R_far_horizontal, R_far_bottom];
            element(1,:) = (k_outer-k+1)*p1/k_outer + (k-1)*p4/k_outer;
            element(2,:) = (k_outer-k+1)*p2/k_outer + (k-1)*p3/k_outer;
            element(3,:) = (k_outer-k)*p2/k_outer + (k)*p3/k_outer;
            element(4,:) = (k_outer-k)*p1/k_outer + (k)*p4/k_outer;

            elements(end+1, :, :) = element;
            if k == 1
                boundaries(end+1, :) = [size(elements, 1); 6;];
            end
            if k == k_outer
                boundaries(end+1, :) = [size(elements, 1); 7;];
            end
        end
    end
    for k = 1:k_outer
        for i = 1:n_segs
            p1 = [0, -R*cosd(-45 + (i)/n_segs*90), R*sind(-45 + (i)/n_segs*90)];
            p2 = [0, -R*cosd(-45 + (i-1)/n_segs*90), R*sind(-45 + (i-1)/n_segs*90)];
            p3 = [0, -R_far_horizontal, (i-1)/n_segs * (R_far_top-R_far_bottom) + R_far_bottom];
            p4 = [0, -R_far_horizontal, (i)/n_segs * (R_far_top-R_far_bottom) + R_far_bottom];
            element(1,:) = (k_outer-k+1)*p1/k_outer + (k-1)*p4/k_outer;
            element(2,:) = (k_outer-k)*p1/k_outer + (k)*p4/k_outer;
            element(3,:) = (k_outer-k)*p2/k_outer + (k)*p3/k_outer;
            element(4,:) = (k_outer-k+1)*p2/k_outer + (k-1)*p3/k_outer;

            elements(end+1, :, :) = element;
            if k == 1
                boundaries(end+1, :) = [size(elements, 1); 8;];
            end
            if k == k_outer
                boundaries(end+1, :) = [size(elements, 1); 9;];
            end
        end
    end

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
            
            if cross_prod(1) <= 0
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
        plot(y, z, 'k-', 'LineWidth', 0.5);

    end
    axis equal;
    hold off;
end

