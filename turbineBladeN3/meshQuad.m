function [elements, boundaries] = meshQuad(u1, u2, v1, v2, bc, fu, fv)
    nu = size(u1, 1);
    nv = size(v1, 1);

    elements = [];
    boundaries = [];
    if fu==0
        fu = linspace(0, 1, nu);
    end
    if fv==0
        fv = linspace(0, 1, nv);
    end
    
    function point = interp(iu,jv)
        v = fv(jv);
        u = fu(iu);
        point = (1-v)*u1(iu, :) + v*u2(iu, :) +...
                (1-u)*v1(jv, :) + u*v2(jv, :) -...
                (1-v)*(1-u)*u1(1, :) -...
                (1-v)*u*v2(1, :) -...
                v*(1-u)*u2(1, :) -...
                u*v*v2(end, :);
        % point = u2(end, :) +...
        %         v2(end, :) -...
        %         v1(end, :);
    end


    for i=2:nu
    for j=2:nv
        element = zeros(4,3);
        element(1, 1:3) = interp(i,j);
        element(2, 1:3) = interp(i,j-1);
        element(3, 1:3) = interp(i-1,j-1);
        element(4, 1:3) = interp(i-1,j);
        elements(end+1, :, :) = element;
        if i == nu && bc == 3
            boundaries(end+1, :) = [size(elements, 1); 9;];
        end
    end
    end
    % plotElements(elements, []);
end
% function plotElements(elements, boundaries)
%     figure;
%     hold on;
%     for k = 1:size(elements, 1)
%         element = elements(k, :, :);
%         element = squeeze(element);
%
%         x = element(:, 1);
%         y = element(:, 2);
%
%         x = [x; x(1)];
%         y = [y; y(1)];
%         % if ismember(k, boundaries(:, 1))
%         %     plot(y, z, 'r-', 'LineWidth', 0.5);
%         % else
%         %     plot(y, z, 'k-', 'LineWidth', 0.5);
%         % end
%         plot(x, y, 'k-', 'LineWidth', 0.5);
%
%         centroid_x = mean(x(1:end-1));  % exclude the repeated first point
%         centroid_y = mean(y(1:end-1));
%
%         % Plot element number
%         % text(centroid_x, centroid_y, num2str(k), 'FontSize', 14, 'Color', 'b');
%
%     end
%     axis equal;
%     hold off;
% end
