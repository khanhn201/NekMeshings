
function plotElementsSym(elements, boundaries)
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
        plot(-x, y, 'b-', 'LineWidth', 0.5);
        
        centroid_x = mean(x(1:end-1));  % exclude the repeated first point
        centroid_y = mean(y(1:end-1));
        
        % Plot element number
        % text(centroid_x, centroid_y, num2str(k), 'FontSize', 14, 'Color', 'b');

    end
    axis equal;
    hold off;
end
