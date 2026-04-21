% elements = (N, 4, 3)
function plotElements2D(elements, boundaries)
    figure;
    hold on;
    
    for k = 1:size(elements,1)
        el = squeeze(elements(k,:,:));
        x = [el(:,1); el(1,1)];
        y = [el(:,2); el(1,2)];
        plot(x, y, 'k-', 'LineWidth', 0.5);
        centroid_x = mean(x(1:end));
        centroid_y = mean(y(1:end));
        
        % Plot element number
        text(centroid_x, centroid_y, num2str(k), 'FontSize', 14, 'Color', 'b');
    end

    % now plot boundaries by face
    for b = 1:size(boundaries,1)
        e  = boundaries(b,1);
        f  = boundaries(b,2);
        bc = boundaries(b,3);

        el = squeeze(elements(e,:,:));

        switch f
            case 1, edge = el([1 2],:);
            case 2, edge = el([2 3],:);
            case 3, edge = el([3 4],:);
            case 4, edge = el([4 1],:);
        end

        plot(edge(:,1), edge(:,2), '-', 'Color', 'r', 'LineWidth', 1);
    end
end
