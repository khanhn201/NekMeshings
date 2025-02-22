for i = 52:52:52
% for i = 156:52:156
% for i = 103:1:105
% for i = 156:104:156
    filename = sprintf('slices/slice%04d.txt', i);
    slice = readSliceFile(filename);
    x = slice(1,1)
    flipped = false;
    if x > 0
        flipped = true;
    end
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice, flipped);

    % figure;
    % hold on;
    % t_values = linspace(pp.breaks(1), pp.breaks(end), 1000);
    % spline_points = ppval(pp, t_values);
    % plot(slice(:,2), slice(:,3), 'ro-', 'MarkerSize', 3, 'DisplayName', 'Original Points');
    % plot(spline_points(2, :), spline_points(3, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline');
    % legend;
    % axis equal;
    % grid on;
    % title('Fitted Spline to the Given Slice1');
    % hold off;

    [elementsOuter, boundariesOuter, pp_coarse] = meshOuterOMesh(pp, arc_length, arc_length_at_max_y, flipped);
    % [elementsOuter, boundariesOuter, pp_coarse] = meshOuter(pp, arc_length, arc_length_at_max_y, flipped);
    [elementsInner, boundariesInner] = meshInner(pp, arc_length, arc_length_at_max_y, flipped);
    elements = [elementsOuter; elementsInner;];
    boundaries = boundariesOuter;
    % plotElements(elements, boundaries);
    plotElements(elements, boundaries);

    figure;
    hold on;
    t_values = linspace(pp_coarse.breaks(1), pp_coarse.breaks(end), 1000);
    spline_points = ppval(pp_coarse, t_values);
    plot(slice(:,2), slice(:,3), 'ro-', 'MarkerSize', 3, 'DisplayName', 'Original Points');
    plot(spline_points(2, :), spline_points(3, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline');
    legend;
    axis equal;
    grid on;
    title('Fitted Spline to the Given Slice2');
    hold off;
end
