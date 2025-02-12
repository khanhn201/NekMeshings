% for i = 52:104:156
for i = 52:104:52
    filename = sprintf('slices/slice%04d.txt', i);
    slice = readSliceFile(filename);
    [pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
    [elementsOuter, boundariesOuters, pp_coarse] = meshOuter(pp, arc_length, arc_length_at_max_y);
    [elementsInner, boundariesInner] = meshInner(pp, arc_length, arc_length_at_max_y);
    elements = [elementsOuter; elementsInner;];
    boundaries = boundariesOuter;
    plotElements(elements, boundaries);

    figure;
    hold on;
    t_values = linspace(pp.breaks(1), pp.breaks(end), 1000);
    spline_points = ppval(pp, t_values);
    plot(slice(:,2), slice(:,3), 'ro-', 'MarkerSize', 3, 'DisplayName', 'Original Points');
    plot(spline_points(2, :), spline_points(3, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline');
    legend;
    axis equal;
    grid on;
    title('Fitted Spline to the Given Slice');
    hold off;

    figure;
    hold on;
    t_values = linspace(pp_coarse.breaks(1), pp_coarse.breaks(end), 1000);
    spline_points = ppval(pp_coarse, t_values);
    plot(slice(:,2), slice(:,3), 'ro-', 'MarkerSize', 3, 'DisplayName', 'Original Points');
    plot(spline_points(2, :), spline_points(3, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline');
    legend;
    axis equal;
    grid on;
    title('Fitted Spline to the Given Slice');
    hold off;
end
