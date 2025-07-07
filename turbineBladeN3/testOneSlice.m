slicesCoord = readSlices('nrel5mw2.mat');
slice = squeeze(slicesCoord(end, :, :));
[pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
[elements, boundaries, pp_coarse] = meshOuterElliptic(pp, arc_length, arc_length_at_max_y);
[elementsInner, boundariesInner] = meshInnerRec(pp, arc_length);
elements = [elements; elementsInner;];
% plotElements(elements, []);
[elements, boundaries, pp_coarse] = meshHub();
plotElementsSym(elements, []);

[elements, boundaries] = wrapFan2([0, 4, 6, 8, 10, 12, 14, 16, 18, 20]');


% figure;
% hold on;
% t_values = linspace(pp_coarse.breaks(1), pp_coarse.breaks(end), 1000);
% spline_points = ppval(pp, t_values);
% plot(slice(:,1), slice(:,2), 'ro-', 'MarkerSize', 3, 'DisplayName', 'Original Points');
% plot(spline_points(1, :), spline_points(2, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Spline');
% legend;
% axis equal;
% grid on;
% title('Fitted Spline to the Given Slice1');
% hold off;
