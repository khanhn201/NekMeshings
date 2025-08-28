slicesCoord = readSlices('nrel5mw2.mat');
config
zss = slicesCoord(:, 1, 3);
zs = linspace(z_shift, zss(end), n_slices);
slicesCoord = interp1(zss, slicesCoord, zs);

slice = squeeze(slicesCoord(5, :, :));
% slice = squeeze(slicesCoord(1, :, :));
[pp, arc_length, arc_length_at_max_y] = fitSpline(slice);
[elements, boundaries, pp_coarse] = meshOuterElliptic(pp, arc_length, arc_length_at_max_y);
[elementsInner, boundariesInner] = meshInnerRec(pp, arc_length, arc_length_at_max_y);
elements = [elements; elementsInner;];
plotElements(elements, []);
config;
da = R_a/sqrt(3);
projAngle = atan2(da,R_a);
i=5;
da = (i-1)/(hub_layers)*(z_shift-R_a/sqrt(3)) + R_a/sqrt(3);
db = (i-1)/(hub_layers)*z_shift;
    projAngle = atan2(da - db,R_a);
[elements, boundaries, pp_coarse] = meshHub(projAngle, db);
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
