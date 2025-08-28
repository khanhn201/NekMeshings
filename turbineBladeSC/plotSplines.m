
figure;
hold on;
grid on;
xlabel('X-axis');  % Axial coordinate (slice index)
ylabel('Y-axis');  % First spatial coordinate
zlabel('Z-axis');  % Second spatial coordinate
title('3D Visualization of Spline Surfaces');

num_points = 1000; % Number of points to sample along each spline

% Plot the slice splines
for i = 1:length(sliceSplines)
    s_vals = linspace(sliceSplines{i}.breaks(1), sliceSplines{i}.breaks(end), num_points);
    slice_points = ppval(sliceSplines{i}, s_vals);
    plot3(slice_points(1, :), slice_points(2, :),slice_points(3, :), 'b', 'LineWidth', 1.5);
end
for i = 1:length(connectingSplines)
    s_vals = linspace(connectingSplines{i}.breaks(1), connectingSplines{i}.breaks(end), num_points);
    slice_points = ppval(connectingSplines{i}, s_vals);
    plot3(slice_points(1, :), slice_points(2, :),slice_points(3, :), 'r', 'LineWidth', 1.5);
end
axis equal;
legend('Slice Splines', 'Connecting Splines');
view(3); % Set to 3D view
hold off;
